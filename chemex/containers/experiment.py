import abc
import contextlib as cl
import functools as ft
import itertools as it

import lmfit as lf
import matplotlib.backends.backend_pdf as pdf
import numpy as np

import chemex.containers.plot as ccp
import chemex.experiments as ce
import chemex.helper as ch
import chemex.nmr.rates as cnr
import chemex.parameters.helper as cph


class Experiments:
    def __init__(self):
        self._experiments = {}
        self._chisq_ref = 1e32
        self.verbose = False

    def add(self, experiment: "Experiment"):
        self._experiments[experiment.filename] = experiment

    def residuals(self, params, threshold=1e-3):
        residuals = np.asarray(
            list(
                it.chain.from_iterable(
                    experiment.residuals(params)
                    for experiment in self._experiments.values()
                )
            )
        )
        if self.verbose:
            chisq = (residuals ** 2).sum()
            change = (chisq - self._chisq_ref) / self._chisq_ref
            if change < -threshold:
                nvarys = len([param for param in params.values() if param.vary])
                redchi = chisq / (len(residuals) - nvarys)
                print(f"  - {chisq:.3e} / {redchi:.3e}")
                self._chisq_ref = chisq
        return residuals

    def write(self, params, path):
        path_ = path / "Data"
        path_.mkdir(parents=True, exist_ok=True)
        for experiment in self._experiments.values():
            experiment.write(params, path_)

    def plot(self, params, path, simulation=False):
        for experiment in self._experiments.values():
            experiment.plot(params, path, simulation)

    def select(self, selection=None):
        if selection is None:
            selection = {}
        include, exclude = (selection.get(key) for key in ("include", "exclude"))
        if include is None and exclude is None:
            return
        print("\nSelecting profiles...")
        for experiment in self._experiments.values():
            experiment.select(include, exclude)
        print(f"  - Profile(s): {len(self)}")

    def filter(self, params=None):
        for experiment in self._experiments.values():
            experiment.filter(params)

    def merge_same_profiles(self):
        for experiment in self._experiments.values():
            experiment.merge_same_profiles()

    def monte_carlo(self, params):
        experiments_mc = Experiments()
        for experiment in self._experiments.values():
            experiments_mc.add(experiment.monte_carlo(params))
        return experiments_mc

    def bootstrap(self):
        experiments_bs = Experiments()
        for experiment in self._experiments.values():
            experiments_bs.add(experiment.bootstrap())
        return experiments_bs

    @property
    def params(self):
        return cph.merge(experiment.params for experiment in self._experiments.values())

    def select_params(self, params):
        pnames = self.params.keys()
        selected = lf.Parameters(usersyms=cnr.rate_functions)
        selected.add_many(*(params[pname] for pname in pnames))
        return selected

    def get_relevant_subset(self, pnames):
        relevant_subset = Experiments()
        for experiment in self._experiments.values():
            subset = experiment.get_relevant_subset(pnames)
            if subset:
                relevant_subset.add(subset)
        return relevant_subset

    @property
    def pname_sets(self):
        return list(
            it.chain.from_iterable(
                experiment.pname_sets for experiment in self._experiments.values()
            )
        )

    def get_cluster_name(self):
        names = [
            experiment.get_cluster_name() for experiment in self._experiments.values()
        ]
        return ft.reduce(lambda a, b: a & b, names)

    def __len__(self):
        return sum([len(experiment) for experiment in self._experiments.values()])

    def __bool__(self):
        return bool(len(self))

    def __iter__(self):
        yield from self._experiments.values()


class Experiment(abc.ABC):
    print_name = ""

    def __init__(self, config, profiles, verbose=True):
        self.config = config
        self.filename = config["filename"]
        self.name = config["experiment"]["name"]
        self._profiles = sorted(profiles)
        self._filtered = []
        if verbose:
            print(f"  - Experiment: {self.name}")
            print(f"  - {self.print_name}: {len(profiles)}")

    def residuals(self, params):
        return list(
            it.chain.from_iterable(
                profile.residuals(params) for profile in self._profiles
            )
        )

    @abc.abstractmethod
    def plot(self, params, path, simulation=False):
        pass

    def write(self, params, path):
        filename = (path / self.filename.name).with_suffix(".dat")
        with filename.open("w") as file_dat:
            for profile in sorted(self._profiles):
                file_dat.write(profile.print(params))

    def select(self, include=None, exclude=None):
        profiles_all = self._profiles + self._filtered
        if isinstance(include, str) and include.lower() in ("all", "*"):
            self._profiles = profiles_all
            self._filtered = []
            return
        profiles, filtered = [], []
        for profile in profiles_all:
            included = include is None or profile.name.part_of(include)
            excluded = exclude is not None and profile.name.part_of(exclude)
            if included and not excluded:
                profiles.append(profile)
            else:
                filtered.append(profile)
        self._profiles = profiles
        self._filtered = filtered

    def filter(self, params=None):
        for profile in self._profiles:
            profile.filter(params)

    def merge_same_profiles(self):
        self._profiles = _merge_same_profiles(self._profiles)
        self._filtered = _merge_same_profiles(self._filtered)

    def monte_carlo(self, params):
        profiles = [profile.monte_carlo(params) for profile in self._profiles]
        return type(self)(self.config, profiles, verbose=False)

    @abc.abstractmethod
    def bootstrap(self):
        pass

    def estimate_noise(self, kind):
        implemented = ("file", "scatter", "duplicates")
        if kind not in implemented:
            print(
                f"Warning: Experiment {self.filename.name}: The method '{kind}' is not "
                f"implemented. Please choose one of the following methods: "
                f"{implemented}"
            )
            kind = "file"
        if kind == "duplicates" and not self._any_duplicate():
            print(
                f"Warning: Experiment {self.filename.name}: Some profiles have no "
                f"duplicate points: Uncertainties are not estimated and directly taken "
                f"from the files."
            )
            kind = "file"
        if kind == "file" or not self._profiles:
            return
        noise_variance_values = [
            profile.estimate_noise_variance(kind) for profile in self._profiles
        ]
        noise_mean = np.sqrt(np.mean(noise_variance_values))
        for profile in self._profiles:
            profile.set_noise(noise_mean)

    @property
    def params(self):
        return cph.merge(profile.params for profile in self._profiles)

    @property
    def pname_sets(self):
        return list(set(profile.params) for profile in self._profiles)

    def get_relevant_subset(self, pnames):
        profiles = [
            profile for profile in self._profiles if set(profile.params) & set(pnames)
        ]
        return type(self)(self.config, profiles=profiles, verbose=False)

    def get_cluster_name(self):
        names = [profile.name for profile in self._profiles]
        return ft.reduce(lambda a, b: a & b, names)

    def __len__(self):
        return len(self._profiles)

    def __bool__(self):
        return bool(len(self._profiles))

    def __iter__(self):
        yield from self._profiles

    def _any_duplicate(self):
        return any(profile.any_duplicate() for profile in self._profiles)


class RelaxationExperiment(Experiment):
    print_name = "Profiles"

    def plot(self, params, path, simulation=False):
        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")
        name_exp = basename.with_suffix(".exp")
        name_fit = basename.with_suffix(".fit")
        print(f"  - {name_pdf} [.fit, .exp]")
        with cl.ExitStack() as stack:
            file_pdf = stack.enter_context(pdf.PdfPages(str(name_pdf)))
            file_fit = stack.enter_context(name_fit.open("w"))
            file_exp = None if simulation else stack.enter_context(name_exp.open("w"))
            for profile in sorted(self._profiles):
                profile.plot(params, file_pdf, file_exp, file_fit, simulation)

    def bootstrap(self):
        profiles = [profile.bootstrap() for profile in self._profiles]
        return type(self)(self.config, profiles, verbose=False)


class ShiftExperiment(Experiment):
    print_name = "Shifts"

    def plot(self, params, path, simulation=False):
        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")
        print(f"  - {name_pdf}")
        fit, exp, err = [], [], []
        for profile in self._profiles:
            fit.append(profile.calculate(params))
            exp.append(profile.data["shift"])
            err.append(profile.data["error"])
        ccp.shift(name_pdf, self.filename.name, fit, exp, err)

    def bootstrap(self):
        profiles = np.random.choice(self._profiles, len(self._profiles))
        return type(self)(self.config, profiles, verbose=False)


def read(filenames=None, model=None, selection=None, defaults=None):
    if not filenames:
        return None
    ch.header1("Reading experimental data")
    experiments = Experiments()
    for filename in filenames:
        print(f"\nReading '{filename}'...")
        experiment = ce.read(filename, model, selection, defaults)
        experiments.add(experiment)
    return experiments


def _merge_same_profiles(profiles):
    merged = []
    profile_sets = {}
    for profile in profiles:
        profile_sets.setdefault(profile.name, []).append(profile)
    for profile_set in profile_sets.values():
        if len(profile_set) > 0:
            merged.append(np.sum(profile_set))
        else:
            merged.append(profile_set.pop())
    return merged
