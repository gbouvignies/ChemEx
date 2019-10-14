import contextlib as cl
import functools as ft
import itertools as it

import matplotlib.backends.backend_pdf as pdf
import numpy as np

import chemex.experiments as ce
import chemex.helper as ch
import chemex.parameters.helper as cph


class Experiments:
    def __init__(self):
        self._experiments = {}
        self._chisq_ref = 1e32
        self.verbose = False

    def add(self, experiment):
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

    @property
    def par_name_sets(self):
        return list(
            it.chain.from_iterable(
                experiment.par_name_sets for experiment in self._experiments.values()
            )
        )

    def get_relevant_subset(self, par_names):
        relevant_subset = Experiments()
        for experiment in self._experiments.values():
            subset = experiment.get_relevant_subset(par_names)
            if subset:
                relevant_subset.add(subset)
        return relevant_subset

    def write(self, params, path):
        path_ = path / "Data"
        path_.mkdir(parents=True, exist_ok=True)
        for experiment in self._experiments.values():
            experiment.write(params, path_)

    def plot(self, params, path):
        for experiment in self._experiments.values():
            experiment.plot(params, path)

    def select(self, selection=None, discard=False):
        for experiment in self._experiments.values():
            experiment.select(selection, discard)

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

    def set_params(self, params, model_free):
        for experiment in self._experiments.values():
            experiment.set_params(params, model_free)

    @property
    def params_default(self):
        return cph.merge(
            profile.params_default for profile in self._experiments.values()
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


class RelaxationExperiment:
    def __init__(self, filename, exp_type, profiles, rates=None, verbose=True):
        self.filename = filename
        self.exp_type = exp_type
        self._profiles = sorted(profiles)
        self._filtered = []
        self._rates = rates
        if verbose:
            print(f"  - Experiment: {exp_type}")
            print(f"  - Profiles: {len(profiles)}")

    def residuals(self, params):
        return list(
            it.chain.from_iterable(
                profile.residuals(params) for profile in self._profiles
            )
        )

    @property
    def par_name_sets(self):
        return list(set(profile.params_default.keys()) for profile in self._profiles)

    def get_relevant_subset(self, par_names):
        profiles = []
        for profile in self._profiles:
            if set(profile.params_default) & set(par_names):
                profiles.append(profile)
        return RelaxationExperiment(
            self.filename, self.exp_type, profiles=profiles, verbose=False
        )

    def plot(self, params, path):
        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")
        name_exp = basename.with_suffix(".exp")
        name_fit = basename.with_suffix(".fit")
        print(f"  - {name_pdf} [.fit, .exp]")
        with cl.ExitStack() as stack:
            file_pdf = stack.enter_context(pdf.PdfPages(str(name_pdf)))
            file_exp = stack.enter_context(name_exp.open("w"))
            file_fit = stack.enter_context(name_fit.open("w"))
            for profile in sorted(self._profiles):
                profile.plot(params, file_pdf)
                profile.write_plot(params, file_exp, file_fit)

    def write(self, params, path):
        filename = (path / self.filename.name).with_suffix(".dat")
        with filename.open("w") as file_dat:
            for profile in sorted(self._profiles):
                file_dat.write(profile.print(params))

    def select(self, selection=None, discard=False):
        if selection is None:
            return
        profiles = self._profiles + self._filtered
        if isinstance(selection, str) and selection.lower() in ("all", "*"):
            selected = profiles
        else:
            selected = []
            for profile in profiles:
                for name_incl in selection:
                    if profile.name & name_incl == name_incl:
                        selected.append(profile)
                        break
        self._profiles = sorted(selected)
        if not discard:
            self._filtered = sorted(
                profile for profile in profiles if profile not in selected
            )
        else:
            self._filtered = []

    def filter(self, params=None):
        for profile in self._profiles:
            profile.filter(params)

    def merge_same_profiles(self):
        self._profiles = _merge_same_profiles(self._profiles)
        self._filtered = _merge_same_profiles(self._filtered)

    def monte_carlo(self, params):
        profiles = []
        for profile in self._profiles:
            profiles.append(profile.monte_carlo(params))
        return RelaxationExperiment(
            self.filename, self.exp_type, profiles, verbose=False
        )

    def bootstrap(self):
        profiles = []
        for profile in self._profiles:
            profiles.append(profile.bootstrap())
        return RelaxationExperiment(
            self.filename, self.exp_type, profiles, verbose=False
        )

    def estimate_noise(self, kind):
        implemented = ("file", "scatter", "duplicates")
        if kind not in implemented:
            print(
                f"Warning: Experiment {self.filename.name}: The method '{kind}' is not "
                f"implemented. Please choose one of the following methods: "
                f"{implemented}"
            )
            kind = "file"
        if kind == "file":
            return
        noise_variance_values = []
        for profile in self._profiles:
            noise_variance_values.append(profile.estimate_noise_variance(kind))
        noise_mean = np.sqrt(np.mean(noise_variance_values))
        for profile in self._profiles:
            profile.set_noise(noise_mean)

    def set_params(self, params, model_free):
        rates = {}
        if self._rates is not None:
            rates = self._rates.calculate(model_free)
        for profile in self._profiles:
            profile.set_params(params, rates)

    @property
    def params_default(self):
        return cph.merge(profile.params_default for profile in self._profiles)

    def get_cluster_name(self):
        names = [profile.name for profile in self._profiles]
        return ft.reduce(lambda a, b: a & b, names)

    def __len__(self):
        return len(self._profiles)


def read(filenames=None, model=None):
    if not filenames:
        return None
    ch.header1("Reading experimental data")
    experiments = Experiments()
    for filename in filenames:
        print(f"\nReading '{filename}'...")
        experiment = ce.read(filename, model)
        experiments.add(experiment)
    return experiments


def _merge_same_profiles(profiles):
    merged = []
    profile_sets = {}
    for profile in profiles:
        profile_sets.setdefault(profile.name, []).append((profile))
    for profile_set in profile_sets.values():
        if len(profile_set) > 0:
            merged.append(np.sum(profile_set))
        else:
            merged.append(profile_set.pop())
    return merged
