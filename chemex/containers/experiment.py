import contextlib as cl
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

    def add(self, experiment):
        self._experiments[experiment.filename] = experiment

    def residuals(self, params, verbose=True, threshold=1e-3):
        residuals = np.asarray(
            list(
                it.chain.from_iterable(
                    experiment.residuals(params)
                    for experiment in self._experiments.values()
                )
            )
        )
        if verbose:
            chisq = sum(np.asarray(residuals) ** 2)
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
            relevant_subset.add(experiment.get_relevant_subset(par_names))
        return relevant_subset

    def write(self, path, params):
        for experiment in self._experiments.values():
            experiment.write(path, params)

    def plot(self, path, params):
        for experiment in self._experiments.values():
            experiment.plot(path, params)

    def select(self, included=None, excluded=None):
        for experiment in self._experiments.values():
            experiment.select(included, excluded)

    def filter(self, params=None):
        for experiment in self._experiments.values():
            experiment.filter(params)

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
    def params_default(self):
        return cph.merge(
            profile.params_default for profile in self._experiments.values()
        )

    def __len__(self):
        return sum([len(experiment) for experiment in self._experiments.values()])

    def __bool__(self):
        return bool(len(self))


class RelaxationExperiment:
    def __init__(self, filename, exp_type, profiles, verbose=True):
        self.filename = filename
        self.exp_type = exp_type
        self._profiles = profiles
        self._filtered = {}
        if verbose:
            print(f"  - Experiment: {exp_type}")
            print(f"  - Profiles: {len(profiles)}")

    def residuals(self, params):
        return list(
            it.chain.from_iterable(
                profile.residuals(params) for profile in self._profiles.values()
            )
        )

    @property
    def par_name_sets(self):
        return list(
            set(profile.params_default.keys()) for profile in self._profiles.values()
        )

    def get_relevant_subset(self, par_names):
        profiles = {}
        for filename, profile in self._profiles.items():
            if set(profile.params_default) & set(par_names):
                profiles[filename] = profile
        return RelaxationExperiment(
            self.filename, self.exp_type, profiles=profiles, verbose=False
        )

    def plot(self, path, params):
        basename = path / self.filename.name
        name_pdf = basename.with_suffix(".pdf")
        name_exp = basename.with_suffix(".exp")
        name_fit = basename.with_suffix(".fit")
        print(f"  - {name_pdf} [.fit, .exp]")
        with cl.ExitStack() as stack:
            file_pdf = stack.enter_context(pdf.PdfPages(str(name_pdf)))
            file_exp = stack.enter_context(name_exp.open("w"))
            file_fit = stack.enter_context(name_fit.open("w"))
            for profile in self._profiles.values():
                profile.plot(file_pdf, params)
                profile.write_plot(file_exp, file_fit, params)

    def write(self, path, params):
        filename = (path / self.filename.name).with_suffix(".dat")
        print(f"  - {filename}")
        with filename.open("w") as file_dat:
            profiles = self._profiles.values()
            profiles_sorted = sorted(profiles, key=lambda profile_: profile_.name)
            for profile in profiles_sorted:
                file_dat.write(profile.print(params=params))

    def select(self, included=None, excluded=None):
        # TODO: Clean this
        if included is None:
            included = list(str(profile.name) for profile in self._profiles.values())
        if excluded is None:
            excluded = []
        included = [name.lower() for name in included]
        excluded = [name.lower() for name in excluded]
        for name, profile in self._profiles.copy().items():
            if str(profile.name) not in included or str(profile.name) in excluded:
                self._filtered[name] = self._profiles.pop(name)

    def filter(self, params=None):
        for profile in self._profiles.values():
            profile.filter(params)

    def monte_carlo(self, params):
        profiles = {}
        for name, profile in self._profiles.items():
            profiles[name] = profile.monte_carlo(params)
        return RelaxationExperiment(
            self.filename, self.exp_type, profiles, verbose=False
        )

    def bootstrap(self):
        profiles = {}
        for name, profile in self._profiles.items():
            profiles[name] = profile.bootstrap()
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
        for profile in self._profiles.values():
            noise_variance_values.append(profile.estimate_noise_variance(kind))
        noise_mean = np.sqrt(np.mean(noise_variance_values))
        for profile in self._profiles.values():
            profile.set_noise(noise_mean)

    @property
    def params_default(self):
        return cph.merge(profile.params_default for profile in self._profiles.values())

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
