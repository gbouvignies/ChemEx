import copy
import functools as ft
import sys

import numpy as np

import chemex.containers.helper as cch
import chemex.containers.noise as ccn
import chemex.parameters.name as cpn
import chemex.plot as cp


CPMG_SCHEMA = {
    "type": "object",
    "properties": {
        "data": {
            "type": "object",
            "properties": {
                "error": {
                    "type": "string",
                    "enum": ["file", "duplicates"],
                    "default": "file",
                },
                "filter_planes": {
                    "type": "array",
                    "items": {"type": "integer"},
                    "default": [],
                },
                "path": {"type": "string", "default": "./"},
                "profiles": {
                    "type": "array",
                    "items": {
                        "type": "array",
                        "minItems": 2,
                        "maxItems": 2,
                        "items": {"type": "string"},
                    },
                },
            },
            "required": ["profiles"],
        }
    },
}


@ft.total_ordering
class CpmgProfile:
    def __init__(self, name, data, pulse_seq, par_names, params_default):
        self.name = name
        self.data = data
        self._pulse_seq = pulse_seq
        self._par_names = par_names
        self.params_default = params_default
        self._plot = cp.cpmg

    @classmethod
    def from_file(cls, path, config, pulse_seq, par_names, params_default):
        name = config["spin_system"]["spin_system"]
        data = CpmgData.from_file(
            path,
            filter_planes=config["data"]["filter_planes"],
            time_t2=config["experiment"]["time_t2"],
        )
        return cls(name, data, pulse_seq, par_names, params_default)

    def residuals(self, params):
        data = self.data.points[self.data.mask]
        residuals = (self.calculate(params) - data["intensities"]) / data["errors"]
        return residuals

    def calculate(self, params, ncycs=None):
        data = self.data.points[self.data.mask]
        par_values = self._get_parvals(params)
        calculated = self._pulse_seq.calculate(tuple(data["ncycs"]), par_values)
        scale = cch.get_scale(data["intensities"], data["errors"], calculated)
        if ncycs is not None:
            calculated = self._pulse_seq.calculate(tuple(ncycs), par_values)
        return scale * calculated

    def estimate_noise_variance(self, kind):
        return self.data.estimate_noise_variance(kind)

    def set_noise(self, value):
        self.data.points["errors"] = value

    def print(self, params):
        output = f"[{self.name}]\n"
        output += f"# {'NCYC':>12s}  {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s} {'INTENSITY (CALC)':>17s}\n"
        values = self.calculate(params)
        for point, mask, value in zip(self.data.points, self.data.mask, values):
            ncyc, intensity, error = point
            output += "#" if not mask else " "
            output += f" {ncyc:12d}  {intensity:17.8e} {error:17.8e} {value:17.8e}"
            output += " # NOT USED IN THE FIT\n" if not mask else "\n"
        return output + "\n\n"

    def filter(self, params):
        pass

    def plot(self, params, file_pdf, file_exp, file_fit, simulation=False):
        data_exp = self.data.get_r2_exp(simulation)
        intst_fit = self.calculate(params)
        data_fit = self.data.get_r2_fit(intst_fit, simulation)
        self._plot(file_pdf, self.name, data_exp, data_fit)
        output_fit = self._format_data_fit(data_fit)
        file_fit.write(output_fit + "\n\n")
        if not simulation:
            output_exp = self._format_data_exp(data_exp)
            file_exp.write(output_exp + "\n\n")

    def monte_carlo(self, params):
        intensities_ref = self.calculate(params)
        profile = copy.copy(self)
        profile.data = profile.data.monte_carlo(intensities_ref)
        return profile

    def bootstrap(self):
        """Make a profile for boostrap analysis."""
        profile = copy.copy(self)
        profile.data = profile.data.bootstrap()
        return profile

    def set_params(self, params, rates):
        for name1, name2 in self._par_names.items():
            name = cpn.remove_state(name1)
            if name in rates:
                params[name2].value = rates[name]

    def _get_parvals(self, params):
        parvals = tuple(
            (name1, params[name2].value) for name1, name2 in self._par_names.items()
        )
        return parvals

    def _format_data_exp(self, data_exp):
        result = f"[{self.name}]\n"
        result += f"# {'NU_CPMG':>12s}  {'R2 (EXP)':>17s} {'ERROR DOWN (EXP)':>17s} {'ERROR UP (EXP)':>17s}\n"
        for point in data_exp:
            nu_cpmgs = point["nu_cpmgs"]
            r2 = point["r2"]
            errors = point["errors"]
            result += (
                f"  {nu_cpmgs:12.3f}  {r2:17.8e} {errors[0]:17.8e} {errors[1]:17.8e}"
            )
            result += " # NOT USED IN THE FIT\n" if not point["mask"] else "\n"
        return result

    def _format_data_fit(self, data_fit):
        result = f"[{self.name}]\n"
        result += f"# {'NU_CPMG':>12s}  {'R2 (CALC)':>17s}\n"
        for point in data_fit:
            result += "  {nu_cpmgs:12.3f}  {r2:17.8e}\n".format_map(point)
        return result

    def any_duplicate(self):
        return self.data.any_duplicate()

    def __add__(self, other: "CpmgProfile"):
        data = self.data + other.data
        return CpmgProfile(
            self.name, data, self._pulse_seq, self._par_names, self.params_default
        )

    def __eq__(self, other: "CpmgProfile"):
        return self.name == other.name

    def __lt__(self, other: "CpmgProfile"):
        return self.name < other.name


class CpmgData:
    randm1 = np.random.randn(10000, 1)
    randm2 = np.random.randn(10000, 1)
    dtype = np.dtype([("ncycs", "i4"), ("intensities", "f8"), ("errors", "f8")])

    def __init__(self, points, refs, mask, time_t2):
        self.points = points
        self.refs = refs
        self.mask = mask
        self._time_t2 = time_t2

    @classmethod
    def from_file(cls, path, filter_planes, time_t2):
        try:
            points = np.loadtxt(path, dtype=cls.dtype)
        except OSError as err:
            sys.exit(f"\nerror: {err}")
        else:
            refs = points["ncycs"] == 0
            mask = np.array([True] * len(points))
            planes_to_filter = [
                index for index in filter_planes if 0 <= index < len(points)
            ]
            mask[planes_to_filter] = False
            return cls(points, refs, mask, time_t2)

    def estimate_noise_variance(self, kind):
        return ccn.estimate_noise_variance[kind](self.points)

    def monte_carlo(self, intensities_ref):
        noise = np.random.randn(len(self.points["intensities"])) * self.points["errors"]
        data = copy.deepcopy(self)
        data.points["intensities"] = intensities_ref + noise
        return data

    def bootstrap(self):
        indexes = np.arange(self.points["intensities"].size)
        pool1 = indexes[self.refs & self.mask]
        pool2 = indexes[~self.refs & self.mask]
        bs_indexes = []
        if pool1.size:
            bs_indexes.extend(np.random.choice(pool1, pool1.size))
        bs_indexes.extend(np.random.choice(pool2, pool2.size))
        bs_indexes = sorted(bs_indexes)
        data = copy.deepcopy(self)
        data.points = self.points[bs_indexes]
        return data

    def get_r2_exp(self, simulation=False):
        dtype = [
            ("nu_cpmgs", "f8"),
            ("r2", "f8"),
            ("errors", "f8", (2,)),
            ("mask", "?"),
        ]
        if simulation:
            return np.rec.array([[], [], [], []], dtype=dtype)
        points = self.points[~self.refs]
        points_ref = self.points[self.refs]
        nu_cpmgs = self._ncycs_to_nu_cpmg(points["ncycs"])
        r2 = self._intst_to_r2(points["intensities"], points_ref["intensities"])
        intst_ens = points["intensities"] + points["errors"] * self.randm1
        intst_ref_ens = points_ref["intensities"] + points_ref["errors"] * self.randm2
        r2_ens = self._intst_to_r2(intst_ens, intst_ref_ens)
        errors = np.percentile(r2_ens - r2, [15.9, 84.1], axis=0).transpose()
        mask = self.mask[~self.refs]
        r2_exp = np.rec.array([nu_cpmgs, r2, errors, mask], dtype=dtype)
        return np.sort(r2_exp, order="nu_cpmgs")

    def get_r2_fit(self, intst_fit, simulation=False):
        refs = self.refs
        intst = intst_fit[~refs]
        if simulation:
            intst_ref = intst_fit[refs]
        else:
            intst_ref = self.points[refs]["intensities"]
        nu_cpmgs = self._ncycs_to_nu_cpmg(self.points[~refs]["ncycs"])
        r2 = self._intst_to_r2(intst, intst_ref)
        r2_fit = np.rec.array([nu_cpmgs, r2], names=["nu_cpmgs", "r2"])
        return np.sort(r2_fit, order="nu_cpmgs")

    def _ncycs_to_nu_cpmg(self, ncycs=None):
        ncycs_ = np.array(ncycs, dtype=np.float)
        ncycs_[ncycs_ == -1.0] = 0.5
        return ncycs_[ncycs_ != 0.0] / self._time_t2

    def _intst_to_r2(self, intst, intst_ref):
        intst_norm = intst / np.mean(intst_ref, axis=-1, keepdims=True)
        r2 = np.full_like(intst, np.inf)
        neg = intst_norm <= 0.0
        r2[~neg] = -np.log(intst_norm[~neg]) / self._time_t2
        return r2

    def any_duplicate(self):
        return np.unique(self.points["ncycs"]).size != self.points.size

    def __add__(self, other: "CpmgData"):
        points = self.points.copy()
        points["intensities"] = self.points["intensities"] + other.points["intensities"]
        points["errors"] = np.sqrt(
            self.points["errors"] ** 2 + other.points["errors"] ** 2
        )
        refs = self.refs.copy()
        mask = self.mask.copy()
        return CpmgData(points, refs, mask)
