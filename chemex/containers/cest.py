import copy
import functools as ft
import sys

import numpy as np

import chemex.containers.helper as cch
import chemex.containers.noise as ccn
import chemex.parameters.name as cpn
import chemex.plot as cp


CEST_SCHEMA = {
    "type": "object",
    "properties": {
        "data": {
            "type": "object",
            "properties": {
                "error": {
                    "type": "string",
                    "enum": ["file", "scatter"],
                    "default": "file",
                },
                "filter_offsets": {
                    "type": "array",
                    "items": {
                        "type": "array",
                        "minItems": 2,
                        "maxItems": 2,
                        "items": {"type": "number"},
                    },
                    "default": [[0.0, 0.0]],
                },
                "filter_planes": {
                    "type": "array",
                    "items": {"type": "integer"},
                    "default": [],
                },
                "filter_ref_planes": {"type": "boolean", "default": False},
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
class CestProfile:
    def __init__(self, name, data, pulse_seq, par_names, params_default):
        self.name = name
        self.data = data
        self._pulse_seq = pulse_seq
        self._par_names = par_names
        self.params_default = params_default
        self._plot = cp.cest

    @classmethod
    def from_file(cls, path, config, pulse_seq, par_names, params_default):
        name = config["spin_system"]["spin_system"]
        data = CestData.from_file(
            path,
            filter_offsets=config["data"]["filter_offsets"],
            filter_planes=config["data"]["filter_planes"],
            filter_ref_planes=config["data"]["filter_ref_planes"],
        )
        return cls(name, data, pulse_seq, par_names, params_default)

    def residuals(self, params):
        data = self.data.points[self.data.mask]
        residuals = (self.calculate(params) - data["intensities"]) / data["errors"]
        return residuals

    def calculate(self, params, offsets=None):
        data = self.data.points[self.data.mask]
        par_values = self._get_parvals(params)
        calculated = self._pulse_seq.calculate(tuple(data["offsets"]), par_values)
        scale = cch.get_scale(data["intensities"], data["errors"], calculated)
        if offsets is not None:
            calculated = self._pulse_seq.calculate(tuple(offsets), par_values)
        return scale * calculated

    def estimate_noise_variance(self, kind):
        return self.data.estimate_noise_variance(kind)

    def set_noise(self, value):
        self.data.points["errors"] = value

    def print(self, params):
        output = f"[{self.name}]\n"
        output += f"# {'OFFSET (HZ)':>12s}  {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s} {'INTENSITY (CALC)':>17s}\n"
        values = self.calculate(params, self.data.points["offsets"])
        for point, mask, value in zip(self.data.points, self.data.mask, values):
            offset, intensity, error = point
            output += "#" if not mask else " "
            output += (
                f" {offset: 12.2f}  {intensity: 17.8e} {error: 17.8e} {value: 17.8e}"
            )
            output += " # NOT USED IN THE FIT\n" if not mask else "\n"
        return output + "\n\n"

    def filter(self, params):
        cs_values = self._get_cs_values(params)
        cs_offset = self._pulse_seq.ppms_to_offsets(cs_values)[0]
        sw_dante = getattr(self._pulse_seq, "sw_dante", None)
        self.data.filter(cs_offset, sw_dante)

    def plot(self, params, file_pdf):
        data_exp, data_fit = self._get_plot_data(params)
        cs_values = self._get_cs_values(params)
        self._plot(file_pdf, self.name, data_exp, data_fit, cs_values)

    def write_plot(self, params, file_exp, file_fit):
        data_exp, data_fit = self._get_plot_data(params)
        output_exp = f"[{self.name}]\n"
        output_exp += (
            f"# {'CS (PPM)':>12s}  {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s}\n"
        )
        for point in data_exp:
            nu_cpmgs = point["ppms"]
            intensities = point["intensities"]
            errors = point["errors"]
            output_exp += f"  {nu_cpmgs:12.2f}  {intensities:17.8e} {errors[1]:17.8e}"
            output_exp += " # NOT USED IN THE FIT" if not point["mask"] else "\n"
        file_exp.write(output_exp + "\n\n")
        output_fit = f"[{self.name}]\n"
        output_fit += f"# {'CS (PPM)':>12s}  {'INTENSITY (CALC)':>17s}\n"
        for point in data_fit:
            output_fit += "  {ppms: 12.2f}  {intensities: 17.8e}\n".format_map(point)
        file_fit.write(output_fit + "\n\n")

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

    def _get_cs_values(self, params):
        names = (f"cs_i_{state}" for state in "abcd")
        fnames = (self._par_names[name] for name in names if name in self._par_names)
        return [params[fname] for fname in fnames]

    def _get_plot_data(self, params):
        data = self.data.points[~self.data.refs]
        data_refs = self.data.points[self.data.refs]
        mask = self.data.mask[~self.data.refs]
        intstref = np.mean(data_refs["intensities"])
        ppms = self._pulse_seq.offsets_to_ppms(data["offsets"])
        intensities = data["intensities"] / intstref
        errors = data["errors"] / abs(intstref)
        errors = np.array([-errors, errors]).transpose()
        offsets_fit = cp.get_grid(data["offsets"], 500, 0.02)
        ppms_fit = self._pulse_seq.offsets_to_ppms(offsets_fit)
        intst_fit = self.calculate(params, offsets_fit) / intstref
        data_exp = np.rec.array(
            [ppms, intensities, errors, mask],
            dtype=[
                ("ppms", "f8"),
                ("intensities", "f8"),
                ("errors", "f8", (2,)),
                ("mask", "?"),
            ],
        )
        data_fit = np.rec.array([ppms_fit, intst_fit], names=["ppms", "intensities"])
        data_exp = np.sort(data_exp, order="ppms")
        data_fit = np.unique(np.sort(data_fit, order="ppms"))
        return data_exp, data_fit

    def any_duplicate(self):
        return self.data.any_duplicate()

    def __add__(self, other: "CestProfile"):
        data = self.data + other.data
        return CestProfile(
            self.name, data, self._pulse_seq, self._par_names, self.params_default
        )

    def __eq__(self, other: "CestProfile"):
        return self.name == other.name

    def __lt__(self, other: "CestProfile"):
        return self.name < other.name


class CestData:
    dtype = np.dtype([("offsets", "f8"), ("intensities", "f8"), ("errors", "f8")])

    def __init__(self, points, refs, mask, filter_offsets):
        self.points = points
        self.refs = refs
        self.mask = mask
        self._filter_offsets = filter_offsets

    @classmethod
    def from_file(cls, path, filter_offsets, filter_planes, filter_ref_planes):
        try:
            points = np.loadtxt(path, dtype=cls.dtype)
        except OSError as err:
            sys.exit(f"\nerror: {err}")
        else:
            refs = abs(points["offsets"]) >= 1.0e4
            mask = np.array([True] * len(points))
            planes_to_filter = [
                index for index in filter_planes if 0 <= index < len(points)
            ]
            mask[planes_to_filter] = False
            if filter_ref_planes:
                mask[refs] = False
            return cls(points, refs, mask, filter_offsets)

    def estimate_noise_variance(self, kind):
        return ccn.estimate_noise_variance[kind](self.points)

    def filter(self, cs_offset, sw_dante=None):
        offsets = self.points["offsets"] - cs_offset
        for filter_offset, filter_bandwidth in self._filter_offsets:
            offsets_ = offsets - filter_offset
            if sw_dante is not None:
                offsets_ = (offsets_ + 0.5 * sw_dante) % sw_dante - 0.5 * sw_dante
            mask_filter = abs(offsets_) < filter_bandwidth * 0.5
            self.mask[mask_filter] = False

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

    def any_duplicate(self):
        return np.unique(self.points["offsets"]).size != self.points.size

    def __add__(self, other: "CestData"):
        points = self.points.copy()
        points["intensities"] = self.points["intensities"] + other.points["intensities"]
        points["errors"] = np.sqrt(
            (self.points["errors"] ** 2 + other.points["errors"] ** 2)
        )
        refs = self.refs.copy()
        mask = self.mask.copy()
        filter_offsets = self._filter_offsets
        return CestData(points, refs, mask, filter_offsets)
