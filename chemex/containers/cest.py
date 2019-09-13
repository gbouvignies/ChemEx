import copy
import sys

import numpy as np

import chemex.containers.helper as cch
import chemex.containers.noise as ccn
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
        data = CestData.from_file(path, filter_offsets=config["data"]["filter_offsets"])
        return cls(name, data, pulse_seq, par_names, params_default)

    def residuals(self, params):
        data = self.data.points
        mask = self.data.mask
        residuals = (self.calculate(params) - data["intensities"]) / data["errors"]
        return residuals[mask]

    def calculate(self, params, offsets=None):
        data = self.data.points
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
        output += (
            "# {'OFFSET (HZ)':>12s}  {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s} "
            "{'INTENSITY (CALC)':>17s}\n"
        )
        values = self.calculate(params)
        for point, value in zip(self.data.points, values):
            offset, intensity, error = point
            output += (
                f"  {offset: 12.2f}  {intensity: 17.8e} {error: 17.8e} {value: 17.8e}\n"
            )

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
            output_exp += f"  {nu_cpmgs:12.2f}  {intensities:17.8e} {errors[1]:17.8e}\n"
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
        intensity_ref = np.mean(data_refs["intensities"])
        ppms = self._pulse_seq.offsets_to_ppms(data["offsets"])
        intensities = data["intensities"] / intensity_ref
        errors = data["errors"] / abs(intensity_ref)
        errors = np.array([-errors, errors]).transpose()
        offsets_fit = cp.get_grid(data["offsets"], 500, 0.02)
        ppms_fit = self._pulse_seq.offsets_to_ppms(offsets_fit)
        intensities_fit = self.calculate(params, offsets_fit) / intensity_ref
        data_exp = np.rec.array(
            [ppms, intensities, errors, mask],
            dtype=[
                ("ppms", "f8"),
                ("intensities", "f8"),
                ("errors", "f8", (2,)),
                ("mask", "?"),
            ],
        )
        data_fit = np.rec.array(
            [ppms_fit, intensities_fit], names=["ppms", "intensities"]
        )
        data_exp = np.sort(data_exp, order="ppms")
        data_fit = np.unique(np.sort(data_fit, order="ppms"))
        return data_exp, data_fit


class CestData:
    dtype = np.dtype([("offsets", "f8"), ("intensities", "f8"), ("errors", "f8")])

    def __init__(self, points, refs, mask, filter_offsets):
        self.points = points
        self.refs = refs
        self.mask = mask
        self._filter_offsets = filter_offsets

    @classmethod
    def from_file(cls, path, filter_offsets):
        try:
            points = np.loadtxt(path, dtype=cls.dtype)
        except OSError as err:
            sys.exit(f"\nerror: {err}")
        else:
            refs = abs(points["offsets"]) >= 1.0e4
            mask = np.array([True] * len(points))
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
