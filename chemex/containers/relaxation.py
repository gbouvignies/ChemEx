import copy
import functools as ft
import sys

import numpy as np

import chemex.containers.helper as cch
import chemex.containers.noise as ccn
import chemex.containers.plot as ccp


@ft.total_ordering
class RelaxationProfile:
    def __init__(self, name, data, pulse_seq, pnames, params, params_mf):
        self.name = name
        self.data = data
        self._pulse_seq = pulse_seq
        self._pnames = pnames
        self.params = params
        self.params_mf = params_mf
        self._plot = ccp.relaxation

    @classmethod
    def from_file(cls, path, config, pulse_seq, pnames, params, params_mf):
        name = config["spin_system"]
        data = RelaxationData.from_file(
            path, filter_planes=config["data"]["filter_planes"]
        )
        return cls(name, data, pulse_seq, pnames, params, params_mf)

    def residuals(self, params):
        data = self.data.points[self.data.mask]
        return (self.calculate(params) - data["intensities"]) / data["errors"]

    def calculate(self, params, times=None):
        data = self.data.points[self.data.mask]
        par_values = self._get_parvals(params)
        calculated = self._pulse_seq.calculate(tuple(data["times"]), par_values)
        scale = cch.get_scale(data["intensities"], data["errors"], calculated)
        if times is not None:
            calculated = self._pulse_seq.calculate(tuple(times), par_values)
        return scale * calculated

    def estimate_noise_variance(self, kind):
        return self.data.estimate_noise_variance(kind)

    def set_noise(self, value):
        self.data.points["errors"] = value

    def print(self, params):
        output = f"[{self.name}]\n"
        output += (
            f"# {'TIMES (S)':>12s}  "
            f"{'INTENSITY (EXP)':>17s} "
            f"{'ERROR (EXP)':>17s} "
            f"{'INTENSITY (CALC)':>17s}\n"
        )
        values = self.calculate(params, self.data.points["times"])
        for point, mask, value in zip(self.data.points, self.data.mask, values):
            offset, intensity, error = point
            output += "#" if not mask else " "
            output += (
                f" {offset: 12.2f}  {intensity: 17.8e} {error: 17.8e} {value: 17.8e}"
            )
            output += " # NOT USED IN THE FIT\n" if not mask else "\n"
        return output + "\n\n"

    def filter(self, params):
        pass

    def plot(self, params, file_pdf, file_exp, file_fit, simulation=False):
        data_exp = self._get_plot_data_exp(simulation)
        data_fit = self._get_plot_data_fit(params, simulation)
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
        """Make a profile for bootstrap analysis."""
        profile = copy.copy(self)
        profile.data = profile.data.bootstrap()
        return profile

    def _get_parvals(self, params):
        return tuple(
            (name1, params[name2].value) for name1, name2 in self._pnames.items()
        )

    def _get_plot_data_exp(self, simulation=False):
        dtype = [
            ("times", "f8"),
            ("intensities", "f8"),
            ("errors", "f8", (2,)),
            ("mask", "?"),
        ]
        if simulation:
            return np.rec.array([[], [], [], []], dtype=dtype)
        points = self.data.points
        times = points["times"]
        intst_ref = points["intensities"][np.argmax(np.abs(points["intensities"]))]
        intensities = points["intensities"] / intst_ref
        errors = points["errors"] / abs(intst_ref)
        errors = np.array([-errors, errors]).transpose()
        mask = self.data.mask
        data_exp = np.rec.array([times, intensities, errors, mask], dtype=dtype)
        return np.sort(data_exp, order="times")

    def _get_plot_data_fit(self, params, simulation=False):
        data = self.data.points
        if simulation:
            intst_calc = self.calculate(params)
            intst_ref = intst_calc[np.argmax(np.abs(intst_calc))]

        else:
            points = self.data.points
            intst_ref = points["intensities"][np.argmax(np.abs(points["intensities"]))]
        times_fit = cch.get_grid(data["times"], 100, 0.02)
        intensities = self.calculate(params, times_fit) / intst_ref
        data_fit = np.rec.array(
            [times_fit, intensities], names=["times", "intensities"]
        )
        return np.sort(data_fit, order="times")

    def _format_data_exp(self, data_exp):
        result = f"[{self.name}]\n"
        result += (
            f"# {'TIMES (S)':>12s}  {'INTENSITY (EXP)':>17s} {'ERROR (EXP)':>17s}\n"
        )
        for point in data_exp:
            times = point["times"]
            intensities = point["intensities"]
            errors = point["errors"]
            result += f"  {times:12.2f}  {intensities:17.8e} {errors[1]:17.8e}"
            result += " # NOT USED IN THE FIT" if not point["mask"] else "\n"
        return result

    def _format_data_fit(self, data_fit):
        result = f"[{self.name}]\n"
        result += f"# {'TIMES (S)':>12s}  {'INTENSITY (CALC)':>17s}\n"
        for point in data_fit:
            result += "  {times: 12.2f}  {intensities: 17.8e}\n".format_map(point)
        return result

    def any_duplicate(self):
        return self.data.any_duplicate()

    def __add__(self, other: object):
        if not isinstance(other, type(self)):
            return NotImplemented
        data = self.data + other.data
        return RelaxationProfile(
            self.name, data, self._pulse_seq, self._pnames, self.params, self.params_mf
        )

    def __eq__(self, other: object):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name == other.name

    def __lt__(self, other: object):
        if not isinstance(other, type(self)):
            return NotImplemented
        return self.name < other.name


class RelaxationData:
    dtype = np.dtype([("times", "f8"), ("intensities", "f8"), ("errors", "f8")])

    def __init__(self, points, mask):
        self.points = points
        self.mask = mask

    @classmethod
    def from_file(cls, path, filter_planes):
        try:
            points = np.loadtxt(path, dtype=cls.dtype)
        except OSError as err:
            sys.exit(f"\nerror: {err}")
        else:
            mask = np.array([True] * len(points))
            planes_to_filter = [
                index for index in filter_planes if 0 <= index < len(points)
            ]
            mask[planes_to_filter] = False
            return cls(points, mask)

    def estimate_noise_variance(self, kind):
        return ccn.estimate_noise_variance[kind](self.points)

    def monte_carlo(self, intensities_ref):
        noise = np.random.randn(len(self.points["intensities"])) * self.points["errors"]
        data = copy.deepcopy(self)
        data.points["intensities"][self.mask] = intensities_ref + noise[self.mask]
        return data

    def bootstrap(self):
        indexes = np.arange(self.points["intensities"].size)
        pool = indexes[self.mask]
        bs_indexes = np.random.choice(pool, pool.size)
        bs_indexes = sorted(bs_indexes)
        data = copy.deepcopy(self)
        data.points[self.mask] = self.points[bs_indexes]
        return data

    def any_duplicate(self):
        return np.unique(self.points["times"]).size != self.points.size

    def __add__(self, other: "RelaxationData"):
        points = self.points.copy()
        points["intensities"] = self.points["intensities"] + other.points["intensities"]
        points["errors"] = np.sqrt(
            self.points["errors"] ** 2 + other.points["errors"] ** 2
        )
        mask = self.mask.copy()
        return RelaxationData(points, mask)
