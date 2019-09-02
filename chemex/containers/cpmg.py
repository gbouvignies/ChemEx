import copy

import numpy as np

import chemex.containers.helper as cch
import chemex.containers.noise as ccn
import chemex.nmr.helper as cnsn
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
        data = CpmgData.from_file(path)
        return cls(name, data, pulse_seq, par_names, params_default)

    def residuals(self, params):
        data = self.data.points
        mask = self.data.mask
        residuals = (self.calculate(params) - data["intensities"]) / data["errors"]
        return residuals[mask]

    def calculate(self, params, ncycs=None):
        data = self.data.points
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
        output = f"[{str(self.name).upper()}]\n"
        output += "# {:>12s}  {:>17s} {:>17s} {:>17s}\n".format(
            "NCYC", "INTENSITY (EXP)", "ERROR (EXP)", "Intensity (calc)"
        )
        values = self.calculate(params)
        for point, value in zip(self.data.points, values):
            ncyc, intensity, error = point
            output += f"  {ncyc:12d}  {intensity:17.8e} {error:17.8e} {value:17.8e}\n"

        return output + "\n\n"

    def filter(self, params):
        pass

    def plot(self, file_pdf, params):
        data_exp, data_fit = self._get_plot_data(params)
        self._plot(file_pdf, self.name, data_exp, data_fit)

    def write_plot(self, file_exp, file_fit, params):
        data_exp, data_fit = self._get_plot_data(params)
        output_exp = f"[{self.name}]\n"
        output_exp += (
            "# {'NU_CPMG':>12s}  {'R2 (EXP)':>17s} {'ERROR DOWN (EXP)':>17s} "
            "{'ERROR UP (EXP)':>17s}\n"
        )
        for point in data_exp:
            nu_cpmgs = point["nu_cpmgs"]
            r2 = point["r2"]
            errors = point["errors"]
            output_exp += (
                f"  {nu_cpmgs:12.3f}  {r2:17.8e} {errors[0]:17.8e} {errors[1]:17.8e}\n"
            )
        file_exp.write(output_exp + "\n\n")
        output_fit = f"[{self.name}]\n"
        output_fit += "# {'NU_CPMG':>12s}  {'R2 (CALC)':>17s}\n"
        for point in data_fit:
            output_fit += "  {nu_cpmgs:12.3f}  {r2:17.8e}\n".format_map(point)
        file_fit.write(output_fit + "\n\n")

    def monte_carlo(self, params):
        intensities_ref = self.calculate(params)
        profile = copy.deepcopy(self)
        profile.data = profile.data.monte_carlo(intensities_ref)
        return profile

    def bootstrap(self):
        """Make a profile for boostrap analysis."""
        profile = copy.deepcopy(self)
        profile.data = profile.data.bootstrap()
        return profile

    def _get_parvals(self, params):
        parvals = tuple(
            (name1, params[name2].value) for name1, name2 in self._par_names.items()
        )
        return parvals

    def _get_plot_data(self, params):
        time_t2 = self._pulse_seq.time_t2
        nu_cpmgs, r2_exp, r2_err, mask = self.data.to_r2(time_t2)
        intst_fit = self.calculate(params)
        refs = self.data.refs
        intst_ref = self.data.points[refs]["intensities"]
        r2_fit = intensities_to_r2(intst_fit[~refs], intst_ref, time_t2)
        data_exp = np.rec.array(
            [nu_cpmgs, r2_exp, r2_err, mask],
            dtype=[
                ("nu_cpmgs", "f8"),
                ("r2", "f8"),
                ("errors", "f8", (2,)),
                ("mask", "?"),
            ],
        )
        data_fit = np.rec.array([nu_cpmgs, r2_fit], names=["nu_cpmgs", "r2"])
        data_exp = np.sort(data_exp, order="nu_cpmgs")
        data_fit = np.unique(np.sort(data_fit, order="nu_cpmgs"))
        return data_exp, data_fit


class CpmgData:
    randm1 = np.random.randn(10000, 1)
    randm2 = np.random.randn(10000, 1)
    dtype = np.dtype([("ncycs", "i4"), ("intensities", "f8"), ("errors", "f8")])

    def __init__(self, points, refs, mask):
        self.points = points
        self.refs = refs
        self.mask = mask

    @classmethod
    def from_file(cls, path):
        try:
            points = np.loadtxt(path, dtype=cls.dtype)
        except OSError as err:
            exit(f"\nerror: {err}")
        else:
            refs = points["ncycs"] == 0
            mask = np.array([True] * len(points))
            return cls(points, refs, mask)

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

    def to_r2(self, time_t2):
        points = self.points[~self.refs]
        points_ref = self.points[self.refs]
        intst = points["intensities"]
        intst_ref = points_ref["intensities"]
        nu_cpmgs = ncycs_to_nu_cpmgs(points["ncycs"], time_t2)
        r2 = intensities_to_r2(intst, intst_ref, time_t2)
        intst_ens = intst + points["errors"] * self.randm1
        intst_ref_ens = intst_ref + points_ref["errors"] * self.randm2
        r2_ens = intensities_to_r2(intst_ens, intst_ref_ens, time_t2) - r2
        r2_err = np.percentile(r2_ens, [15.9, 84.1], axis=0).transpose()
        return nu_cpmgs, r2, r2_err, self.mask[~self.refs]


def intensities_to_r2(intst, intst_ref, time_t2):
    r2 = np.zeros_like(intst)
    intst_norm = intst / np.mean(intst_ref, axis=-1, keepdims=True)
    neg = intst_norm <= 0.0
    r2[neg] = np.inf
    r2[~neg] = -np.log(intst_norm[~neg]) / time_t2
    return r2


def ncycs_to_nu_cpmgs(ncycs, time_t2):
    return np.asarray(ncycs) / time_t2
