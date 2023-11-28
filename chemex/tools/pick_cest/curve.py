import numpy as np
from scipy.interpolate import CubicSpline

from chemex.containers.profile import Profile
from chemex.plotters.cest import create_plot_data_exp


class Curve:
    def __init__(self, profile: Profile, sw: float | None = None) -> None:
        data = create_plot_data_exp(profile)
        settings = getattr(profile.pulse_sequence, "settings", None)

        delta_ppm = 0.0
        if settings is not None:
            spectrometer = profile.spectrometer
            cos_n = getattr(settings, "cos_n", None)
            sw = getattr(settings, "sw", None)
            if cos_n is not None and cos_n % 2 == 0 and sw is not None:
                shifts_ppm = spectrometer.offsets_to_ppms(np.array([sw / 2.0, 0.0]))
                delta_ppm = shifts_ppm[0] - shifts_ppm[1]

        ppms = data.metadata - delta_ppm
        xmin, xmax = min(ppms), max(ppms)
        self.xrange = xmax - xmin
        self.xcentre = 0.5 * (xmin + xmax)
        self.x = ppms
        self.y = data.exp
        _, indices = np.unique(ppms, return_index=True)
        x = ppms[indices]
        y = data.exp[indices]
        if sw is not None:
            y[0] = y[-1] = 0.5 * (y[0] + y[-1])
            self.spline = CubicSpline(x, y, bc_type="periodic", extrapolate="periodic")
        else:
            self.spline = CubicSpline(x, y)

    def get_xrange(self, sw: float | None = None):
        if sw is None:
            sw = 1.0
        return self.xcentre + sw * self.xrange * np.array([-0.5, 0.5])
