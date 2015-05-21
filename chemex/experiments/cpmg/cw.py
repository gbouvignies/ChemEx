"""Pure In-phase Nitrogen CPMG (2-state)

Analyzes 15N chemical exchange in the presence of high power 1H CW decoupling
during the CPMG block. This keeps the spin system purely in-phase throughout,
and is calculated using the 6x6, single spin matrix:

[ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b) ]

Notes
-----

Off resonance effects are taken into account.

The calculation is designed specifically to analyze the experiment found in
the reference:

Journal of Physical Chemistry B (2008), 112, 5898-5904
"""

import functools

import lmfit
import numpy as np
from scipy import linalg

from chemex import constants
from chemex.bases.two_states import single_spin
from chemex.experiments import base_profile, util, sputil
from chemex.experiments.cpmg import plotting

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

reduce = functools.reduce
dot = np.dot
matrix_power = np.linalg.matrix_power
calculate_shift_ex_2st = util.calculate_shift_ex_2st
compute_liouvillian = single_spin.compute_liouvillian
check_par = base_profile.check_par
compute_propagator_from_time_series = single_spin.compute_propagators_from_time_series

two_pi = 2.0 * np.pi

attributes_exp = {
    'h_larmor_frq': float,
    'temperature': float,
    'carrier': float,
    'pw': float,
    'time_t2': float,
    'time_equil': float,
}

params_exp = {
    'pb': {
        'value': 0.05,
        'vary': True,
        'min': 0.0,
        'max': 1.0,
        'expr': None,
    },
    'kex': {
        'value': 100.0,
        'vary': True,
        'min': 0.0,
        'max': None,
        'expr': None,
    },
    'cs': {
        'value': 0.0,
        'vary': False,
        'min': None,
        'max': None,
        'expr': None,
    },
    'dw': {
        'value': 0.0,
        'vary': True,
        'min': None,
        'max': None,
        'expr': None,
    },
    'r_ixy': {
        'value': 10.0,
        'vary': True,
        'min': 0.0,
        'max': None,
        'expr': None,
    },
    'dr_ixy': {
        'value': 0.0,
        'vary': False,
        'min': None,
        'max': None,
        'expr': None,
    },
    'r_iz': {
        'value': 1.0,
        'vary': False,
        'min': 0.0,
        'max': None,
        'expr': None,
    },
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.ncycs = measurements['ncycs']
        self.val = measurements['intensities']
        self.err = measurements['intensities_err']

        for name, convert in attributes_exp.items():
            setattr(self, name, check_par(exp_details, name, convert))

        self.experiment_name = check_par(exp_details, 'experiment_name')

        self._calculate_profile = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = sputil.Peak(self.profile_name)
        self.resonance = self.peak.resonances[0]

        self.ppm_to_rads = (
            self.h_larmor_frq * two_pi *
            constants.xi_ratio[self.resonance.atom.nucleus]
        )

        temperature_str = str(self.temperature).replace('.', '_p_')
        h_larmor_frq_str = str(self.h_larmor_frq).replace('.', '_p_')
        resonance_str = self.resonance.name.replace('-', '_m_')

        param_names = (
            ('pb', temperature_str),
            ('kex', temperature_str),
            ('cs', resonance_str, temperature_str),
            ('dw', resonance_str, temperature_str),
            ('r_ixy', resonance_str, temperature_str, h_larmor_frq_str),
            ('dr_ixy', resonance_str, temperature_str, h_larmor_frq_str),
            ('r_iz', resonance_str, temperature_str, h_larmor_frq_str),
        )

        self.map_names = {name[0]: '__'.join(name) for name in param_names}

    def make_default_parameters(self):

        parameters = lmfit.Parameters()

        for name, settings in params_exp.items():
            parameters.add(self.map_names[name], **settings)

        return parameters

    def _calculate_profile(self, pb, kex, dw, r_iz, r_ixy, dr_ixy, cs):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb : float
            Fractional population of state B.
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_nz : float
            Longitudinal relaxation rate of states A and B in /s.
        r_nxy : float
            Transverse relaxation rate of state A in /s.
        dr_nxy : float
            Transverse relaxation rate difference between states A and B in /s.
        cs : float
            Resonance position in ppm.

        Returns
        -------
        out : float
            Intensity after the CEST block
        """

        dw_ = dw * self.ppm_to_rads
        w = (cs - self.carrier) * self.ppm_to_rads

        pa = 1.0 - pb
        mag_eq = np.array([[0.0, 0.0, pa,
                            0.0, 0.0, pb]]).T

        profile = []

        w1 = 2.0 * np.pi / (4.0 * self.pw)

        # Calculation of the different liouvillians

        kwargs = {'pb': pb, 'kex': kex, 'dw': dw_, 'r_ixy': r_ixy,
                  'dr_ixy': dr_ixy, 'r_iz': r_iz, 'w': w}

        l_free = single_spin.compute_liouvillian(**kwargs)
        l_pw1x = single_spin.compute_liouvillian(w1x=+w1, **kwargs)
        l_pw1y = single_spin.compute_liouvillian(w1y=+w1, **kwargs)
        l_mw1x = single_spin.compute_liouvillian(w1x=-w1, **kwargs)

        # Calculation of all the needed propagators

        p_90px = linalg.expm(l_pw1x * self.pw)
        p_90py = linalg.expm(l_pw1y * self.pw)
        p_90mx = linalg.expm(l_mw1x * self.pw)
        p_180pmx = 0.5 * (matrix_power(p_90px, 2) +
                          matrix_power(p_90mx, 2))
        p_180py = matrix_power(p_90py, 2)

        t_neg = -2.0 * self.pw / np.pi
        time_series = [t_neg, self.time_equil]
        tau_cps = self.time_t2 / (4.0 * self.ncycs) - self.pw
        time_series.extend(tau_cps)

        propagators = compute_propagator_from_time_series(l_free, time_series)

        p_equil = propagators[self.time_equil]
        p_neg = propagators[t_neg]

        # Simulate the CPMG block as function of ncyc

        for ncyc in self.ncycs:

            if ncyc == 0:

                mag = reduce(dot, [p_equil, p_90px, p_180pmx, p_90px, mag_eq])

            else:

                tau_cp = self.time_t2 / (4.0 * ncyc) - self.pw
                p_free = propagators[tau_cp]
                p_cp = matrix_power(dot(dot(p_free, p_180py), p_free), int(ncyc))
                mag = reduce(dot, [p_equil, p_90px, p_neg, p_cp, p_180pmx, p_cp, p_neg, p_90px, mag_eq])

            profile.append(mag[2, 0])

        return np.asarray(profile)

    def calculate_profile(self, params, ncycs=None):

        kwargs = {
            short_name: params[long_name].value
            for short_name, long_name in self.map_names.items()
            }

        values = self._calculate_profile(**kwargs)
        scale = self._calculate_scale(values)

        if ncycs is not None:
            self.ncycs, ncycs_orig = ncycs, self.ncycs
            values = self._calculate_profile.__wrapped__(**kwargs)
            self.ncycs = ncycs_orig

        return values * scale

    def _calculate_scale(self, cal):

        scale = (
            sum(cal * self.val / self.err ** 2) /
            sum((cal / self.err) ** 2)
        )

        return scale

    def calculate_residuals(self, params):
        """Calculates the residual between the experimental and
        back-calculated values.
        """

        values = self.calculate_profile(params)

        return (self.val - values) / self.err

    def ncycs_to_nu_cpmgs(self, ncycs=None):

        if ncycs is None:
            ncycs = self.ncycs

        return ncycs / self.time_t2

    def filter_points(self, params):
        """Evaluate some criteria to know whether the point should be considered
        in the calculation or not.

        Returns 'True' if the point should NOT be considered.
        """

        return False

    def print_profile(self, params=None):
        """Print the data point"""

        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.val

        iter_vals = zip(self.ncycs, self.val, self.err, values)

        for ncyc, val, err, cal in iter_vals:

            line = (
                "{0.profile_name:10s} "
                "{0.h_larmor_frq:8.1f} "
                "{0.time_t2:8.1e} "
                "{0.temperature:5.1f} "
                "{1:15.8e} "
                "{2:15.8e} "
                "{3:15.8e} "
                    .format(self, ncyc, val, err)
            )

            if params is not None:
                line += "{:15.8e}".format(cal)

            output.append(line)

        output.append("")

        return "\n".join(output).upper()
