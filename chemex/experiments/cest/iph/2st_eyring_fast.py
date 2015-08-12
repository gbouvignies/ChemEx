"""Pure In-phase CEST (2-state, fast)

Analyzes chemical exchange in the presence of 1H composite decoupling
during the CEST block. This keeps the spin system purely in-phase throughout,
and is calculated using the 6x6, single spin matrix:

[ Ix(a), Iy(a), Iz(a), Ix(b), Iy(b), Iz(b) ]

Notes
-----
The calculation is designed specifically to analyze the experiment found in
the reference:

J Am Chem Soc (2012), 134, 8148-61

"""

import lmfit
import numpy as np
from scipy import linalg, constants as constants_sp

from chemex import constants, parameters, peaks
from chemex.bases import iph_2st, util
from chemex.experiments import base_profile
from chemex.experiments.cest import plotting

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

calculate_shift_2st = util.calculate_shift_ex_2st
compute_liouvillian = iph_2st.compute_liouvillian
check_par = base_profile.check_par
ParameterName = parameters.ParameterName

two_pi = 2.0 * np.pi
h_planck = constants_sp.h
kb = constants_sp.k
r_gas = constants_sp.R

attributes_exp = {
    'h_larmor_frq': float,
    'temperature': float,
    'carrier': float,
    'b1_frq': float,
    'time_t1': float,
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.b1_offsets = measurements['b1_offsets']
        self.val = measurements['intensities']
        self.err = measurements['intensities_err']

        self.h_larmor_frq = check_par(exp_details, 'h_larmor_frq', float)
        self.temperature = check_par(exp_details, 'temperature', float)
        self.carrier = check_par(exp_details, 'carrier', float)
        self.b1_frq = check_par(exp_details, 'b1_frq', float)
        self.time_t1 = check_par(exp_details, 'time_t1', float)

        self.experiment_name = check_par(exp_details, 'experiment_name')
        self.on_resonance_filter = check_par(exp_details, 'on_resonance_filter', convert=float, default=0.0)

        self._calculate_profile = lru_cache(5)(self._calculate_profile)
        self.plot_data = plotting.plot_data

        self.peak = peaks.Peak(self.profile_name)
        self.resonance = self.peak.resonances[0]

        self.ppm_to_rads = (
            self.h_larmor_frq * two_pi *
            constants.xi_ratio[self.resonance['atom']]
        )

        kwargs1 = {'temperature': self.temperature, 'nuclei': self.resonance['name']}
        kwargs2 = {'temperature': self.temperature, 'nuclei': self.resonance['name'], 'h_larmor_frq': self.h_larmor_frq}

        self.map_names = {
            'dh_b': ParameterName('dh_b').to_full_name(),
            'ds_b': ParameterName('ds_b').to_full_name(),
            'dh_ab': ParameterName('dh_ab').to_full_name(),
            'ds_ab': ParameterName('ds_ab').to_full_name(),
            'cs_i_a': ParameterName('cs_a', **kwargs1).to_full_name(),
            'dw_i_ab': ParameterName('dw_ab', **kwargs1).to_full_name(),
            'lambda_i_a': ParameterName('lambda_a', **kwargs2).to_full_name(),
            'dlambda_i_ab': ParameterName('dlambda_ab', **kwargs2).to_full_name(),
            'rho_i_a': ParameterName('rho_a', **kwargs2).to_full_name(),
        }

    def create_default_parameters(self):

        parameters = lmfit.Parameters()

        parameters.add_many(
            # Name, Value, Vary, Min, Max, Expr
            (self.map_names['dh_b'], 6000.0, True, None, None, None),
            (self.map_names['ds_b'], 0.0, False, None, None, None),
            (self.map_names['dh_ab'], 60000.0, True, None, None, None),
            (self.map_names['ds_ab'], 0.0, False, None, None, None),
            (self.map_names['cs_i_a'], 0.0, False, None, None, None),
            (self.map_names['dw_i_ab'], 0.0, True, None, None, None),
            (self.map_names['lambda_i_a'], 10.0, True, 0.0, None, None),
            (self.map_names['dlambda_i_ab'], 0.0, True, None, None, None),
            (self.map_names['rho_i_a'], 1.0, True, 0.0, None, None),
        )

        return parameters

    def _calculate_profile(self, dh_b, ds_b, dh_ab, ds_ab, dw_i_ab, rho_i_a, lambda_i_a, dlambda_i_ab, cs_i_a):
        """Calculate the intensity in presence of exchange after a CEST block.

        Parameters
        ----------
        pb : float
            Fractional population of state B.
        kex_ab : float
            Exchange rates between states A and B in /s.
        dw_i_ab : float
            Chemical shift difference between states A and B in rad/s.
        rho_i_a : float
            Longitudinal relaxation rate of states A in /s.
        lambda_i_a : float
            Transverse relaxation rate of state A in /s.
        dlambda_i_ab : float
            Transverse relaxation rate difference between states A and B in /s.
        cs_i_a : float
            Resonance position of state A in ppm.

        Returns
        -------
        out : float
            Intensity after the CEST block
        """

        t_kelvin = self.temperature + 273.15
        dg_b = dh_b - t_kelvin * ds_b
        dg_ab = dh_ab - t_kelvin * ds_ab
        dg = np.asarray([dg_ab, dg_ab - dg_b])
        kab, kba = kb * t_kelvin / h_planck * np.exp(-dg / (r_gas * t_kelvin))
        kex_ab = kab + kba
        pb = kab / kex_ab

        omega_i_a_array = (cs_i_a - self.carrier) * self.ppm_to_rads - two_pi * self.b1_offsets
        domega_i_ab = dw_i_ab * self.ppm_to_rads
        omega1x_i = two_pi * self.b1_frq
        lambda_i_b = lambda_i_a + dlambda_i_ab

        # Correct chemical shift against exchange induced shift

        shift_ex, _ = calculate_shift_2st(
            pb=pb,
            kex_ab=kex_ab,
            domega_i_ab=domega_i_ab,
            lambda_i_a=lambda_i_a,
            lambda_i_b=lambda_i_b
        )

        omega_i_a_array -= shift_ex

        magz_eq = np.array([[1 - pb], [pb]])

        profile = []

        for b1_offset, omega_i_a in zip(self.b1_offsets, omega_i_a_array):

            if b1_offset <= -1.0e+04:

                magz_a = magz_eq[0, 0]

            else:

                liouvillian = compute_liouvillian(
                    pb=pb, kex_ab=kex_ab,
                    lambda_i_a=lambda_i_a, rho_i_a=rho_i_a, omega_i_a=omega_i_a,
                    lambda_i_b=lambda_i_b, rho_i_b=rho_i_a, omega_i_b=omega_i_a + domega_i_ab,
                    omega1x_i=omega1x_i
                )

                s, vr = linalg.eig(liouvillian)
                vri = linalg.inv(vr)

                sl1 = [2, 5]
                sl2 = [i for i, omega_i_a_array in enumerate(s.imag) if abs(omega_i_a_array) < 0.9 * omega1x_i]
                sl3 = [2]

                vri = vri[np.ix_(sl2, sl1)]
                t = np.diag(np.exp(s[sl2] * self.time_t1))
                vr = vr[np.ix_(sl3, sl2)]

                magz_a = np.dot(np.dot(np.dot(vr, t), vri), magz_eq)[0, 0]
                magz_a = magz_a.real

            profile.append(magz_a)

        return np.asarray(profile)

    def calculate_profile(self, params, b1_offsets=None):

        kwargs = {
            short_name: params[long_name].value
            for short_name, long_name in self.map_names.items()
            }

        values = self._calculate_profile(**kwargs)
        scale = self._calculate_scale(values)

        if b1_offsets is not None:
            self.b1_offsets, b1_orig = b1_offsets, self.b1_offsets
            values = self._calculate_profile.__wrapped__(**kwargs)
            self.b1_offsets = b1_orig

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

    def b1_offsets_to_ppm(self, b1_offsets=None):

        if b1_offsets is None:
            b1_offsets = self.b1_offsets

        return two_pi * b1_offsets / self.ppm_to_rads + self.carrier

    def filter_points(self, params):
        """Evaluate some criteria to know whether the point should be considered
        in the calculation or not.

        Returns 'True' if the point should NOT be considered.
        """

        cs = params[self.map_names['cs_i_a']].value
        nu_offsets = (
            (cs - self.carrier) * self.ppm_to_rads / (2.0 * np.pi) - self.b1_offsets
        )

        mask = abs(nu_offsets) > self.on_resonance_filter * 0.5

        self.b1_offsets = self.b1_offsets[mask]
        self.val = self.val[mask]
        self.err = self.err[mask]

    def print_profile(self, params=None):
        """Print the data point"""

        output = []

        if params is not None:
            values = self.calculate_profile(params)
        else:
            values = self.val

        iter_vals = zip(self.b1_offsets, self.val, self.err, values)

        for b1_offset, val, err, cal in iter_vals:

            line = (
                "{0.profile_name:10s} "
                "{0.h_larmor_frq:8.1f} "
                "{0.time_t1:8.1e} "
                "{0.b1_frq:10.1f} "
                "{0.temperature:5.1f} "
                "{1:15.8e} "
                "{2:15.8e} "
                "{3:15.8e} "
                    .format(self, b1_offset, val, err)
            )

            if params is not None:
                line += "{:15.8e}".format(cal)

            output.append(line)

        output.append("")

        return "\n".join(output).upper()
