"""Pure In-phase CEST (2-state, fast)

Analyzes 15N chemical exchange in the presence of 1H composite decoupling
during the CEST block. This keeps the spin system purely in-phase throughout,
and is calculated using the 6x6, single spin matrix:

[ Nx(a), Ny(a), Nz(a), Nx(b), Ny(b), Nz(b) ]

Notes
-----
The calculation is designed specifically to analyze the experiment found in
the reference:

J Am Chem Soc (2012), 134, 8148-61

"""

import lmfit
import numpy as np
from scipy import linalg

from chemex import constants
from chemex.bases.two_states import single_spin
from chemex.experiments import base_profile, util, sputil
from chemex.experiments.cest import plotting

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

calculate_shift_ex_2st = util.calculate_shift_ex_2st
compute_liouvillian = single_spin.compute_liouvillian
check_par = base_profile.check_par

two_pi = 2.0 * np.pi

attributes_exp = {
    'h_larmor_frq': float,
    'temperature': float,
    'carrier': float,
    'b1_frq': float,
    'time_t1': float,
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
        'vary': True,
        'min': None,
        'max': None,
        'expr': None,
    },
    'r_iz': {
        'value': 1.0,
        'vary': True,
        'min': 0.0,
        'max': None,
        'expr': None,
    },
}


class Profile(base_profile.BaseProfile):
    def __init__(self, profile_name, measurements, exp_details):

        self.profile_name = profile_name
        self.b1_offsets = measurements['b1_offsets']
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

        shift_ex, _ = calculate_shift_ex_2st(
            pb=pb,
            kex=kex,
            dw=dw_,
            r_ixy=r_ixy,
            dr_ixy=dr_ixy
        )

        cs_offsets = (
            (cs - self.carrier) * self.ppm_to_rads -
            shift_ex - two_pi * self.b1_offsets
        )

        w1 = two_pi * self.b1_frq

        magz_eq = np.array([[1 - pb], [pb]])

        profile = []

        for b1_offset, cs_offset in zip(self.b1_offsets, cs_offsets):

            if b1_offset <= -1.0e+04:

                magz_a = magz_eq[0, 0]

            else:

                liouvillian = compute_liouvillian(
                    pb=pb,
                    kex=kex,
                    dw=dw_,
                    r_ixy=r_ixy,
                    dr_ixy=dr_ixy,
                    r_iz=r_iz,
                    w=cs_offset,
                    w1x=w1
                )

                s, vr = linalg.eig(liouvillian)
                vri = linalg.inv(vr)

                sl1 = [2, 5]
                sl2 = [i for i, cs_offsets in enumerate(s.imag) if abs(cs_offsets) < 0.1 * w1]
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

        return False

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
