import scipy as sp
import scipy.linalg as la
import lmfit as lf

from chemex import constants
from chemex import parsing
from chemex.experiments import abc_profile
from chemex.experiments import utils
from chemex.caching import lru_cache
from chemex.bases.two_states import single_spin
from ..plotting import plot_data


inv = la.inv
eig = la.eig
dot = sp.dot
diag = sp.diag
exp = sp.exp
calculate_shift_ex_2st = utils.calculate_shift_ex_2st
compute_liouvillian = single_spin.compute_liouvillian
check_par = abc_profile.check_par

RATIO_N = constants.xi_ratio['N']
TWO_PI = 2.0 * sp.pi


class Profile(abc_profile.ABCProfile):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):

        self.val = sp.array(val)
        self.err = sp.array(err)

        self.cal = [None] * len(self.val)
        self.scale = 0.0

        self.experiment_name = check_par(par, 'experiment_name')
        self.h_larmor_frq = check_par(par, 'h_larmor_frq', float)
        self.temperature = check_par(par, 'temperature', float)
        self.carrier = check_par(par, 'carrier', float)
        self.b1_frq = check_par(par, 'b1_frq', float)
        self.time_t1 = check_par(par, 'time_t1', float)
        self.resonance_id = check_par(par, 'resonance_id')
        self.b1_offsets = check_par(par, 'b1_offsets', sp.array)
        self.experiment_name = check_par(par, 'experiment_name')

        self.ppm_to_rads = self.h_larmor_frq * RATIO_N * TWO_PI

        self.calc_profile = lru_cache(maxsize=5)(self.calc_profile)
        self.plot_data = plot_data

        self.assignment = parsing.parse_assignment(self.resonance_id)
        self.nucleus_name = parsing.assignment_name(self.assignment[:1])

        self.map_names = {
            'pb': '__'.join([
                'pb',
                str(self.temperature),
            ]),
            'kex': '__'.join([
                'kex',
                str(self.temperature),
            ]),
            'cs': '__'.join([
                'cs',
                self.nucleus_name,
                str(self.temperature),
            ]),
            'r_nxy': '__'.join([
                'r_nxy',
                self.nucleus_name,
                str(self.h_larmor_frq),
                str(self.temperature),
            ]),
            'dr_nxy': '__'.join([
                'dr_nxy',
                self.nucleus_name,
                str(self.h_larmor_frq),
                str(self.temperature),
            ]),
            'r_nz': '__'.join([
                'r_nz',
                self.nucleus_name,
                str(self.h_larmor_frq),
                str(self.temperature),
            ])
        }

    def make_default_parameters(self):

        parameters = lf.Parameters()

        parameters.add_many(
            # name, value, vary, min, max, expr
            (self.map_names['pb'], 0.05, True, 0.0, 1.0, None),
            (self.map_names['kex'], 100.0, True, 0.0, None, None),
            (self.map_names['cs'], 0.0, False, None, None, None),
            (self.map_names['r_nxy'], 10.0, True, 0.0, 1.0, None),
            (self.map_names['dr_nxy'], 0.0, True, None, None, None),
            (self.map_names['r_nz'], 1.0, True, 0.0, None, None),
        )

        return parameters

    def calculate_profile(self, pb=0.0, kex=0.0, dw=0.0, r_nz=1.5, r_nxy=0.0,
                          dr_nxy=0.0, cs=0.0):
        """Calculate the intensity in presence of exchange after a CEST block
        assuming initial intensity of 1.0.

        Parameters
        ----------
        pb : float
            Fractional population of state B,
            0.0 for 0%, 1.0 for 100%
        kex : float
            Exchange rate between state A and B in /s.
        dw : float
            Chemical shift difference between states A and B in rad/s.
        r_nz : float
            Longitudinal relaxation rate of state {a,b} in /s.
        r_nxy : float
            Transverse relaxation rate of state a in /s.
        dr_nxy : float
            Transverse relaxation rate difference between states a and b in /s.
        cs : float
            Resonance position in rad/s.

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
            r_ixy=r_nxy,
            dr_ixy=dr_nxy
        )

        w = (
            (cs - self.carrier) * self.ppm_to_rads -
            shift_ex -
            TWO_PI * self.b1_offsets
        )

        w1 = TWO_PI * self.b1_frq

        magz_eq = sp.array([[1 - pb, pb]]).transpose()

        magz_a_list = []

        for b1_offset, cs_offset in zip(self.b1_offsets, w):

            if b1_offset <= -1.0e+04:

                magz_a = magz_eq[0, 0]

            else:

                liouvillian = compute_liouvillian(
                    pb=pb,
                    kex=kex,
                    dw=dw_,
                    r_ixy=r_nxy,
                    dr_ixy=dr_nxy,
                    r_iz=r_nz,
                    w=cs_offset,
                    w1x=w1
                )

                s, vr = eig(liouvillian)
                vri = inv(vr)

                sl1 = [2, 5]
                sl2 = [i for i, w in enumerate(s.imag) if abs(w) < 1.0e-6]
                sl3 = [2]

                vri = vri[sp.ix_(sl2, sl1)].real
                t = diag(exp(s[sl2].real * self.time_t1))
                vr = vr[sp.ix_(sl3, sl2)].real

                magz_a = dot(dot(dot(vr, t), vri), magz_eq)[0, 0]

            magz_a_list.append(magz_a)

        return sp.asarray(magz_a_list)

    def back_calculate(self, parameters=None, scaling=True):

        kwargs = {
            short_name: parameters[long_name].value
            for short_name, long_name in self.map_names.items()
        }

        values = self.calc_profile(**kwargs)

        if scaling:
            values *= self.calculate_scale(values)

        self.cal = values

        return values

    def calculate_scale(self, cal):

        scale = (
            sum(cal * self.val / self.err ** 2) /
            sum((cal / self.err) ** 2)
        )

        self.scale = scale

        return scale

    def calculate_residuals(self, parameters=None):
        """Calculates the residual between the experimental and
        back-calculated values.
        """

        values = self.back_calculate(parameters)

        return (self.val - values) / self.err

    def b1_offsets_to_ppm(self, b1_offsets=None):

        if b1_offsets is None:
            b1_offsets = self.b1_offsets

        return TWO_PI * b1_offsets / self.ppm_to_rads + self.carrier

    def filter_points(self, parameters=None):
        """Evaluate some criteria to know whether the point should be considered
        in the calculation or not.

        Returns 'True' if the point should NOT be considered.
        """

        return False

    def print_profile(self):
        """Print the data point"""

        output = []

        iter_vals = zip(self.b1_offsets, self.val, self.err, self.cal)

        for b1_offset, val, err, cal in iter_vals:

            line = (
                "{0.resonance_id:10s} "
                "{0.h_larmor_frq:8.1f} "
                "{0.time_t1:8.1e} "
                "{0.b1_frq:10.1f} "
                "{0.temperature:5.1f} "
                "{1:15.8e} "
                "{2:15.8e} "
                "{3:15.8e} "
                .format(self, b1_offset, val, err)
            )

            if cal is not None:
                line += "{:15.8e}".format(cal * self.scale)

            output.append(line)

        return "\n".join(output)


