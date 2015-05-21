from inspect import getargspec
from math import pi

from chemex.parsing import parse_assignment
from chemex.experiments.base_profile import BaseDataPoint, get_par
from chemex.constants import xi_ratio
from .back_calculation import make_calc_observable
from chemex.experiments.util import calc_multiplet
from ..plotting import plot_data








# Constants
RATIO_N = xi_ratio['N']
TWO_PI = 2.0 * pi
PAR_DICT = {
    'par_conv': ((str, ('resonance_id',)),
                 (float, (
                 'h_larmor_frq', 'temperature', 'carrier', 'time_t1', 'b1_frq',
                 'b1_offset',)),
                 (int, ())),
    'exp': (
    'resonance_id', 'h_larmor_frq', 'temperature', 'carrier', 'time_t1',
    'b1_frq', 'b1_offset',),
    'fit': ('pb', 'kex', 'dw', 'i0', 'r_nxy', 'dr_nxy', 'r_nz'),
    'fix': ('cs',),
}

J_COUPLINGS = (7.7, 10.7, 14.4)


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'],
                               plot_data)

        self.par['ppm_to_rads'] = TWO_PI * self.par['h_larmor_frq'] * RATIO_N
        self.par['multiplet'] = calc_multiplet(J_COUPLINGS)

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']
        h_larmor_frq = self.par['h_larmor_frq']
        experiment_name = self.par['experiment_name']

        assignment = parse_assignment(resonance_id)
        index, residue_type, nucleus_type = assignment[0]
        nucleus_name = ''.join([residue_type, str(index), nucleus_type])

        self.par['_id'] = tuple((temperature, nucleus_name, h_larmor_frq))

        args = (self.par[arg] for arg in
                getargspec(make_calc_observable.__wrapped__).args)
        self.calc_observable = make_calc_observable(*args)

        self.short_long_par_names = (
            ('pb', ('pb', temperature)),
            ('kex', ('kex', temperature)),
            ('dw', ('dw', nucleus_name)),
            ('cs', ('cs', nucleus_name, temperature)),
            ('i0', ('i0', resonance_id, experiment_name)),
            ('r_nxy', ('r_nxy', nucleus_name, h_larmor_frq, temperature)),
            ('dr_nxy', ('dr_nxy', nucleus_name, h_larmor_frq, temperature)),
            ('r_nz', ('r_nz', nucleus_name, h_larmor_frq, temperature)),
        )

        self.fitting_parameter_names.update(
            long_name
            for short_name, long_name in self.short_long_par_names
            if short_name in PAR_DICT['fit']
        )

        self.fixed_parameter_names.update(
            long_name
            for short_name, long_name in self.short_long_par_names
            if short_name in PAR_DICT['fix']
        )

    def __repr__(self):
        """Print the data point"""

        output = list()
        output.append('{resonance_id:10s}'.format(**self.par))
        output.append('{h_larmor_frq:6.1f}'.format(**self.par))
        output.append('{time_t1:6.1e}'.format(**self.par))
        output.append('{b1_offset:9.3e}'.format(**self.par))
        output.append('{b1_frq:6.1e}'.format(**self.par))
        output.append('{temperature:4.1f}'.format(**self.par))
        output.append('{:8.5f}'.format(self.val))
        output.append('{:8.5f}'.format(self.err))

        if self.cal:
            output.append('{:8.5f}'.format(self.cal))

        return ' '.join(output)

    def update_b1_offset(self, b1_offset):
        """Update b1_offset value"""

        self.par['b1_offset'] = b1_offset
        args = (self.par[arg] for arg in
                getargspec(make_calc_observable.__wrapped__).args)
        self.calc_observable = make_calc_observable(*args)

    def filter(self, par, par_indexes, par_fixed=None):
        filter_range = float(self.par.get('on_resonance_filter', 0.0))

        par_val = dict(
            (short_name, get_par(long_name, par, par_indexes, par_fixed))
            for short_name, long_name in self.short_long_par_names)

        cs = par_val['cs']
        cs_offset_hz = (cs - self.par['carrier']) * self.par['ppm_to_rads'] / (
        2.0 * pi) - self.par['b1_offset']

        return abs(cs_offset_hz) < filter_range * 0.5

