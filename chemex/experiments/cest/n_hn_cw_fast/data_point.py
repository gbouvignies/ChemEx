from inspect import getargspec
from math import pi

from chemex import parsing
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.constants import xi_ratio
from .back_calculation import make_calc_observable
from ..plotting import plot_data


TWO_PI = 2.0 * pi
RATIO_N = xi_ratio['N']

PAR_DICT = {
    'par_conv': (
        (str, ('resonance_id',)),
        (float, ('h_larmor_frq', 'temperature', 'carrier_h', 'carrier',
                 'time_t1', 'b1_frq', 'b1_frq_h', 'b1_offset')),
        (int, ())
    ),
    'exp': ('resonance_id', 'h_larmor_frq', 'temperature', 'carrier',
            'carrier_h', 'time_t1', 'b1_frq', 'b1_frq_h', 'b1_offset'),
    'fit': ('pb', 'kex', 'dw_h', 'dw_n', 'i0', 'r_nxy', 'dr_nxy', 'r_nz',
            'r_2hxynxy', 'r_hxy', 'etaxy', 'etaz'),
    'fix': ('cs_n', 'cs_h', 'r_2hznz', 'r_hz', 'j_hn'),
}


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'],
                               plot_data)

        self.par['ppm_to_rads_h'] = TWO_PI * self.par['h_larmor_frq']
        self.par['ppm_to_rads'] = TWO_PI * self.par['h_larmor_frq'] * RATIO_N

        temperature = self.par['temperature']
        h_larmor_frq = self.par['h_larmor_frq']
        resonance_id = self.par['resonance_id']
        experiment_name = self.par['experiment_name']

        assignment = parsing.parse_assignment(resonance_id)
        index_1, residue_type_1, nucleus_type_1 = assignment[0]
        index_2, residue_type_2, nucleus_type_2 = assignment[1]
        nucleus_name_1 = residue_type_1 + str(index_1) + nucleus_type_1
        nucleus_name_2 = residue_type_2 + str(index_2) + nucleus_type_2
        nuclei_name = parsing.assignment_name(assignment)

        self.par['_id'] = tuple((temperature, nucleus_name_2, h_larmor_frq))

        args = (
            self.par[arg]
            for arg in getargspec(make_calc_observable.__wrapped__).args
        )

        self.calc_observable = make_calc_observable(*args)

        self.short_long_par_names = (
            ('i0', ('i0', resonance_id, experiment_name)),
            ('pb', ('pb', temperature)),
            ('kex', ('kex', temperature)),
            ('dw_h', ('dw', nucleus_name_2)),
            ('dw_n', ('dw', nucleus_name_1)),
            ('cs_h', ('cs', nucleus_name_2)),
            ('cs_n', ('cs', nucleus_name_1)),
            ('r_nxy', ('r_nxy', nucleus_name_1, h_larmor_frq, temperature)),
            ('dr_nxy', ('dr_nxy', nucleus_name_1, h_larmor_frq, temperature)),
            ('r_nz', ('r_nz', nucleus_name_1, h_larmor_frq, temperature)),
            ('r_hxy', ('r_hxy', nucleus_name_2, h_larmor_frq, temperature)),
            ('r_hz', ('r_hz', nucleus_name_2, h_larmor_frq, temperature)),
            ('r_2hznz', ('r_2hznz', nuclei_name, h_larmor_frq, temperature)),
            (
            'r_2hxynxy', ('r_2hxynxy', nuclei_name, h_larmor_frq, temperature)),
            ('etaxy', ('etaxy', nuclei_name, h_larmor_frq, temperature)),
            ('etaz', ('etaz', nuclei_name, h_larmor_frq, temperature)),
            ('j_hn', ('j_hn', nuclei_name, temperature)),
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

