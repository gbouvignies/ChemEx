"""
Created on Aug 5, 2011

@author: guillaume
"""

from inspect import getargspec
from scipy import pi

from chemex.parsing import parse_assignment
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.constants import xi_ratio
from .back_calculation import make_calc_observable
from ..plotting import plot_data


# Constants
TWO_PI = 2.0 * pi
RATIO_N = xi_ratio['N']

PAR_DICT = {
    'par_conv': (
        (str, ('resonance_id',)),
        (float, (
        'h_larmor_frq', 'temperature', 'carrier', 'time_t2', 'pw', 'taub',
        'time_equil')),
        (int, ('ncyc',))
    ),
    'exp': (
        'resonance_id', 'h_larmor_frq', 'temperature', 'carrier',
        'time_t2', 'time_equil', 'pw', 'taub', 'ncyc',
    ),
    'fit': (
        'pb', 'pc', 'kex_ab', 'kex_bc', 'dw_ab', 'dw_ac', 'i0', 'r_nxy'
    ),
    'fix': (
        'kex_ac', 'r_nz', 'dr_nxy_ab', 'dr_nxy_ac', 'r_2hznz',
        'etaxy', 'etaz', 'cs', 'j_hn', 'dj_hn_ab', 'dj_hn_ac',
    ),
}


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'],
                               plot_data)

        self.par['ppm_to_rads'] = TWO_PI * self.par['h_larmor_frq'] * RATIO_N

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']
        h_larmor_frq = self.par['h_larmor_frq']
        experiment_name = self.par['experiment_name']

        assignment = parse_assignment(resonance_id)
        index_1, residue_type_1, nucleus_type_1 = assignment[0]
        index_2, residue_type_2, nucleus_type_2 = assignment[1]
        nucleus_name_1 = residue_type_1 + str(index_1) + nucleus_type_1
        nucleus_name_2 = residue_type_2 + str(index_2) + nucleus_type_2

        self.par['_id'] = tuple((temperature, nucleus_name_1, h_larmor_frq))

        args = (
            self.par[arg]
            for arg in getargspec(make_calc_observable.__wrapped__).args
        )

        self.calc_observable = make_calc_observable(*args)

        self.kwargs_default = {'ncyc': self.par['ncyc']}

        self.short_long_par_names = (
            ('pb', ('pb', temperature)),
            ('pc', ('pc', temperature)),
            ('kex_ab', ('kex_ab', temperature)),
            ('kex_bc', ('kex_bc', temperature)),
            ('kex_ac', ('kex_ac', temperature)),
            ('dw_ab', ('dw_ab', nucleus_name_1)),
            ('dw_ac', ('dw_ac', nucleus_name_1)),
            ('cs', ('cs', nucleus_name_1, temperature)),
            ('i0', ('i0', resonance_id, experiment_name)),
            ('r_nxy', ('r_nxy', nucleus_name_1, h_larmor_frq, temperature)),
            ('dr_nxy_ab',
             ('dr_nxy_ab', nucleus_name_1, h_larmor_frq, temperature)),
            ('dr_nxy_ac',
             ('dr_nxy_ac', nucleus_name_1, h_larmor_frq, temperature)),
            ('r_nz', ('r_nz', nucleus_name_1, h_larmor_frq, temperature)),
            ('r_2hznz', (
            'r_2hznz', nucleus_name_1, nucleus_name_2, h_larmor_frq,
            temperature)),
            ('etaxy', ('etaxy', nucleus_name_1, h_larmor_frq, temperature)),
            ('etaz', ('etaz', nucleus_name_1, h_larmor_frq, temperature)),
            ('j_hn', ('j_hn', nucleus_name_1, nucleus_name_2, temperature)),
            ('dj_hn_ab',
             ('dj_hn_ab', nucleus_name_1, nucleus_name_2, temperature)),
            ('dj_hn_ac',
             ('dj_hn_ac', nucleus_name_1, nucleus_name_2, temperature)),
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
        output.append('{resonance_id:6s}'.format(**self.par))
        output.append('{h_larmor_frq:6.1f}'.format(**self.par))
        output.append('{time_t2:6.1e}'.format(**self.par))
        output.append('{ncyc:4d}'.format(**self.par))
        output.append('{temperature:4.1f}'.format(**self.par))
        output.append('{:8.5f}'.format(self.val))
        output.append('{:8.5f}'.format(self.err))

        if self.cal:
            output.append('{:8.3f}'.format(self.cal))

        return ' '.join(output)

