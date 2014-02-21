'''
Created on Feb 23, 2012

@author: Mike Latham
@author: Guillaume Bouvignies
'''

from inspect import getargspec
from scipy import pi

from chemex.tools import parse_assignment
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.constants import gamma_ratio
from .back_calculation import make_calc_observable
from ..plotting import plot_data

#Constants
TWO_PI = 2.0 * pi
RATIO = gamma_ratio['C']

PAR_DICT = {
    'par_conv': (
        (str, ('resonance_id', 'smallflg')),
        (float, ('h_larmor_frq', 'temperature', 'time_t2')),
        (int, ('ncyc',))
    ),
    'exp': (
        'resonance_id', 'smallflg', 'h_larmor_frq', 'temperature', 'time_t2', 'ncyc'
    ),
    'fit': (
        'i0', 'pb', 'kex', 'dwc', 'dwh', 'r_2hxycxy'
    ),
    'fix': (
        'dr_2hxycxy',
    ),
}


class DataPoint(BaseDataPoint):
    '''Intensity measured during a cpmg pulse train of frequency frq'''

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'], plot_data)

        self.par['ppm_to_rads_h'] = TWO_PI * self.par['h_larmor_frq']
        self.par['ppm_to_rads_c'] = TWO_PI * self.par['h_larmor_frq'] * RATIO

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
            ('kex', ('kex', temperature)),
            ('dwc', ('dw', nucleus_name_1, temperature)),
            ('dwh', ('dw', nucleus_name_2, temperature)),
            ('i0', ('i0', resonance_id, experiment_name)),
            ('r_2hxycxy', ('r_2hxycxy', nucleus_name_1, nucleus_name_2, h_larmor_frq, temperature)),
            ('dr_2hxycxy', ('dr_2hxycxy', nucleus_name_1, nucleus_name_2, h_larmor_frq, temperature)),
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
        '''Print the data point'''

        output = list()
        output.append('{resonance_id:6s}'.format(**self.par))
        output.append('{h_larmor_frq:6.1f}'.format(**self.par))
        output.append('{time_t2:6.1e}'.format(**self.par))
        output.append('{ncyc:4d}'.format(**self.par))
        output.append('{temperature:4.1f}'.format(**self.par))
        output.append('{smallflg:3s}'.format(**self.par))
        output.append('{:10.5f}'.format(self.val))
        output.append('{:10.5f}'.format(self.err))

        if self.cal:
            output.append('{:10.5f}'.format(self.cal))

        return ' '.join(output)
