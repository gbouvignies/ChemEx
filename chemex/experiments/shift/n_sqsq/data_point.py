from scipy import pi

from ....parsing import parse_assignment
from ....experiments.base_data_point import BaseDataPoint
from ....constants import xi_ratio
from .back_calculation import calc_observable
from ..plotting import plot_data


TWO_PI = 2.0 * pi
RATIO = xi_ratio['N']

PAR_DICT = {
    'par_conv': (
        (str, ('resonance_id',)),
        (float, ('h_larmor_frq_1', 'h_larmor_frq_2', 'temperature',)),
        (int, ())
    ),
    'exp': ('resonance_id', 'h_larmor_frq_1', 'h_larmor_frq_2', 'temperature',),
    'fit': (),
    'fix': ('pb', 'kex', 'dw_n',),
}


class DataPoint(BaseDataPoint):
    '''Intensity measured during a cpmg pulse train of frequency frq'''

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'],
                               plot_data)

        self.par['ppm_to_rads_n_1'] = \
            TWO_PI * self.par['h_larmor_frq_1'] * RATIO
        self.par['ppm_to_rads_n_2'] = \
            TWO_PI * self.par['h_larmor_frq_2'] * RATIO

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']

        assignment = parse_assignment(resonance_id)
        index, residue_type, nucleus_type = assignment[0]
        nucleus_name = ''.join([residue_type, str(index), nucleus_type])

        self.calc_observable = calc_observable

        self.kwargs_default = {
            'ppm_to_rads_n_1': self.par['ppm_to_rads_n_1'],
            'ppm_to_rads_n_2': self.par['ppm_to_rads_n_2'],
        }

        self.short_long_par_names = (
            ('pb', ('pb', temperature)),
            ('kex', ('kex', temperature)),
            ('dw_n', ('dw', nucleus_name)),
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
        output.append('{h_larmor_frq_1:6.1f}'.format(**self.par))
        output.append('{h_larmor_frq_2:6.1f}'.format(**self.par))
        output.append('{temperature:4.1f}'.format(**self.par))
        output.append('{:10.5f}'.format(self.val))
        output.append('{:10.5f}'.format(self.err))

        if self.cal:
            output.append('{:10.5f}'.format(self.cal))

        return ' '.join(output)
    

