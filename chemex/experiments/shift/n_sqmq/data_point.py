from scipy import pi

from chemex.parsing import parse_assignment
from chemex.experiments.base_data_point import BaseDataPoint
from back_calculation import calc_observable
from chemex.constants import xi_ratio
from ..plotting import plot_data


TWO_PI = 2.0 * pi
RATIO = xi_ratio['N']

PAR_DICT = {
    'par_conv': (
        (str, ('resonance_id',)),
        (float, ('h_larmor_frq', 'temperature',)),
        (int, ())
    ),
    'exp': ('resonance_id', 'h_larmor_frq', 'temperature',),
    'fit': (),
    'fix': ('pb', 'kex', 'dw_h', 'dw_n',),
}


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'],
                               plot_data)

        self.par['ppm_to_rads_h'] = TWO_PI * self.par['h_larmor_frq']
        self.par['ppm_to_rads_n'] = TWO_PI * self.par['h_larmor_frq'] * RATIO

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']

        assignment = parse_assignment(resonance_id)
        index_1, residue_type_1, nucleus_type_1 = assignment[0]
        index_2, residue_type_2, nucleus_type_2 = assignment[1]
        nucleus_name_1 = residue_type_1 + str(index_1) + nucleus_type_1
        nucleus_name_2 = residue_type_2 + str(index_2) + nucleus_type_2

        self.calc_observable = calc_observable

        self.kwargs_default = {
            'ppm_to_rads_h': self.par['ppm_to_rads_h'],
            'ppm_to_rads_n': self.par['ppm_to_rads_n'],
        }

        self.short_long_par_names = (
            ('pb', ('pb', temperature)),
            ('kex', ('kex', temperature)),
            ('dw_n', ('dw', nucleus_name_1)),
            ('dw_h', ('dw', nucleus_name_2)),
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
        output.append('{temperature:4.1f}'.format(**self.par))
        output.append('{:10.5f}'.format(self.val))
        output.append('{:10.5f}'.format(self.err))

        if self.cal:
            output.append('{:10.5f}'.format(self.cal))

        return ' '.join(output)
    

