"""
Created on Aug 5, 2011

@author: guillaume
"""

from inspect import getargspec
from scipy import pi

from chemex.tools import parse_assignment
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.constants import xi_ratio
from .back_calculation import make_calc_observable
from ..plotting import plot_data


TWO_PI = 2.0 * pi
RATIO_C = xi_ratio['C']

PAR_DICT = {

    # experimental requirements (used in __init__ and check_parameters)
    'par_conv': (
        (str, (
            'resonance_id',
        )),
        (float, (
            'h_larmor_frq',
            'temperature',
            'carrier',
            'time_t2',
            'pw',
            'time_equil',
        )),
        (int, (
            'ncyc',
        ))
    ),
    'exp': (
        'resonance_id',
        'h_larmor_frq',
        'temperature',
        'carrier',
        'time_t2',
        'time_equil',
        'pw',
        'ncyc',
    ),
    'fit': (
        'pb',
        'kex',
        'dw',
        'i0',
        'r_cxy',
    ),
    'fix': (
        'r_cz',
        'dr_cxy',
        'cs',
    ),

}


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'], plot_data)

        self.par['ppm_to_rads'] = TWO_PI * self.par['h_larmor_frq'] * RATIO_C

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']
        h_larmor_frq = self.par['h_larmor_frq']
        experiment_name = self.par['experiment_name']

        assignment = parse_assignment(resonance_id)
        index, residue_type, nucleus_type = assignment[0]
        nucleus_name = residue_type + str(index) + nucleus_type

        self.par['_id'] = tuple((temperature, nucleus_name, h_larmor_frq))

        args = (self.par[arg] for arg in getargspec(make_calc_observable.__wrapped__).args)
        self.calc_observable = make_calc_observable(*args)

        self.kwargs_default = {'ncyc': self.par['ncyc']}

        self.short_long_par_names = (
            ('i0', ('i0', resonance_id, experiment_name)),
            ('pb', ('pb', temperature)),
            ('kex', ('kex', temperature)),
            ('dw', ('dw', nucleus_name)),
            ('cs', ('cs', nucleus_name, temperature)),
            ('r_cxy', ('r_cxy', nucleus_name, h_larmor_frq, temperature)),
            ('dr_cxy', ('dr_cxy', nucleus_name, h_larmor_frq, temperature)),
            ('r_cz', ('r_cz', nucleus_name, h_larmor_frq, temperature)),
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
        output.append('{:10.5f}'.format(self.val))
        output.append('{:10.5f}'.format(self.err))

        if self.cal:
            output.append('{:10.5f}'.format(self.cal))

        return ' '.join(output)

