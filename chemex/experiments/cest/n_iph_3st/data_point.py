"""
Created on Sept 20, 2013

@author: guillaume
"""

# Standard imports
from inspect import getargspec
from math import pi

# Local imports
from chemex.tools import parse_assignment
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.constants import gamma_ratio
from .back_calculation import make_calc_observable
from ..plotting import plot_data

# Constants
RATIO_N = gamma_ratio['N']
TWO_PI = 2.0 * pi

# first dictionary version
PAR_DICT = {

    # experimental requirements (used in __init__ and check_parameters)
    'par_conv'  :((str  , ('resonance_id',)),

                  (float, ('h_larmor_frq',
                           'temperature',
                           'carrier',
                           'time_t1',
                           'B1_frq',
                           'B1_offset',
                           'B1_inh',)),

                  (int  , ('B1_inh_res',))),

    # Some stuff to get a nice help output
    'exp' : ('resonance_id', 'h_larmor_frq', 'temperature', 'carrier',
             'time_t1', 'B1_frq', 'B1_offset', 'B1_inh', 'B1_inh_res'),

    'fit' : ('pb', 'pc', 'kex_ab', 'kex_bc', 'dw_ab', 'dw_ac',
             'I0', 'r_Nxy', 'dr_Nxy_ab', 'dr_Nxy_ac', 'r_Nz'),

    'fix' : ('cs', 'kex_ac',),

}

class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):

        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'], plot_data)

        self.par['B1_offset'] *= -1.0
        self.par['ppm_to_rads'] = TWO_PI * self.par['h_larmor_frq'] * RATIO_N

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']
        h_larmor_frq = self.par['h_larmor_frq']
        experiment_name = self.par['experiment_name']

        assignment = parse_assignment(resonance_id)
        index, residue_type, nucleus_type = assignment[0]
        nucleus_name = residue_type + str(index) + nucleus_type

        self.par['_id'] = tuple((temperature, nucleus_name, h_larmor_frq))

        args = (self.par[arg] for arg in getargspec(make_calc_observable.__wrapped__).args)  # @UndefinedVariable
        self.calc_observable = make_calc_observable(*args)

        self.short_long_par_names = (
            ('pb'        , ('pb', temperature)),
            ('pc'        , ('pc', temperature)),
            ('kex_ab'    , ('kex_ab', temperature)),
            ('kex_bc'    , ('kex_bc', temperature)),
            ('kex_ac'    , ('kex_ac', temperature)),
            ('dw_ab'     , ('dw_ab', nucleus_name)),
            ('dw_ac'     , ('dw_ac', nucleus_name)),
            ('cs'        , ('cs', nucleus_name, temperature)),
            ('I0'        , ('I0', resonance_id, experiment_name)),
            ('r_Nxy'     , ('r_Nxy', nucleus_name, h_larmor_frq, temperature)),
            ('dr_Nxy_ab' , ('dr_Nxy_ab', nucleus_name, h_larmor_frq, temperature)),
            ('dr_Nxy_ac' , ('dr_Nxy_ac', nucleus_name, h_larmor_frq, temperature)),
            ('r_Nz'      , ('r_Nz', nucleus_name, h_larmor_frq, temperature)),
        )

        self.fitting_parameter_names.update(long_name
                                            for short_name, long_name in self.short_long_par_names
                                            if short_name in PAR_DICT['fit'])

        self.fixed_parameter_names.update(long_name
                                          for short_name, long_name in self.short_long_par_names
                                          if short_name in PAR_DICT['fix'])


    def __repr__(self):
        """Print the data point"""

        output = list()
        output.append('{resonance_id:10s}'.format(**self.par))
        output.append('{h_larmor_frq:6.1f}'.format(**self.par))
        output.append('{time_t1:6.1e}'.format(**self.par))
        output.append('{B1_offset:9.3e}'.format(**self.par))
        output.append('{B1_frq:6.1e}'.format(**self.par))
        output.append('{temperature:4.1f}'.format(**self.par))
        output.append('{:8.5f}'.format(self.val))
        output.append('{:8.5f}'.format(self.err))

        if self.cal:
            output.append('{:8.5f}'.format(self.cal))

        return ' '.join(output)

    def update_B1_offset(self, B1_offset):
        """Update B1_offset value"""

        self.par['B1_offset'] = B1_offset
        args = (self.par[arg] for arg in getargspec(make_calc_observable.__wrapped__).args)  # @UndefinedVariable
        self.calc_observable = make_calc_observable(*args)

