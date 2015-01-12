"""
Created on Aug 5, 2011

@author: guillaume
"""

# Standard imports
from inspect import getargspec
from math import pi

# Local imports
from chemex.parsing import parse_assignment
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.constants import xi_ratio
from chemex.experiments.misc import calc_multiplet
from .back_calculation import make_calc_observable
from ..plotting import plot_data

# Constants
RATIO_C = xi_ratio['C']
TWO_PI = 2.0 * pi

# first dictionary version
PAR_DICT = {
    'par_conv': (
        (str, (
            'resonance_id',
        )),
        (float, (
            'h_larmor_frq',
            'temperature',
            'carrier',
            'time_t1',
            'b1_frq',
            'b1_offset',
            'b1_inh',
        )),
        (int, (
            'b1_inh_res',
        ))
    ),
    'exp': (
        'resonance_id',
        'h_larmor_frq',
        'temperature',
        'carrier',
        'time_t1',
        'b1_frq',
        'b1_offset',
        'b1_inh',
        'b1_inh_res',
    ),
    'fit': (
        'pb',
        'kex',
        'dw',
        'i0',
        'r_cxy',
        'dr_cxy',
        'r_cz'
    ),
    'fix': (
        'cs',
    ),

}

# This is the dictionary that contains scalar coupling values affecting each nucleus
J_COUPLINGS = {
    'a': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0,),
    },
    'c': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0,),
    },
    'd': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 51.0,),
    },
    'e': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0, 51.0,),
    },
    'f': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 51.0,),
    },
    'g': {
        'c': (52.0, 14.4,),
        'ca': (52.0,),
    },
    'h': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 51.0,),
        'cg': (51.0, 72.0,),
        'cd2': (72.0,),
        'ce1': (),
    },
    'i': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0, 35.0,),
        'cg1': (35.0, 35.0,),
        'cg2': (35.0,),
        'cd1': (35.0,),
    },
    'k': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0, 35.0,),
        'cd': (35.0, 35.0,),
        'ce': (35.0,),
    },
    'l': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0, 35.0, 35.0,),
        'cd1': (35.0,),
        'cd2': (35.0,),
    },
    'm': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0,),
        'ce': (),
    },
    'n': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 51.0,),
    },
    'p': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0, 35.0,),
        'cd': (35.0,),
    },
    'q': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0, 51.0,),
    },
    'r': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg': (35.0, 35.0,),
        'cd': (35.0,),
    },
    's': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0,),
    },
    't': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0,),
        'cg2': (35.0,),
    },
    'v': {
        'ca': (52.0, 35.0,),
        'cb': (35.0, 35.0, 35.0,),
        'cg1': (35.0,),
        'cg2': (35.0,),
    },
    'w': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 51.0,),
    },
    'y': {
        'c': (52.0, 14.4,),
        'ca': (52.0, 35.0,),
        'cb': (35.0, 51.0,),
    },
    '': {
        'c': (52.0, 14.4,),
    },
}


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):
        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'], plot_data)

        self.par['ppm_to_rads'] = TWO_PI * self.par['h_larmor_frq'] * RATIO_C

        self.kwargs_default = dict()

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']
        h_larmor_frq = self.par['h_larmor_frq']
        experiment_name = self.par['experiment_name']

        assignment = parse_assignment(resonance_id)
        index, residue_type, nucleus_type = assignment[0]
        nucleus_name = residue_type + str(index) + nucleus_type

        couplings = J_COUPLINGS[residue_type][nucleus_type]
        self.par['multiplet'] = calc_multiplet(couplings)

        self.par['_id'] = tuple((temperature, nucleus_name, h_larmor_frq))

        args = (self.par[arg] for arg in getargspec(make_calc_observable.__wrapped__).args)
        self.calc_observable = make_calc_observable(*args)

        self.short_long_par_names = (
            ('pb', ('pb', temperature)),
            ('kex', ('kex', temperature)),
            ('dw', ('dw', nucleus_name)),
            ('cs', ('cs', nucleus_name, temperature)),
            ('i0', ('i0', resonance_id, experiment_name)),
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
        """update b1_offset value"""

        self.par['b1_offset'] = b1_offset
        args = (self.par[arg] for arg in getargspec(make_calc_observable.__wrapped__).args)
        self.calc_observable = make_calc_observable(*args)

