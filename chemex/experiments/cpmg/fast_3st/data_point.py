"""
Created on Aug 5, 2011

@author: guillaume
"""

from inspect import getargspec

from scipy import pi

from chemex.constants import gamma_ratio
from chemex.experiments.base_data_point import BaseDataPoint
from chemex.tools import parse_assignment
from ..plotting import plot_data
from .back_calculation import make_calc_observable


PAR_DICT = {
    # experimental requirements (used in __init__ and check_parameters)
    'par_conv': (
        (str, ('resonance_id',)),
        (float, ('h_larmor_frq', 'temperature', 'carrier', 'time_t2',)),
        (int, ('ncyc',))
    ),
    # Some stuff to get a nice help output
    'exp': ('resonance_id', 'h_larmor_frq', 'temperature', 'time_t2', 'ncyc',),
    'fit': ('pb', 'pc', 'kex_ab', 'kex_bc', 'kex_ac', 'dw_ab', 'dw_ac', 'i0', 'r_ixy',),
    'fix': ('dr_ixy_ab', 'dr_ixy_ac',),
}


class DataPoint(BaseDataPoint):
    """Intensity measured during a cpmg pulse train of frequency frq"""

    def __init__(self, val, err, par):

        BaseDataPoint.__init__(self, val, err, par, PAR_DICT['par_conv'], plot_data)

        temperature = self.par['temperature']
        resonance_id = self.par['resonance_id']
        h_larmor_frq = self.par['h_larmor_frq']
        experiment_name = self.par['experiment_name']

        index, residue_type, nucleus_type = parse_assignment(resonance_id)[0]
        nucleus_name = residue_type + str(index) + nucleus_type

        try:
            self.par['ppm_to_rads'] = (
                2.0 * pi * self.par['h_larmor_frq'] * gamma_ratio[nucleus_type[0].upper()]
            )
        except KeyError:
            exit(
                "Unknown nucleus type \"{}\" for peak \"{}\" in experiment \"{}\""
                .format(nucleus_type, resonance_id, experiment_name)
            )

        self.par['_id'] = ((temperature, nucleus_name, h_larmor_frq),)

        args = (
            self.par[arg]
            for arg in getargspec(make_calc_observable.__wrapped__).args
        )
        self.calc_observable = make_calc_observable(*args)

        self.kwargs_default = {'ncyc': self.par['ncyc']}

        self.short_long_par_names = (
            ('i0', ('i0', resonance_id, experiment_name)),
            ('pb', ('pb', temperature)),
            ('pc', ('pc', temperature)),
            ('kex_ab', ('kex_ab', temperature)),
            ('kex_bc', ('kex_bc', temperature)),
            ('kex_ac', ('kex_ac', temperature)),
            ('dw_ab', ('dw_ab', nucleus_name)),
            ('dw_ac', ('dw_ac', nucleus_name)),
            ('r_ixy', ('r_ixy', nucleus_name, h_larmor_frq, temperature)),
            ('dr_ixy_ab', ('dr_ixy_ab', nucleus_name, h_larmor_frq, temperature)),
            ('dr_ixy_ac', ('dr_ixy_ac', nucleus_name, h_larmor_frq, temperature)),
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



