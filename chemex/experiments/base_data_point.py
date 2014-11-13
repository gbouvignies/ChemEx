"""
Created on 2012-02-21

@author: guillaume
"""


class BaseDataPoint(object):
    """Base class defining an experimental point."""

    def __init__(self, val=0.0, err=0.0, par=None, par_conv=None, plot_data=None, calc_observable=None):
        """Constructor"""

        self.val = float(val)
        self.err = float(err)
        self.par = par.copy()
        self.cal = None
        self.short_long_par_names = None
        self.fitting_parameter_names = set()
        self.fixed_parameter_names = set()
        self.kwargs_default = dict()
        self.calc_observable = calc_observable
        self.plot_data = plot_data

        self.check_parameters(par_conv)


    def __repr__(self):
        """Prints the data point"""

        output = list()
        output.append('{:10.5f}'.format(self.val))
        output.append('{:10.5f}'.format(self.err))

        if self.cal:
            output.append('{:10.5f}'.format(self.cal))

        return ' '.join(output)

    def check_parameters(self, par_conv):
        """Checks that model parameters are provided and convert them to the right type"""

        missing_par = set()

        for func, par_names in par_conv:
            for par_name in par_names:
                if par_name in self.par:
                    self.par[par_name] = func(self.par[par_name])
                else:
                    missing_par.add(self.par['resonance_id'] + " : " + par_name)

        if missing_par:
            msg = 'Missing parameters detected for {:s}'.format(', '.join(missing_par))
            exit(msg)

        return None

    def calc_val(self, par, par_indexes, par_fixed=None):

        kwargs = dict((short_name, get_par(long_name, par, par_indexes, par_fixed))
                      for short_name, long_name in self.short_long_par_names)

        kwargs.update(self.kwargs_default)

        self.cal = self.calc_observable(**kwargs)

    def calc_residual(self, par, par_indexes, par_fixed=None):
        """Calculates the residual between the experimental and back-calculated values."""

        self.calc_val(par, par_indexes, par_fixed)

        return (self.val - self.cal) / self.err

    def get_fitting_parameter_names(self):
        """Provide the parameters that are needed to back-calculate the experimental value."""

        return self.fitting_parameter_names

    def get_fixed_parameter_names(self):
        """Provide the parameters that are needed to back-calculate the experimental value."""

        return self.fixed_parameter_names

    def filter(self, par, par_indexes, par_fixed=None):
        """
        Evaluate some criteria to know whether the point
        should be considered in the calculation or not.

        Returns 'True' if the point should NOT be considered.
        """

        return False


# Functions

def get_par(par_name, par, par_indexes, par_fixed=list()):
    if par_name in par_indexes:
        return par[par_indexes[par_name]]
    else:
        return par_fixed[par_name]


