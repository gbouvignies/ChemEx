import abc


class BaseProfile(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def print_profile(self):
        pass

    @abc.abstractmethod
    def calculate_profile(self, parameters, *args, **kwargs):
        pass

    @abc.abstractmethod
    def calculate_residuals(self, parameters):
        pass

    @abc.abstractmethod
    def make_default_parameters(self, parameters=None):
        pass

    @abc.abstractmethod
    def filter_points(self, parameters=None):
        pass


def check_par(parameters=None, name=None, func=None):
    """Sets experimental parameter and converts it to the right type."""

    value = None

    try:
        value = parameters[name]
    except KeyError:
        exit(
            "Missing experimental parameter detected. Please set: {:s}"
                .format(name)
        )

    if func is not None:
        try:
            value = func(value)
        except ValueError:
            exit(
                "Experimental parameter of wrong type detected. Please make"
                " sure that {:s} is a {:s}".format(name, func)
            )

    return value

