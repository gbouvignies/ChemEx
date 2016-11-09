"""TODO: module docstring."""

import abc


class BaseProfile(metaclass=abc.ABCMeta):
    """TODO: class docstring."""

    @abc.abstractmethod
    def print_profile(self):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def calculate_profile(self, params):
        """TODO: method docstring."""
        pass

    def create_default_parameters(self):
        """TODO: method docstring."""
        pass

    @abc.abstractmethod
    def filter_points(self, params=None):
        """TODO: method docstring."""
        pass

    def calculate_residuals(self, params):
        """Calculate the residuals between the experimental and
        back-calculated values.
        """
        values = self.calculate_profile(params)

        return (self.val - values) / self.err


def check_par(parameters=None, name=None, convert=None, default=None, required=True):
    """Check for experimental parameters and converts them to their
    appropriate type.
    """
    value = parameters.get(name, default)

    if required and value is None:
        exit("Missing experimental parameter detected. Please set: {:s}".format(name))

    if convert is not None and value is not None:
        try:
            value = convert(value)
        except ValueError:
            exit(
                "Experimental parameter of wrong type detected. Please make"
                " sure that {:s} is a {:s}".format(name, convert)
            )

    return value
