import abc


class BaseProfile(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def print_profile(self):
        pass

    @abc.abstractmethod
    def calculate_profile(self, params):
        pass

    def create_default_parameters(self):
        pass

    @abc.abstractmethod
    def filter_points(self, params=None):
        pass

    def calculate_residuals(self, params):
        """Calculates the residual between the experimental and
        back-calculated values.
        """

        values = self.calculate_profile(params)

        return (self.val - values) / self.err


def check_par(parameters=None, name=None, convert=None, default=None, required=True):
    """Sets experimental parameter and converts it to the right type."""

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
