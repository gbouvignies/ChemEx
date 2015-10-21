from __future__ import absolute_import
import abc

import six


class BaseProfile(six.with_metaclass(abc.ABCMeta, object)):
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
    def create_default_parameters(self):
        pass

    @abc.abstractmethod
    def filter_points(self, parameters=None):
        pass


def check_par(parameters=None, name=None, convert=None, default=None):
    """Sets experimental parameter and converts it to the right type."""

    value = None

    if name not in parameters and default is not None:
        parameters[name] = default

    try:
        value = parameters[name]
    except KeyError:
        exit(
            "Missing experimental parameter detected. Please set: {:s}"
                .format(name)
        )

    if convert is not None:
        try:
            value = convert(value)
        except ValueError:
            exit(
                "Experimental parameter of wrong type detected. Please make"
                " sure that {:s} is a {:s}".format(name, convert)
            )

    return value
