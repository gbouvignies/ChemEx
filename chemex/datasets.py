"""The dataset module contains the code for handling the experimental data."""

import configparser
import copy
import importlib
import operator
import os
import os.path
import sys

import numpy as np

from chemex import util
from chemex.experiments import base_profile


class DataSet(object):
    """DataSet class for handling experimental data."""

    def __init__(self, other=None):
        self.data = []
        self.ndata = 0
        self.chisq_ref = 1e32

        if isinstance(other, DataSet):
            self.data = copy.deepcopy(other.data)
            self.ndata = other.ndata

        elif isinstance(other, base_profile.BaseProfile):
            self.data.append(other)
            self.ndata = len(other)

    def __iter__(self):
        for some_data in self.data:
            yield some_data

    def __add__(self, other):
        data_sum = DataSet(self)
        data_sum.data.extend(other.data)
        data_sum.ndata = self.ndata + other.ndata
        return data_sum

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        self.data.extend(other.data)
        self.ndata += other.ndata
        return self

    def append(self, profile):
        """Append data to an exisiting dataset."""
        if not isinstance(profile, base_profile.BaseProfile):
            raise TypeError
        self.data.append(profile)
        self.ndata += len(profile.val)

    def calculate_residuals(self, params, verbose=True, threshold=1e-3):
        """Calculate the residuals."""

        residuals = np.concatenate(
            [profile.calculate_residuals(params) for profile in self.data]
        )

        if verbose:
            chisq = sum(residuals ** 2)

            change = (chisq - self.chisq_ref) / self.chisq_ref

            if change < -threshold:
                nvarys = len([param for param in params.values() if param.vary])
                redchi = chisq / (self.ndata - nvarys)

                print("  * {:.3e} / {:.3e}".format(chisq, redchi))

                self.chisq_ref = chisq

        return residuals

    def write_to(self, params=None, output_dir="./"):
        """Write experimental and fitted profiles to a file."""
        datasets = dict()

        for profile in self.data:
            experiment_name = profile.experiment_name
            datasets.setdefault(experiment_name, list()).append(profile)

        for experiment_name, data in datasets.items():
            filename = "".join([experiment_name, ".dat"])
            filename = os.path.join(output_dir, filename)

            print("  * {}".format(filename))

            with open(filename, "w") as f:
                for profile in sorted(data, key=operator.attrgetter("peak")):
                    f.write(profile.print_profile(params=params))

    def add_dataset_from_file(self, filename, model=None, included=None, excluded=None):
        """Add data from a file to the dataset."""
        print("{:<45s} ".format(filename), end="")

        # Get the directory of the input file
        working_dir = os.path.dirname(filename)

        # Parse the experiment configuration file
        config = util.read_cfg_file(filename)

        try:
            # Read the experiment information
            details = dict(config.items("experiment"))
            experiment_type = details["type"]
            experiment_class = experiment_type.split(".")[0]

            # Read the experimental parameters
            details.update(
                {
                    key.lower(): val
                    for key, val in config.items("experimental_parameters")
                }
            )

            # Read the profile information (name, filename)
            filenames = {key.lower(): val for key, val in config.items("data")}

        except configparser.NoSectionError as e:
            sys.exit("    Reading aborted: {}".format(e))

        except KeyError as e:
            sys.exit(
                "\nIn the section 'experiment' of {}, '{}' must be provided!".format(
                    filename, e
                )
            )

        try:
            # Read (optional) additional parameters
            details.update(
                {key.lower(): val for key, val in config.items("extra_parameters")}
            )

        except configparser.NoSectionError:
            pass

        path = util.normalize_path(working_dir, details.get("path", "./"))

        try:
            reading = importlib.import_module(
                ".".join(["chemex.experiments", experiment_class, "reading"])
            )

        except ImportError:
            sys.exit(
                "The experiment '{}', referred in '{}' is not implemented.".format(
                    experiment_type, filename
                )
            )

        data, ndata = reading.read_profiles(
            path, filenames, details, model, included, excluded
        )

        self.data.extend(data)
        self.ndata += ndata

        print("{:<25s} {:<25d}".format(experiment_type, len(data)))

        return data

    def make_bs_dataset(self):
        """Create a new dataset to run a bootstrap simulation."""

        data_bs = DataSet()

        for profile in self.data:
            data_bs.append(profile.make_bs_profile())

        return data_bs

    def make_mc_dataset(self, params):
        """Create a new dataset to run a Monte-Carlo simulation."""

        data_mc = DataSet()

        for profile in self.data:
            data_mc.append(profile.make_mc_profile(params=params))

        return data_mc


def read_data(args):
    """Read experimental setup and data."""
    util.header1("Reading Experimental Data")

    data = DataSet()

    if args.experiments:
        print(("{:<45s} {:<25s} {:<25s}".format("File Name", "Experiment", "Profiles")))
        print(("{:<45s} {:<25s} {:<25s}".format("---------", "----------", "--------")))

        for filename in args.experiments:
            data.add_dataset_from_file(
                filename, args.model, args.res_incl, args.res_excl
            )

    if not data.data:
        sys.exit("\nNo data to fit!\n")

    return data
