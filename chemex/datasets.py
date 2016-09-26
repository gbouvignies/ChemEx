# -*- coding: utf-8 -*-


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
from scipy import stats


class DataSet(object):
    def __init__(self, other=None):

        self.data = []
        self.ndata = 0
        self.chisq_ref = 1e32

        if isinstance(other, DataSet):
            self.data = copy.deepcopy(other.data)
            self.ndata = other.ndata

        elif isinstance(other, base_profile.BaseProfile):
            self.data.append(other)
            self.ndata = len(other.val)

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
        if not isinstance(profile, base_profile.BaseProfile):
            raise TypeError
        self.data.append(profile)
        self.ndata += len(profile.val)

    def calculate_residuals(self, params, verbose=True, threshold=1e-3):

        residuals = np.concatenate([profile.calculate_residuals(params) for profile in self.data])

        if verbose:

            chisq = sum(residuals ** 2)

            change = (chisq - self.chisq_ref) / self.chisq_ref

            if change < -threshold:
                nvarys = len([param for param in params.values() if param.vary])
                redchi = chisq / (self.ndata - nvarys)

                print("  * {:.3e} / {:.3e}".format(chisq, redchi))

                self.chisq_ref = chisq

        return residuals

    def calculate_chisq(self, params):

        residuals = self.calculate_residuals(params, verbose=False)

        return sum(residuals ** 2)

    def calculate_redchi(self, params):

        chisq = self.calculate_chisq(params)
        nvarys = len([param for param in params.values() if param.vary])

        return chisq / (self.ndata - nvarys)

    def write_to(self, params=None, output_dir='./'):
        """Write dispersion profiles into a file"""

        datasets = dict()

        for profile in self.data:
            experiment_name = profile.experiment_name
            datasets.setdefault(experiment_name, list()).append(profile)

        for experiment_name, data in datasets.items():

            filename = ''.join([experiment_name, '.dat'])
            filename = os.path.join(output_dir, filename)

            print("  * {}".format(filename))

            with open(filename, 'w') as f:
                for profile in sorted(data, key=operator.attrgetter('peak')):
                    f.write(profile.print_profile(params=params))

    def write_chi2_to(self, params, path='./'):
        """Write reduced chi2
        """

        residuals = self.calculate_residuals(params, verbose=False)
        chisq = sum(residuals ** 2)
        nvarys = len([param for param in params.values() if param.vary])
        nfree = self.ndata - nvarys
        redchi = chisq / nfree
        ks_value, ks_p_value = stats.kstest(residuals, 'norm')
        chi2_p_value = 1.0 - stats.chi2.cdf(chisq, nfree)

        filename = os.path.join(path, 'chi2.fit')

        with open(filename, 'w') as f:
            print("  * {}".format(filename))

            f.write("# Number of values: {}\n".format(self.ndata))
            f.write("# Number of varying parameters: {}\n".format(nvarys))
            f.write("chisq         = {: .5e}\n".format(chisq))
            f.write("redchi        = {: .5e}\n".format(redchi))
            f.write("chisq_p-value = {: .5e} # Chi-squared test\n".format(chi2_p_value))
            f.write("ks_p-value    = {: .5e} # Kolmogorov-Smirnov test for goodness of fit\n".format(ks_p_value))

    def add_dataset_from_file(self, filename, model=None, res_incl=None, res_excl=None):

        print("{:<45s} ".format(filename), end='')

        # Get the directory of the input file
        working_dir = os.path.dirname(filename)

        # Parse the config file
        config = util.read_cfg_file(filename)

        try:

            # Reads experiment information
            experiment_details = dict(config.items('experiment'))
            experiment_type = experiment_details['type']
            experiment_class = experiment_type.split('.')[0]

            # Reads experimental parameters
            experiment_details.update(
                {key.lower(): val for key, val in config.items('experimental_parameters')}
            )

            # Reads profile information (name, filename)
            profile_filenames = {key.lower(): val for key, val in config.items('data')}

            experiment_details['model'] = model

        except configparser.NoSectionError as e:
            sys.exit("    Reading aborted: {}".format(e))

        except KeyError as e:
            sys.exit("\nIn the section 'experiment' of {}, '{}' must be provided!".format(filename, e))

        try:
            # Reads additional parameters
            experiment_details.update(
                {key.lower(): val for key, val in config.items('extra_parameters')}
            )

        except configparser.NoSectionError:
            pass

        path = util.normalize_path(working_dir, experiment_details.get('path', './'))

        try:
            reading = importlib.import_module('.'.join(['chemex.experiments', experiment_class, 'reading']))

        except ImportError:
            sys.exit("The experiment '{}', referred in '{}' is not implemented.".format(experiment_type, filename))

        data, ndata = reading.read_profiles(path, profile_filenames, experiment_details, res_incl, res_excl)

        self.data.extend(data)
        self.ndata += ndata

        print("{:<25s} {:<25d}".format(experiment_type, len(data)))

        return data
