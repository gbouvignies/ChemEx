"""
Created on Oct 11, 2013

@author: guillaume
"""


def plot_data(data, par, par_names, par_fixed, output_dir='./'):
    """ Plot all data types """

    subsets = dict()

    for data_point in data:
        subsets.setdefault(data_point.plot_data, []).append(data_point)

    for plot, dataset in subsets.iteritems():
        plot(dataset, par, par_names, par_fixed, output_dir)

    return


