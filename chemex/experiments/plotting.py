import matplotlib as mpl


dark_gray = '0.13'

base_context = {
    'figure.figsize': (8, 5.5),
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,

    'grid.linewidth': 1.0,
    'patch.linewidth': 0.3,
    'lines.linewidth': 1.5,

    'lines.markersize': 6.0,
    'lines.markeredgewidth': 0.0,

    'xtick.major.width': 1,
    'xtick.minor.width': .5,
    'xtick.major.pad': 7,

    'ytick.major.width': 1,
    'ytick.minor.width': .5,
    'ytick.major.pad': 7,
}

# Common parameters
style_dict = {
    'text.color': dark_gray,

    'xtick.color': dark_gray,
    'xtick.major.size': 0,
    'xtick.minor.size': 0,

    'ytick.color': dark_gray,
    'ytick.major.size': 0,
    'ytick.minor.size': 0,

    'grid.color': 'black',
    'grid.alpha': 0.12,
    'grid.linestyle': '-',

    'axes.labelcolor': dark_gray,
    'axes.facecolor': 'white',  # axes background color
    'axes.edgecolor': 'white',  # axes edge color
    'axes.grid': True,
    'axes.axisbelow': True,
    'axes.hold': True,


    'lines.solid_capstyle': 'round',

    'font.family': ['sans-serif'],
    'font.sans-serif': ['Arial', 'Liberation Sans',
                        'Bitstream Vera Sans', 'sans-serif'],
}

mpl.rcParams.update(base_context)
mpl.rcParams.update(style_dict)
mpl.use('Agg')


def plot_data(data, params, output_dir='./'):
    """ Plot all data types """

    subsets = dict()

    for profile in data:
        subsets.setdefault(profile.plot_data, []).append(profile)

    for plot, dataset in subsets.items():
        plot(dataset, params, output_dir)

    return



