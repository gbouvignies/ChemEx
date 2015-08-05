import matplotlib as mpl

dark_gray = '0.13'
medium_gray = '0.44'

base_context = {
    'figure.figsize': (8, 5.5),
    'figure.subplot.bottom': 0.13,
    'figure.subplot.top': 0.93,
    'figure.subplot.left': 0.10,
    'figure.subplot.right': 0.95,
    'figure.subplot.hspace': 0.10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,

    'grid.linewidth': 0.5,
    'patch.linewidth': 0.3,
    'lines.linewidth': 1.5,
    'axes.linewidth': 0.5,

    'lines.markersize': 6.0,
    'lines.markeredgewidth': 0.0,

    'xtick.major.width': .5,
    'xtick.minor.width': .5,
    'xtick.major.pad': 7,

    'ytick.major.width': .5,
    'ytick.minor.width': .5,
    'ytick.major.pad': 7,
}

style_dict = {
    'text.color': dark_gray,

    'xtick.color': dark_gray,
    'xtick.major.size': 3,
    'xtick.minor.size': 3,
    'xtick.direction': 'outside',

    'ytick.color': dark_gray,
    'ytick.major.size': 3,
    'ytick.minor.size': 3,
    'ytick.direction': 'outside',

    'grid.color': '0.88',
    'grid.linestyle': '-',

    'axes.labelcolor': dark_gray,
    'axes.facecolor': 'white',  # axes background color
    'axes.edgecolor': dark_gray,  # axes edge color
    'axes.grid': True,
    'axes.axisbelow': True,
    'axes.hold': True,

    'lines.solid_capstyle': 'round',

    'font.family': ['sans-serif'],
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif'],
}

mpl.rcParams.update(base_context)
mpl.rcParams.update(style_dict)


def plot_data(data, params, output_dir='./'):
    """ Plot all data types """

    subsets = dict()

    for profile in data:
        subsets.setdefault(profile.plot_data, []).append(profile)

    for plot, dataset in subsets.items():
        plot(dataset, params, output_dir)

    return
