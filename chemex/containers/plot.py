import numpy as np
import scipy.interpolate as ip


_GREY400 = "#BDBDBD"
_LSTYLES = ("-", "--", "-.", ":")
_RED100 = "#FFCDD2"
_RED300 = "#E57373"
_RED700 = "#D32F2F"


def profile(name, data_exp, data_fit):
    import matplotlib.figure as mf

    xname, yname, ename, *_ = data_exp.dtype.names
    residuals = _get_residuals(data_exp, data_fit)
    range_x = get_grid(data_fit[xname], 2, 0.02)
    fig = mf.Figure()
    ax1, ax2 = fig.subplots(2, 1, sharex="all", gridspec_kw={"height_ratios": [1, 4]})
    fig.align_labels()
    fig.suptitle(f"{str(name).upper()}")
    ax1.set_xlim(*range_x)
    ax1.set_ylabel("Residuals")
    ax1.ticklabel_format(style="sci", scilimits=(0, 0), axis="y", useMathText=True)
    ax1.errorbar(
        data_exp[data_exp.mask][xname],
        residuals[data_exp.mask],
        abs(data_exp[data_exp.mask][ename].T),
        fmt=".",
        color=_RED700,
        zorder=3,
    )
    ax1.errorbar(
        data_exp[~data_exp.mask][xname],
        residuals[~data_exp.mask],
        abs(data_exp[~data_exp.mask][ename].T),
        fmt=".",
        color=_RED100,
        zorder=3,
    )
    ax2.plot(data_fit[xname], data_fit[yname], linestyle="-", color=_RED300)
    ax2.errorbar(
        data_exp[data_exp.mask][xname],
        data_exp[data_exp.mask][yname],
        abs(data_exp[data_exp.mask][ename].T),
        fmt=".",
        color=_RED700,
        zorder=3,
    )
    ax2.errorbar(
        data_exp[~data_exp.mask][xname],
        data_exp[~data_exp.mask][yname],
        abs(data_exp[~data_exp.mask][ename].T),
        fmt=".",
        color=_RED100,
        zorder=3,
    )
    range_y = ax2.get_ylim()
    for axis in (ax1, ax2):
        axis.axhline(0, color="k", linewidth=0.5, zorder=1)
    range_y = ax2.set_ylim(range_y)
    return fig


def cpmg(file_pdf, name, data_exp, data_fit):
    xname, *_ = data_fit.dtype.names
    fig = profile(name, data_exp, data_fit)
    ax2 = fig.axes[1]
    ax2.set_xlabel(r"$\nu_\mathregular{CPMG}$ (Hz)")
    ax2.set_ylabel(r"$R_{2,\mathregular{eff}}$ (s$^{-1}$)")
    ax2.set_xlim(0.0, max(data_fit[xname]) + min(data_fit[xname]))
    file_pdf.savefig(fig)


def relaxation(file_pdf, name, data_exp, data_fit):
    xname, *_ = data_fit.dtype.names
    fig = profile(name, data_exp, data_fit)
    ax2 = fig.axes[1]
    ax2.set_xlabel(r"Time (s)")
    ax2.set_ylabel(r"Intensity")
    ax2.set_xlim(0.0, max(data_fit[xname]) + min(data_fit[xname]))
    file_pdf.savefig(fig)


def cest(file_pdf, name, data_exp, data_fit, cs_values, alias_values):
    import matplotlib.ticker as ticker

    residuals = _get_residuals(data_exp, data_fit)
    sigma = _estimate_sigma(residuals)
    fig = profile(name, data_exp, data_fit)
    ax1, ax2 = fig.axes
    lim = sorted(ax2.get_xlim(), reverse=True)
    ax1.set_xlim(lim)
    ax2.set_xlim(lim)
    ax1.xaxis.set_major_locator(ticker.MaxNLocator(6))
    ax2.xaxis.set_major_locator(ticker.MaxNLocator(6))
    ax2.set_xlabel(r"$B_1$ position (ppm)")
    ax2.set_ylabel(r"$I/I_0$")
    kwargs1 = {"facecolor": (0, 0, 0, 0.1), "edgecolor": "none"}
    ax1.fill_between(ax1.get_xlim(), -1.0 * sigma, 1.0 * sigma, **kwargs1)
    ax1.fill_between(ax1.get_xlim(), -2.0 * sigma, 2.0 * sigma, **kwargs1)
    kwargs2 = {"color": _GREY400, "linewidth": 0.75, "zorder": -1}
    for a_cs, alias, lstyle in zip(cs_values, alias_values, _LSTYLES):
        ax1.axvline(a_cs, linestyle=lstyle, **kwargs2)
        ax2.axvline(a_cs, linestyle=lstyle, **kwargs2)
        if alias:
            x, _ = ax2.transLimits.transform((a_cs, 0))
            ax2.text(x - 0.02, 0.95, "*", transform=ax2.transAxes)
    file_pdf.savefig(fig)


def shift(name_pdf, name, fit, exp, err):
    import matplotlib.figure as mf

    fig = mf.Figure()
    ax = fig.subplots(1, 1)
    fig.align_labels()
    ax.errorbar(fit, exp, yerr=err, fmt=".", color=_RED700)
    val_min = min(ax.get_xlim()[0], ax.get_ylim()[0])
    val_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_aspect("equal", "box")
    ax.plot([val_min, val_max], [val_min, val_max], color="k", linewidth=0.5, zorder=1)
    ax.set_xlabel(r"$Δδ_\mathregular{fit}$ (ppb)")
    ax.set_ylabel(r"$Δδ_\mathregular{exp}$ (ppb)")
    fig.savefig(name_pdf)


def _get_residuals(data_exp, data_fit):
    xname, yname, *_ = data_exp.dtype.names
    data_fit_ = np.unique(np.sort(data_fit, order=xname))
    data_fit_f = ip.interp1d(data_fit_[xname], data_fit_[yname], "cubic")
    return data_exp[yname] - data_fit_f(data_exp[xname])


def _estimate_sigma(values):
    """Estimates standard deviation using median to exclude outliers.

    Up to 50% can be bad.

    Ref :
    Rousseeuw, Peter & Croux, Christophe. (1993). Alternatives to Median Absolute
    Deviation. Journal of the American Statistical Association. 88. 1273 - 1283.
    10.1080/01621459.1993.10476408.

    """
    if not all(values):
        return 0.0
    _values = values.reshape(1, -1)
    return 1.1926 * np.median(np.median(abs(_values - _values.T), axis=0))


def get_grid(values, size=500, extension=0.0):
    value_min = np.min(values)
    value_max = np.max(values)
    extra = (value_max - value_min) * extension
    return np.linspace(value_min - extra, value_max + extra, size)


def write_plots(experiments, params, path, simulation=False):
    """Plot the experimental and fitted data."""
    print("\nPlotting data...")
    path_ = path / "Plots"
    path_.mkdir(parents=True, exist_ok=True)
    try:
        experiments.plot(path=path_, params=params, simulation=simulation)
    except KeyboardInterrupt:
        print("  - Plotting cancelled")
