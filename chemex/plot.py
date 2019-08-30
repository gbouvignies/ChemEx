import matplotlib.figure as mf
import numpy as np
import scipy.interpolate as ip


_GREY400 = "#BDBDBD"
_LSTYLES = ("-", "--", "-.", ":")
_RED100 = "#FFCDD2"
_RED300 = "#E57373"
_RED700 = "#D32F2F"


def profile(name, data_exp, data_fit):
    xname, yname, ename, *_ = data_exp.dtype.names
    residuals = _get_residuals(data_exp, data_fit)
    xrange = get_grid(data_fit[xname], 2, 0.02)
    data_y = np.concatenate(
        [(data_exp[yname] + data_exp[ename].T).flatten(), data_fit[yname]]
    )
    yrange2 = get_grid(data_y, 2, 0.1)
    yrange1 = get_grid((residuals + data_exp[ename].T).flatten(), 2, 0.1)
    fig = mf.Figure()
    ax1, ax2 = fig.subplots(2, 1, sharex="all", gridspec_kw={"height_ratios": [1, 4]})
    fig.align_labels()
    fig.suptitle(f"{str(name).upper()}")
    ax1.set_xlim(*xrange)
    ax1.set_ylim(*yrange1)
    ax2.set_ylim(*yrange2)
    ax1.set_ylabel("Residuals")
    ax1.ticklabel_format(style="sci", scilimits=(0, 0), axis="y", useMathText=True)
    for axis in (ax1, ax2):
        axis.axhline(0, color="k", linewidth=0.5)
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
    return fig


def cpmg(file_pdf, name, data_exp, data_fit):
    xname, *_ = data_exp.dtype.names
    fig = profile(name, data_exp, data_fit)
    ax1, ax2 = fig.axes
    ax2.set_xlabel(r"$\nu_\mathregular{CPMG}$ (Hz)")
    ax2.set_ylabel(r"$R_{2,\mathregular{eff}}$ (s$^{-1}$)")
    ax2.set_xlim(0.0, max(data_exp[xname]) + min(data_exp[xname]))
    file_pdf.savefig(fig)


def cest(file_pdf, name, data_exp, data_fit, cs_values):
    residuals = _get_residuals(data_exp, data_fit)
    sigma = _estimate_sigma(residuals)
    fig = profile(name, data_exp, data_fit)
    ax1, ax2 = fig.axes
    ax2.set_xlabel(r"$B_1$ position (ppm)")
    ax2.set_ylabel(r"$I/I_0$")
    kwargs1 = {"facecolor": (0, 0, 0, 0.1), "edgecolor": "none"}
    ax1.fill_between(ax1.get_xlim(), -1.0 * sigma, 1.0 * sigma, **kwargs1)
    ax1.fill_between(ax1.get_xlim(), -2.0 * sigma, 2.0 * sigma, **kwargs1)
    kwargs2 = {"color": _GREY400, "linewidth": 0.75}
    for a_cs, lstyle in zip(cs_values, _LSTYLES):
        ax1.axvline(a_cs, linestyle=lstyle, zorder=-1, **kwargs2)
        ax2.axvline(a_cs, linestyle=lstyle, zorder=-1, **kwargs2)
    lim = sorted(ax2.get_xlim(), reverse=True)
    ax1.set_xlim(lim)
    ax2.set_xlim(lim)
    file_pdf.savefig(fig)


def _get_residuals(data_exp, data_fit):
    xname, yname, ename, *_ = data_exp.dtype.names
    data_fit_ = np.unique(np.sort(data_fit, order=xname))
    data_fit_f = ip.interp1d(data_fit_[xname], data_fit_[yname], "cubic")
    residuals = data_exp[yname] - data_fit_f(data_exp[xname])
    return residuals


def _estimate_sigma(values):
    """Estimates standard deviation using median to exclude outliers.

    Up to 50% can be bad.

    Ref :
    Rousseeuw, Peter & Croux, Christophe. (1993). Alternatives to Median Absolute
    Deviation. Journal of the American Statistical Association. 88. 1273 - 1283.
    10.1080/01621459.1993.10476408.

    """
    _values = values.reshape(1, -1)
    return 1.1926 * np.median(np.median(abs(_values - _values.T), axis=0))


def get_grid(values, size=500, extension=0.0):
    value_min = np.min(values)
    value_max = np.max(values)
    extra = (value_max - value_min) * extension
    return np.linspace(value_min - extra, value_max + extra, size)
