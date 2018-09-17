"""Read the experimental CEST data."""

import importlib
import os

import numpy as np
from scipy import interpolate, linalg, signal, stats


def read_profiles(path, filenames, details, model, included=None, excluded=None):
    """Read the CEST profiles."""

    experiment_type = details["type"].split(".")
    details["name"] = name_experiment(details)
    experiment_module = importlib.import_module(
        ".".join(
            ["chemex.experiments", experiment_type[0], "profiles", experiment_type[1]]
        )
    )

    error = details.get("error", "file")
    if error not in {"file", "auto"}:
        print("Warning: The 'error' option should either be 'file' or ")
        print("'auto'. Using the default 'file' option.")
        error = "file"

    error_values = []

    if included is None:
        included = filenames.keys()
    else:
        included = [_.lower() for _ in included]

    if excluded is None:
        excluded = []
    else:
        excluded = [_.lower() for _ in excluded]

    Profile = getattr(experiment_module, "Profile")

    dtype = [("b1_offsets", "<f8"), ("intensities", "<f8"), ("intensities_err", "<f8")]

    profiles = []

    for profile_name, filename in filenames.items():

        full_path = os.path.join(path, filename)

        measurements = np.loadtxt(full_path, dtype=dtype)
        measurements.sort(order="b1_offsets")

        profile = Profile(profile_name, measurements, details, model)

        is_included = profile_name.lower() in included
        is_not_excluded = profile_name.lower() not in excluded

        if is_included and is_not_excluded:
            profiles.append(profile)

        if error == "auto":
            excl_reference = np.logical_not(profile.reference)
            error_values.append(estimate_noise(profile.val[excl_reference]))

    error = details.get("error", "file")

    if error == "auto":

        error_value = np.mean(error_values)

        for profile in profiles:
            profile.err[:] = error_value

    ndata = sum(len(profile.val) for profile in profiles)

    return profiles, ndata


def name_experiment(experiment_details=None):
    """Generate a unique name for the experiment."""
    if experiment_details is None:
        experiment_details = dict()

    if "name" in experiment_details:
        name = experiment_details["name"].strip().replace(" ", "_")

    else:
        exp_type = experiment_details["type"].replace(".", "_")
        h_larmor_frq = float(experiment_details["h_larmor_frq"])
        temperature = float(experiment_details["temperature"])
        b1_frq = float(experiment_details["b1_frq"])
        time_t1 = float(experiment_details["time_t1"])

        name = "{:s}_{:.0f}Hz_{:.0f}ms_{:.0f}MHz_{:.0f}C".format(
            exp_type, b1_frq, time_t1 * 1e3, h_larmor_frq, temperature
        ).lower()

    return name


def estimate_noise(x):
    """Estimate the uncertainty in the CEST profile.

    # TODO: add reference to this method...

    """
    n = len(x)

    fda = [
        [1, -1],
        [1, -2, 1],
        [1, -3, 3, -1],
        [1, -4, 6, -4, 1],
        [1, -5, 10, -10, 5, -1],
        [1, -6, 15, -20, 15, -6, 1],
    ]

    fda = [np.array(a_fda) / linalg.norm(a_fda) for a_fda in fda]

    perc = np.array([0.05] + list(np.arange(0.1, 0.40, 0.025)))
    z = stats.norm.ppf(1.0 - perc)

    sigma_est = []

    for fdai in fda:
        noisedata = signal.convolve(x, fdai, mode="valid")
        ntrim = len(noisedata)

        if ntrim >= 2:
            noisedata.sort()

            p = 0.5 + np.arange(1, ntrim + 1)
            p /= ntrim + 0.5

            q = []

            f = interpolate.interp1d(p, noisedata, "linear")

            for a_perc, a_z in zip(perc, z):
                try:
                    val = (f(1.0 - a_perc) - f(a_perc)) / (2.0 * a_z)
                    q.append(val)
                except ValueError:
                    pass

            sigma_est.append(np.median(q))

    noisevar = np.median(sigma_est) ** 2
    noisevar /= 1.0 + 15.0 * (n + 1.225) ** -1.245

    return np.sqrt(noisevar)
