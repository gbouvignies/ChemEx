"""Read the experimental CPMG data."""

import importlib
import math
import os

import numpy as np


def read_profiles(path, filenames, details, model, included=None, excluded=None):
    """Read the CPMG profiles."""
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

    dtype = [("ncycs", "<f8"), ("intensities", "<f8"), ("intensities_err", "<f8")]

    profiles = []

    for profile_name, filename in filenames.items():

        full_path = os.path.join(path, filename)

        measurements = np.loadtxt(full_path, dtype=dtype)

        profile = Profile(profile_name, measurements, details, model)

        is_included = profile_name in included
        is_not_excluded = profile_name not in excluded

        if is_included and is_not_excluded:
            profiles.append(profile)

        if error == "auto":
            error_values.append(estimate_noise(profile))

    if error == "auto":
        error_value = np.mean(error_values)
        for profile in profiles:
            profile.err[:] = error_value

    ndata = sum(len(profile.val) for profile in profiles)

    return profiles, ndata


def name_experiment(experiment_details=None):
    """Generate a unique name for the experiment."""
    if not experiment_details:
        experiment_details = dict()

    if "name" in experiment_details:
        name = experiment_details["name"].strip().replace(" ", "_")

    else:
        exp_type = experiment_details["type"].replace(".", "_")
        h_larmor_frq = float(experiment_details["h_larmor_frq"])
        temperature = float(experiment_details["temperature"])
        time_t2 = float(experiment_details["time_t2"])

        name = "{:s}_{:.0f}ms_{:.0f}MHz_{:.0f}C".format(
            exp_type, time_t2 * 1e3, h_larmor_frq, temperature
        ).lower()

    return name


def correction(n):
    """Calculate correction factor for noise estimate."""
    k = n // 2

    if n == 2 * k:
        factor = (
            math.sqrt(math.pi * (2 * k - 1) / 2.0)
            * math.factorial(2 * k - 2)
            / ((2 ** (2 * k - 2)) * (math.factorial(k - 1) ** 2))
        )
    else:
        factor = (
            math.sqrt(k / math.pi)
            * (2 ** (2 * k - 1))
            * (math.factorial(k - 1) ** 2)
            / math.factorial(2 * k - 1)
        )

    return factor


def estimate_noise(profile):
    """Estimate the uncertainty in CPMG data points (i.e., R2eff)."""
    intensity_dict = {}

    for ncyc, intensity in zip(profile.ncycs, profile.val):
        intensity_dict.setdefault(ncyc, []).append(intensity)

    std_list = []
    for duplicates in intensity_dict.values():
        n_duplicates = len(duplicates)
        if n_duplicates > 1:
            std_list.append(np.std(duplicates, ddof=1) * correction(n_duplicates))

    if std_list:
        error = np.mean(std_list)
    else:
        error = np.mean(profile.err)

    return error
