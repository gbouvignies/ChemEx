"""Read the experimental CEST data."""
import importlib

import numpy as np


def read_profiles(path, filenames, details, model, included=None, excluded=None):
    """Read the CEST profiles."""

    experiment_type = details["type"].split(".")

    details["name"] = name_experiment(details)

    experiment_module = importlib.import_module(
        ".".join(
            ["chemex.experiments", experiment_type[0], "profiles", experiment_type[1]]
        )
    )

    if included is None:
        included = filenames.keys()

    if excluded is None:
        excluded = []

    included = [_.lower() for _ in included]
    excluded = [_.lower() for _ in excluded]

    Profile = getattr(experiment_module, "Profile")

    dtype = [("offsets", "f8"), ("intensities", "f8"), ("errors", "f8")]

    profiles = []

    for name, filename in filenames.items():

        full_path = path / filename

        data = np.loadtxt(full_path, dtype=dtype)

        profile = Profile(name, data, details, model)

        is_included = name.lower() in included
        is_not_excluded = name.lower() not in excluded

        if is_included and is_not_excluded:
            profiles.append(profile)

    error = details.get("error", "file")

    if error not in {"file", "auto"}:
        print("Warning: The 'error' option should either be 'file' or ")
        print("'auto'. Using the default 'file' option.")
        error = "file"

    if error == "auto":

        noise_values = []

        for profile in profiles:
            noise_values.append(profile.estimate_noise())

        noise_mean = np.mean(noise_values)

        for profile in profiles:
            profile.data["errors"] = noise_mean

    return profiles


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
