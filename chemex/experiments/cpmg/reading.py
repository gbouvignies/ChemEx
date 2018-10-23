"""Read the experimental CPMG data."""
import numpy as np

from chemex import experiments


def read_profiles(path, filenames, details, model):
    """Read the CPMG profiles."""

    details["name"] = name_experiment(details)
    Profile = experiments.grab(details["type"])
    dtype = [("ncycs", "i4"), ("intensity", "f8"), ("error", "f8")]

    profiles = []

    for profile_name, filename in filenames.items():
        full_path = path / filename
        data = np.loadtxt(full_path, dtype=dtype)
        profiles.append(Profile(profile_name, data, details, model))

    error = details.get("error", "file")

    if error not in {"file", "duplicates"}:
        print("Warning: The 'error' option should either be 'file' or ")
        print("'duplicates'. Using the default 'file' option.")
        error = "file"

    if error == "duplicates":
        # Estimate the uncertainty using the pooled standard deviation
        # Ref: http://goldbook.iupac.org/html/P/P04758.html

        variances = []

        for profile in profiles:
            variances.append(profile.get_variance_from_duplicates())

        noise_mean = np.sqrt(np.mean(variances))

        for profile in profiles:
            profile.data["error"] = noise_mean

    return profiles


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
