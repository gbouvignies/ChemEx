"""Read the experimental CEST data."""
import numpy as np

from chemex import experiments


def read_profiles(path, filenames, details, model):
    """Read the CEST profiles."""

    details["name"] = name_experiment(details)
    Profile = experiments.grab(details["type"])
    dtype = [("offsets", "f8"), ("intensity", "f8"), ("error", "f8")]

    profiles = []

    for name, filename in filenames.items():
        full_path = path / filename
        data = np.loadtxt(full_path, dtype=dtype)
        profiles.append(Profile(name, data, details, model))

    error = details.get("error", "file")

    if error not in {"file", "scatter"}:
        print("Warning: The 'error' option should either be 'file' or ")
        print("'scatter'. Using the default 'file' option.")
        error = "file"

    if error == "scatter":

        noise_values = []

        for profile in profiles:
            noise_values.append(profile.estimate_noise())

        noise_mean = np.mean(noise_values)

        for profile in profiles:
            profile.data["error"] = noise_mean

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
