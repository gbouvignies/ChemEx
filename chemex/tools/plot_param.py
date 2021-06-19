import configparser
import sys

import chemex.nmr.spin_system as cns
import chemex.parameters.name as cpn


def plot_param(args):
    """Plot values of a parameter versus residue number."""
    import matplotlib.pyplot as plt

    params = configparser.ConfigParser()

    if len(args.parameters) > 1:
        sys.exit(
            "\nError: Multiple parameter files were given. 'chemex plot_param' "
            "should only be run with a single parameter file.\n"
        )

    params.read(str(args.parameters.pop()))
    parname = args.parname
    curves = {}
    print("Plotting...")
    for section in params.sections():
        section_ = section.strip('"')
        short_name = cpn.ParamName().from_section(section_).name
        if parname.lower() in short_name:
            print(f"  - [{section}]")
            points = []
            for key, value in params.items(section):
                res = int(cns.SpinSystem(key).numbers["i"])
                split = value.split()
                value = float(split[0])
                try:
                    error = float(split[2].strip("±"))
                except ValueError:
                    error = 0.0
                points.append((res, value, error))
            curves[section] = zip(*points)
    _, axis = plt.subplots(figsize=(12, 5))
    axis.yaxis.grid(True)
    for section, (res, vals, errors) in curves.items():
        axis.errorbar(res, vals, yerr=errors, label=section, fmt=".", barsabove=True)
    plt.legend()
    plt.show()
    return 0
