"""Plot the CEST profiles."""
import contextlib

import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf


def plot_data(data, params, output_dir):
    """Write experimental and fitted data to a file and plot the CEST profiles.

    - *.exp: contains the experimental data
    - *.fit: contains the fitted data
    - *.pdf: contains the plot of experimental and fitted data

    """
    datasets = dict()

    for data_point in data:
        experiment_name = data_point.experiment_name
        datasets.setdefault(experiment_name, []).append(data_point)

    for experiment_name, dataset in datasets.items():
        basename = output_dir / experiment_name
        name_pdf = basename.with_suffix(".pdf")
        name_exp = basename.with_suffix(".exp")
        name_fit = basename.with_suffix(".fit")

        print((f"  * {name_pdf} [.fit, .exp]"))

        with contextlib.ExitStack() as stack:

            file_pdf = stack.enter_context(backend_pdf.PdfPages(name_pdf))
            file_exp = stack.enter_context(name_exp.open("w"))
            file_fit = stack.enter_context(name_fit.open("w"))

            for profile in dataset:

                fig = profile.get_plot_fig(params)
                file_pdf.savefig(fig)
                plt.close()

                ndata_str, nfit_str = profile.plot_data_string(params)
                file_exp.write(ndata_str)
                file_fit.write(nfit_str)

    return
