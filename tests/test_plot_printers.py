import numpy as np

from chemex.containers.data import Data
from chemex.printers.plot import data_plot_printers


def test_cest_fit_printer_keeps_ppm_offsets_distinct() -> None:
    data = Data(
        exp=np.zeros(2),
        err=np.zeros(2),
        metadata=np.array([14.97012, 14.97067]),
    )
    data.calc = np.array([0.1, 0.2])

    output = data_plot_printers["cest"].print_calc("G23CD2-HD2", data)

    lines = output.splitlines()
    assert lines[2].split()[0] == "14.97012"
    assert lines[3].split()[0] == "14.97067"
