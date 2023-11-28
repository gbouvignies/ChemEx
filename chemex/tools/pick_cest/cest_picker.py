import sys
from argparse import Namespace
from collections.abc import Callable
from typing import Any

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.backend_bases import Event
from matplotlib.figure import Figure
from matplotlib.widgets import Button, Slider

from chemex.configuration.methods import Selection
from chemex.containers.experiment import Experiment
from chemex.containers.experiments import Experiments
from chemex.experiments.builder import build_experiments
from chemex.experiments.loader import register_experiments
from chemex.models import model
from chemex.tools.pick_cest.buttons import Buttons

FigureSizeType = tuple[float, float]
AdjustParamsType = tuple[float, float, float, float]
PositionType = tuple[float, float, float, float]

FIGURE_SIZE = (12.0, 5.0)
ADJUST_PARAMS = (0.07, 0.1, 0.8, 0.9)
POSITION_SCALE = (0.855, 0.425, 0.1, 0.02)
POSITION_PREVIOUS = (0.825, 0.3, 0.075, 0.075)
POSITION_NEXT = (0.9, 0.3, 0.075, 0.075)
POSITION_SWAP = (0.825, 0.2, 0.15, 0.075)
POSITION_CLEAR = (0.825, 0.1, 0.15, 0.075)


def create_figure(
    size: FigureSizeType = FIGURE_SIZE,
    adjust_params: AdjustParamsType = ADJUST_PARAMS,
) -> tuple[Figure, Axes]:
    """Creates and returns a matplotlib figure and axes."""
    fig, axes = plt.subplots()
    fig.set_size_inches(*size)
    fig.subplots_adjust(*adjust_params)
    return fig, axes


def is_cest_experiment(experiment: Experiment) -> bool:
    """Checks if an experiment is a CEST experiment."""
    return experiment.name.startswith(("cest", "dcest", "coscest"))


def initialize_cest(args: Namespace) -> Experiments:
    """Initializes CEST experiments."""
    register_experiments()
    model.set_model("2st")
    return build_experiments(args.experiments, Selection(include=None, exclude=None))


def create_slider(
    fig: Figure,
    position: PositionType,
    label: str,
    val_range: tuple[float, float],
    initial_value: float,
    on_change: Callable[[float], None],
) -> Slider:
    """Creates and returns a matplotlib slider."""
    ax = fig.add_axes(position)
    slider = Slider(ax, label, *val_range, valfmt="%3.1f", valinit=initial_value)
    slider.on_changed(on_change)
    return slider


def create_button(
    fig: Figure,
    position: PositionType,
    label: str,
    on_click: Callable[[Event], Any],
) -> Button:
    """Creates and returns a matplotlib button."""
    ax = fig.add_axes(position)
    button = Button(ax, label)
    button.on_clicked(on_click)
    return button


def pick_cest(args: Namespace) -> None:
    """Pick peak positions in CEST profiles."""
    experiments = initialize_cest(args)
    sw = None

    for experiment in experiments:
        if not is_cest_experiment(experiment):
            sys.exit(
                f"\nError: '{experiment.name}' experiment not supported. "
                "The command 'chemex pick_cest' only works with CEST experiments.\n",
            )
        if experiment.name.startswith(("dcest", "coscest")):
            sw = 4.0

    fig, axes = create_figure()

    callback = Buttons(fig, axes, experiments, args.output, sw)

    if sw is not None:
        # For the Slider to remain responsive you must keep a reference to it.
        _slider = create_slider(
            fig,
            POSITION_SCALE,
            "Scale",
            (1.0, 10.0),
            sw,
            callback.set_sw,
        )

    # For the button to remain responsive you must keep a reference to it.
    _button1 = create_button(fig, POSITION_PREVIOUS, "Previous", callback.previous)
    _button2 = create_button(fig, POSITION_NEXT, "Next", callback.next)
    _button3 = create_button(fig, POSITION_SWAP, "Swap", callback.swap)
    _button4 = create_button(fig, POSITION_CLEAR, "Clear", callback.clear)

    fig.canvas.mpl_connect("button_press_event", callback.set_cs)

    plt.show()
