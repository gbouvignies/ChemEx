from importlib.resources import Package
from importlib.resources import read_text
from importlib.resources import Resource
from pathlib import Path

from rich.markdown import Markdown

from chemex.messages import console


def toml_to_md(text: str) -> str:
    return f"```toml\n{text}\n```"


class ExperimentTomlConfigurations:

    experiment_configurations: dict[str, str] = {}

    def register(self, type: str, package: Package, filename: Resource) -> None:
        self.experiment_configurations[type] = read_text(package, filename)

    def print(self, type: str) -> None:
        console.print(Markdown(toml_to_md(self.experiment_configurations[type])))

    def write(self, type: str, filename: Path) -> None:
        filename.write_text(self.experiment_configurations[type])


configurations = ExperimentTomlConfigurations()
