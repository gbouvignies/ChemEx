from collections.abc import Iterable
from importlib.resources import Package
from importlib.resources import read_text
from importlib.resources import Resource

from rich.markdown import Markdown

from chemex.messages import console


class ExperimentDescriptions:

    experiment_descriptions: dict[str, str] = {}

    def register(self, type: str, package: Package, filename: Resource) -> None:
        self.experiment_descriptions[type] = read_text(package, filename)

    def first_line(self, type: str) -> str:
        return self.experiment_descriptions[type].partition("\n")[0].strip(" #")

    def print(self, type: str) -> None:
        console.print(Markdown(self.experiment_descriptions[type]), width=80)

    def experiment_names(self) -> Iterable[str]:
        return self.experiment_descriptions.keys()


descriptions = ExperimentDescriptions()
