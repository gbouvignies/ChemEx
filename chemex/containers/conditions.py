from __future__ import annotations

import sys
from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace

import jsonschema as js

import chemex.helper as ch


@dataclass(frozen=True, eq=True)
class Conditions:
    h_larmor_frq: float | None = None
    temperature: float | None = None
    p_total: float | None = None
    l_total: float | None = None
    d2o: float | None = None
    label: list[str] = field(default_factory=list)

    def __post_init__(self):
        for name in ("h_larmor_frq", "temperature", "p_total", "l_total", "d2o"):
            value = getattr(self, name)
            if value is not None:
                super().__setattr__(name, float(value))
        super().__setattr__("label", tuple(self.label))

    @classmethod
    def from_dict(cls, dict_):
        return cls(
            dict_.get("h_larmor_frq"),
            dict_.get("temperature"),
            dict_.get("p_total"),
            dict_.get("l_total"),
            dict_.get("d2o"),
            dict_.get("label", []),
        )

    def rounded(self) -> Conditions:
        h_larmor_frq = round(self.h_larmor_frq, 1) if self.h_larmor_frq else None
        temperature = round(self.temperature, 1) if self.temperature else None
        return replace(self, h_larmor_frq=h_larmor_frq, temperature=temperature)

    def match(self, other: Conditions) -> bool:
        return self == self & other

    def __and__(self, other: Conditions) -> Conditions:
        self_dict = asdict(self)
        other_dict = asdict(other)
        intersection = {
            key: value for key, value in self_dict.items() if other_dict[key] == value
        }
        return Conditions.from_dict(intersection)


def parse_conditions(config):
    model = config["model"]
    _schema = {
        "type": "object",
        "properties": {
            "conditions": {
                "type": "object",
                "properties": {
                    "h_larmor_frq": {"type": "number"},
                    "temperature": {"type": "number"},
                    "p_total": {"type": "number"},
                    "l_total": {"type": "number"},
                    "d2o": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                    "label": {
                        "type": "array",
                        "contains": {
                            "type": "string",
                            "enum": [
                                "1H",
                                "2H",
                                "13C",
                                "15N",
                                "1h",
                                "2h",
                                "13c",
                                "15n",
                            ],
                        },
                    },
                },
                "required": ["h_larmor_frq"],
            }
        },
        "required": ["conditions"],
    }

    if "binding" in model:
        _schema["properties"]["conditions"]["required"].extend(["p_total", "l_total"])
        _schema["properties"]["conditions"]["dependencies"] = {
            "p_total": ["l_total"],
            "l_total": ["p_total"],
        }

    if "hd" in model:
        _schema["properties"]["conditions"]["required"].append("d2o")

    if "eyring" in model:
        _schema["properties"]["conditions"]["required"].append("temperature")

    try:
        ch.validate(config, _schema)

    except js.ValidationError as e:
        filename = config["filename"]
        if len(e.path) == 1:
            sys.exit(
                f"\nerror: The experiment file '{filename}' has no section "
                f"'[conditions]'."
            )
        else:
            sys.exit(
                f"\nerror: The model '{model}' requires the condition '{e.instance}' "
                f"to be defined in the experiment file. '{e.instance}' is missing from"
                f"'{filename}'."
            )

    else:
        return Conditions(**config["conditions"])
