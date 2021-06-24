import dataclasses as dc
import sys
from typing import List
from typing import Optional

import jsonschema as js

import chemex.helper as ch


@dc.dataclass(frozen=True, eq=True)
class Conditions:
    h_larmor_frq: Optional[float] = None
    temperature: Optional[float] = None
    p_total: Optional[float] = None
    l_total: Optional[float] = None
    d2o: Optional[float] = None
    label: List[str] = dc.field(default_factory=list)

    def __post_init__(self):
        for field in ("h_larmor_frq", "temperature", "p_total", "l_total", "d2o"):
            value = getattr(self, field)
            if value is not None:
                super().__setattr__(field, float(value))
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

    def rounded(self):
        h_larmor_frq = round(self.h_larmor_frq, 1) if self.h_larmor_frq else None
        temperature = round(self.temperature, 1) if self.temperature else None
        return dc.replace(self, h_larmor_frq=h_larmor_frq, temperature=temperature)

    def __and__(self, other):
        both = [
            value1 if value1 == value2 else None
            for value1, value2 in zip(dc.astuple(self), dc.astuple(other))
        ]
        return Conditions(*both)


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
                            "enum": ["1h", "2h", "13c", "15n"],
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
