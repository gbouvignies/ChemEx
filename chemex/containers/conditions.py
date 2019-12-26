import dataclasses as dc
import sys
from typing import List
from typing import Optional

import jsonschema as js

import chemex.helper as ch


@dc.dataclass(frozen=True, eq=True)
class Conditions:
    h_larmor_frq: float
    temperature: Optional[float] = None
    p_total: Optional[float] = None
    l_total: Optional[float] = None
    label: List[str] = dc.field(default_factory=list)

    def __post_init__(self):
        super().__setattr__("label", tuple(self.label))


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
                    "label": {
                        "type": "array",
                        "contains": {
                            "type": "string",
                            "enum": ["1H", "2H", "13C", "15N"],
                        },
                    },
                },
                "dependencies": {"p_total": ["l_total"], "l_total": ["p_total"]},
                "required": ["h_larmor_frq"],
            }
        },
        "required": ["conditions"],
    }
    if "binding" in model:
        _schema["properties"]["conditions"]["required"].extend(["p_total", "l_total"])

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
