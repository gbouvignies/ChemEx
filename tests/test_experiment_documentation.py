from __future__ import annotations

import re
import tomllib
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import TypeIs

import pytest
from pydantic import BaseModel

from chemex.experiments.experiment_types import (
    ExperimentType,
    registered_experiment_types,
)
from chemex.experiments.loader import register_experiments

ROOT = Path(__file__).parents[1]
EXPERIMENT_DOCS = ROOT / "website" / "docs" / "experiments"
EXPERIMENT_CONFIG = re.compile(
    r'^```toml title="experiment\.toml"\n(?P<config>.*?)^```$',
    flags=re.DOTALL | re.MULTILINE,
)
EXPERIMENT_PAGES = tuple(sorted(EXPERIMENT_DOCS.rglob("*.md")))

type ConfigurationTable = dict[str, object]


def _is_configuration_table(value: object) -> TypeIs[ConfigurationTable]:
    return isinstance(value, dict) and all(isinstance(key, str) for key in value)


def _find_unknown_keys(
    raw: ConfigurationTable,
    validated: BaseModel,
    path: tuple[str, ...] = (),
) -> Iterator[str]:
    fields = type(validated).model_fields
    for key, raw_value in raw.items():
        normalized_key = key.lower()
        key_path = (*path, key)
        if normalized_key not in fields:
            yield ".".join(key_path)
            continue
        validated_value = getattr(validated, normalized_key)
        if _is_configuration_table(raw_value) and isinstance(
            validated_value,
            BaseModel,
        ):
            yield from _find_unknown_keys(raw_value, validated_value, key_path)


@pytest.fixture(scope="module")
def experiment_types() -> Mapping[str, ExperimentType]:
    register_experiments()
    return registered_experiment_types()


def test_experiment_documentation_discovery_is_not_empty() -> None:
    assert EXPERIMENT_PAGES, f"No experiment pages found under {EXPERIMENT_DOCS}"


@pytest.mark.parametrize("page", EXPERIMENT_PAGES, ids=lambda page: page.stem)
def test_experiment_page_configuration_matches_current_schema(
    page: Path,
    experiment_types: Mapping[str, ExperimentType],
) -> None:
    matches = tuple(EXPERIMENT_CONFIG.finditer(page.read_text(encoding="utf-8")))
    assert len(matches) == 1, (
        f"{page.relative_to(ROOT)} must contain exactly one copyable "
        '```toml title="experiment.toml" configuration'
    )

    raw_config = tomllib.loads(matches[0].group("config"))
    experiment_name = raw_config["experiment"]["name"]
    assert experiment_name == page.stem
    assert experiment_name in experiment_types

    config = experiment_types[experiment_name].config_type.model_validate(raw_config)
    assert not tuple(_find_unknown_keys(raw_config, config))


def test_unknown_documented_keys_are_detected_when_pydantic_ignores_them(
    experiment_types: Mapping[str, ExperimentType],
) -> None:
    raw_config: ConfigurationTable = {
        "experiment": {"name": "shift_15n_sqmq"},
        "conditions": {},
        "data": {"shifts": "shifts.dat"},
    }
    config = experiment_types["shift_15n_sqmq"].config_type.model_validate(raw_config)

    assert tuple(_find_unknown_keys(raw_config, config)) == ("data.shifts",)
