from __future__ import annotations

import pytest

from chemex.models.loader import register_kinetic_settings
from chemex.models.model import ModelSpec


def setup_module() -> None:
    register_kinetic_settings()


def test_model_spec_parses_suffix_combinations() -> None:
    spec = ModelSpec.from_name("2st.rs.mf.tc")

    assert spec.name == "2st"
    assert spec.states == "ab"
    assert spec.residue_specific is True
    assert spec.model_free is True
    assert spec.temp_coef is True


def test_model_spec_rejects_unknown_suffix() -> None:
    with pytest.raises(SystemExit):
        ModelSpec.from_name("2st.xyz")
