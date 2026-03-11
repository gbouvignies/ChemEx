import pytest
from pydantic import ValidationError

from chemex.configuration.conditions import ConditionsWithValidations
from chemex.models.model import ModelSpec


def test_hd_model_requires_d2o() -> None:
    with pytest.raises(
        ValidationError, match='To use the "hd" model, d2o must be provided'
    ):
        ConditionsWithValidations.model_validate(
            {},
            context={"model": ModelSpec(name="2st_hd")},
        )


def test_eyring_model_requires_temperature() -> None:
    with pytest.raises(
        ValidationError,
        match='To use the "eyring" model, "temperature" must be provided',
    ):
        ConditionsWithValidations.model_validate(
            {},
            context={"model": ModelSpec(name="2st_eyring")},
        )
