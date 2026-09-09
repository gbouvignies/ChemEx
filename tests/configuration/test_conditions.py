import pytest
from pydantic import ValidationError

from chemex.configuration.conditions import Conditions, ConditionsWithValidations
from chemex.models.model import ModelSpec
from chemex.parameters.name import ParamName


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


@pytest.mark.parametrize(
    "temperature",
    (-273.15, -273.15000000000003, float("nan"), float("inf"), float("-inf")),
)
def test_eyring_model_rejects_temperature_outside_physical_domain(
    temperature: float,
) -> None:
    with pytest.raises(
        ValidationError,
        match="Eyring temperature must be finite and above absolute zero",
    ):
        ConditionsWithValidations.model_validate(
            {"temperature": temperature},
            context={"model": ModelSpec(name="2st_eyring")},
        )


def test_exact_zero_celsius_remains_part_of_condition_identity() -> None:
    zero = Conditions(temperature=0.0).rounded()
    absent = Conditions().rounded()

    assert zero.temperature == 0.0
    assert zero != absent
    assert ParamName("kab", conditions=zero).id_ == "__KAB__0_0C"
    assert (
        ParamName("kab", conditions=zero).id_ != ParamName("kab", conditions=absent).id_
    )
