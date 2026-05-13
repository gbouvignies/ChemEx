import numpy as np
import pytest
from pydantic import ValidationError

from chemex.containers.data import Data


def test_data_rejects_derived_array_inputs() -> None:
    with pytest.raises(ValidationError, match="Derived array fields"):
        Data(
            exp=np.array([1.0]),
            err=np.array([0.1]),
            calc=np.array([1.0]),
        )


def test_data_rejects_mismatched_exp_and_err_shapes() -> None:
    with pytest.raises(ValidationError, match="'err' must match 'exp' shape"):
        Data(exp=np.array([1.0, 2.0]), err=np.array([0.1]))


def test_data_accepts_two_sided_error_bars() -> None:
    data = Data(
        exp=np.array([1.0, 2.0]),
        err=np.array([[0.1, 0.2], [0.3, 0.4]]),
    )

    assert data.err.shape == (2, 2)


def test_data_rejects_mismatched_metadata_length() -> None:
    with pytest.raises(ValidationError, match="'metadata' must be empty"):
        Data(
            exp=np.array([1.0, 2.0]),
            err=np.array([0.1, 0.2]),
            metadata=np.array([10.0]),
        )
