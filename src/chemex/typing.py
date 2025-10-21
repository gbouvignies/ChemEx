from typing import Any

import numpy.typing as npt

# Use a pragmatic, mypy-friendly alias for numpy arrays used across the codebase.
# Many arrays in this project can be float, complex, or integer and have varying
# dtypes/shapes. Using NDArray[Any] keeps typing useful without fighting
# numpy's intricate dtype generics in many places.
Array = npt.NDArray[Any]
