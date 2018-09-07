"""chemex.__main__: executed when the chemex directory is called as script."""

import numpy as np

import chemex

with np.errstate(all="raise"):
    chemex.main()
