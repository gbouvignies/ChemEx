# ChemEx

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

## Overview

ChemEx is an analysis program for chemical exchange detected by NMR. It is designed to analyze various types of NMR data, with a focus on CPMG relaxation dispersion and Chemical Exchange Saturation Transfer techniques.

## Installation

You can install `chemex` using different methods:

### Conda

The recommended way to install `chemex` is via [conda](http://conda.pydata.org):

```bash
conda install -c conda-forge chemex
```

If your version of Python is less than 3.9, you can create a separate conda environment and enforce the use of Python 3.9+:

```bash
conda create -c conda-forge -n chemex python=3.10 chemex
conda activate chemex
```

### PyPI (Python Package Index)

`chemex` is also available on the [Python Package Index](https://pypi.python.org/pypi/chemex) and can be installed using `pip`:

```bash
pip install chemex
```

### Development Version

To install the development version directly from GitHub, you can use `pip`:

```bash
pip install git+https://github.com/gbouvignies/chemex.git
```

Make sure you have Git installed and configured on your system before running this command.

## Contributing

Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or a discussion on the [GitHub repository](https://github.com/gbouvignies/chemex). We appreciate your feedback and involvement in making ChemEx better.

## License

ChemEx is released under the [MIT License](https://github.com/gbouvignies/chemex/blob/master/LICENSE). Please see the LICENSE file for more details.
