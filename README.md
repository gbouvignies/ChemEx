# ChemEx: NMR Chemical Exchange Analysis Tool

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

## Table of Contents

-   [About ChemEx](#about-chemex)
-   [Quick Overview](#quick-overview)
-   [Prerequisites](#prerequisites)
-   [Installation](#installation)
-   [Performance Optimization](#performance-optimization)
-   [Contributing](#contributing)
-   [Support and Documentation](#support-and-documentation)
-   [License](#license)
<!-- -   [Citing ChemEx](#citing-chemex) -->

## About ChemEx

ChemEx is an advanced, open-source software specifically designed for analyzing NMR experimental data to characterize chemical exchange processes. Ideal for researchers and scientists in the field of biochemistry and molecular biology, ChemEx aids in the analysis of NMR experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

Before installing ChemEx, ensure you have Python installed on your system. For beginners and for a seamless setup, we recommend using the [Anaconda Distribution](https://www.anaconda.com/distribution/), which includes Python, Numpy, and other essential scientific computing tools.

## Installation

ChemEx offers several installation methods to suit your specific setup:

### Using conda

```shell
conda create -n chemex
conda activate chemex
conda install python=3.11
conda config --env --add channels conda-forge
conda install chemex
```

### Using pip

```shell
pip install chemex
```

### From source

```shell
pip install git+https://github.com/gbouvignies/ChemEx.git
```

## Performance Optimization

For the best performance, install [NumPy](https://numpy.org) and [SciPy](https://scipy.org) with Intel® Math Kernel Library (Intel® MKL), available via [Anaconda](https://www.anaconda.com/distribution/) or the [Intel® Distribution for Python](https://software.intel.com/en-us/distribution-for-python).

## Contributing

We encourage contributions from the community. Please see our [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to make ChemEx better. For any issues or suggestions, please open an issue or a discussion on our [GitHub repository](https://github.com/gbouvignies/ChemEx).

## Support and Documentation

For additional support, tutorials, and detailed documentation, visit the [ChemEx Documentation](https://gbouvignies.github.io/ChemEx/).

## License

ChemEx is licensed under the [GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html). See the [LICENSE](LICENSE.md) file for more details.

<!-- ## Citing ChemEx

If you use ChemEx in your research, please cite it as follows: [Citation details](#). -->

---

Developed with ❤️ by the [ChemEx Contributors](https://github.com/gbouvignies/ChemEx/graphs/contributors)
