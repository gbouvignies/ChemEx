# ChemEx: NMR Chemical Exchange Analysis Tool

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

## Table of Contents

- [ChemEx: NMR Chemical Exchange Analysis Tool](#chemex-nmr-chemical-exchange-analysis-tool)
  - [Table of Contents](#table-of-contents)
  - [About ChemEx](#about-chemex)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
    - [Quick Start with uv (Recommended)](#quick-start-with-uv-recommended)
    - [Using pip with venv](#using-pip-with-venv)
    - [Using pip (global)](#using-pip-global)
    - [From source](#from-source)
    - [Using conda](#using-conda)
  - [Performance Optimization](#performance-optimization)
  - [Contributing](#contributing)
  - [Support and Documentation](#support-and-documentation)
  - [License](#license)
<!-- -   [Citing ChemEx](#citing-chemex) -->

## About ChemEx

ChemEx is an advanced, open-source software specifically designed for analyzing NMR experimental data to characterize chemical exchange processes. Ideal for researchers and scientists in the field of biochemistry and molecular biology, ChemEx aids in the analysis of NMR experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

Before installing ChemEx, ensure you have **Python 3.13** installed on your system.

> **Note**: ChemEx requires Python 3.13. Python 3.14 was recently released and is being tested for compatibility, but **Python 3.13 is recommended** for production use until the scientific Python ecosystem fully adopts 3.14.

## Installation

ChemEx offers several installation methods to suit your specific setup:

### Quick Start with uv (Recommended)

[uv](https://docs.astral.sh/uv/) is a fast Python package and project manager. If you don't have it installed:

```shell
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

The fastest way to try ChemEx without installation:

```shell
uvx chemex --help
```

Or install it as a tool:

```shell
uv tool install chemex
chemex --help
```

### Using pip with venv

Create an isolated environment and install ChemEx:

```shell
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install chemex
```

### Using pip (global)

```shell
pip install chemex
```

### From source

```shell
pip install git+https://github.com/gbouvignies/ChemEx.git
```

### Using conda

If you prefer conda/mamba:

```shell
conda create -n chemex python=3.13
conda activate chemex
conda config --env --add channels conda-forge
conda install chemex
```

## Performance Optimization

ChemEx performance depends on the underlying numerical libraries (NumPy and SciPy). The default installation provides good performance for most users:

- **pip** (PyPI wheels): Uses OpenBLAS on Linux/Windows, or Apple's Accelerate framework on macOS
- **conda-forge**: Uses OpenBLAS as the BLAS/LAPACK backend
- **Anaconda** (defaults channel): Uses Intel® MKL, which can provide better performance for some operations
- **Intel® Distribution for Python**: Also uses Intel® MKL

For most use cases, the default pip or conda-forge installation is sufficient. If you need maximum performance and are doing intensive numerical computations, consider using Anaconda's defaults channel or Intel's Python distribution.

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
