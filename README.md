# ChemEx: NMR Chemical Exchange Analysis Tool

[![Lint: Ruff](https://img.shields.io/badge/lint-Ruff-D7FF64.svg?logo=ruff)](https://docs.astral.sh/ruff/)

## Table of Contents

- [ChemEx: NMR Chemical Exchange Analysis Tool](#chemex-nmr-chemical-exchange-analysis-tool)
    - [Table of Contents](#table-of-contents)
    - [About ChemEx](#about-chemex)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)
        - [Recommended installation](#recommended-installation)
        - [Try ChemEx without installing it](#try-chemex-without-installing-it)
        - [Updating ChemEx](#updating-chemex)
        - [Reproducible, version-pinned installation](#reproducible-version-pinned-installation)
        - [Alternative: pip](#alternative-pip)
        - [Conda packages](#conda-packages)
    - [Contributing](#contributing)
    - [Support and Documentation](#support-and-documentation)
    - [License](#license)

<!-- -   [Citing ChemEx](#citing-chemex) -->

## About ChemEx

ChemEx is an advanced, open-source software specifically designed for analyzing NMR experimental data to characterize chemical exchange processes. Ideal for researchers and scientists in the field of biochemistry and molecular biology, ChemEx aids in the analysis of NMR experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

ChemEx requires **Python 3.13 or later**. The recommended installer, uv, can
download and manage a compatible Python interpreter when needed.

## Installation

[PyPI](https://pypi.org/project/chemex/) is the authoritative Python package
distribution for ChemEx.

### Recommended installation

[uv](https://docs.astral.sh/uv/) installs ChemEx as an application in an
isolated environment. On macOS, install uv with Homebrew:

```shell
brew install uv
```

On Linux or Windows, follow Astral's current
[uv installation instructions](https://docs.astral.sh/uv/getting-started/installation/).

Then install ChemEx from PyPI and verify that it starts:

```shell
uv tool install --python 3.13 chemex
chemex --version
```

### Try ChemEx without installing it

```shell
uvx --python 3.13 chemex --help
```

### Updating ChemEx

Update an unpinned tool installation with:

```shell
uv tool upgrade chemex
```

### Reproducible, version-pinned installation

To install a specific release:

```shell
uv tool install --python 3.13 "chemex==2026.09.0"
```

### Alternative: pip

If you need conventional Python tooling, use Python 3.13 or later to install
ChemEx from PyPI inside a virtual environment:

```shell
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
python -m pip install chemex
```

### Conda packages

The historical conda-forge package is no longer maintained by the ChemEx
project and may be outdated. Install the current PyPI release with uv or pip.

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
