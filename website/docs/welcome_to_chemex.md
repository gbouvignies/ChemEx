---
sidebar_position: 1
description: Install ChemEx and learn the first steps for analyzing NMR chemical exchange data with CPMG relaxation dispersion and CEST experiments.
---

# Welcome to ChemEx

## Quick Overview

ChemEx is a powerful tool for analyzing NMR experimental data to characterize chemical exchange processes. This guide will help you get started with ChemEx, from installation to initial usage. ChemEx is designed to support experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

ChemEx requires Python 3.13 or later. The recommended installer, uv, can
download and manage a compatible Python interpreter when needed.

## Installation {#installation}

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

uv creates a dedicated environment for ChemEx and automatically downloads a
Python 3.13 interpreter if a suitable interpreter is not already available.

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
