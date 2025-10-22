---
sidebar_position: 1
---

# Welcome to ChemEx

## Quick Overview

ChemEx is a powerful tool for analyzing NMR experimental data to characterize chemical exchange processes. This guide will help you get started with ChemEx, from installation to initial usage. ChemEx is designed to support experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

Before installing ChemEx, ensure Python 3.13 (recommended) or later is installed on your system. It's recommended to create an isolated environment for a clean and conflict-free setup.

## Installation Options {#installation}

ChemEx can be installed using various methods. Choose the one that best suits your workflow.

### Quick Start with uv (Recommended)

[uv](https://docs.astral.sh/uv/) is a fast Python package and project manager written in Rust. It's significantly faster than pip and handles virtual environments automatically.

**Installing uv:**

```shell
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Using ChemEx with uv:**

The fastest way to try ChemEx without installation:

```shell
uvx chemex --help
```

Or install it as a tool for repeated use:

```shell
uv tool install chemex
chemex --help
```

### Using pip with a Virtual Environment

The standard Python approach for an isolated installation:

1. **Create and activate a virtual environment:**

   ```shell
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

2. **Install ChemEx using pip:**

   ```shell
   pip install chemex
   ```

### Using pip (global installation)

For a system-wide installation:

```shell
pip install chemex
```

### Using mamba or conda

If you prefer conda environments:

1. **Create and activate a new environment:**

   ```shell
   conda create -n chemex python=3.13
   conda activate chemex
   ```

2. **Install ChemEx via conda-forge channel:**

   ```shell
   conda install -c conda-forge chemex
   ```

   > **Tip**: For faster package management, consider using [mamba](https://mamba.readthedocs.io/) instead of conda: `mamba install -c conda-forge chemex`

### Installing from GitHub

To install the latest development version of ChemEx directly from GitHub, use the following command:

```shell
pip install git+https://github.com/gbouvignies/ChemEx.git
```

This method provides access to the latest features and updates that may not yet be included in the pip or conda-forge releases.
