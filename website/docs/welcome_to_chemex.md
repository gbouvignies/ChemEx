---
sidebar_position: 1
---

# Welcome to ChemEx

## Quick Overview

ChemEx is a powerful tool for analyzing NMR experimental data to characterize chemical exchange processes. This guide will help you get started with ChemEx, from installation to initial usage. ChemEx is designed to support experiments like Carr-Purcell-Meiboom-Gill (CPMG) relaxation dispersion and Chemical Exchange Saturation Transfer (CEST).

## Prerequisites

Before installing ChemEx, ensure Python 3.11 or 3.12 is installed on your system. It’s recommended to create an isolated environment for a clean and conflict-free setup.

## Installation Options {#installation}

ChemEx can be installed through **pip** or **mamba/conda** from the **conda-forge** channel. Choose the method that best suits your setup.

### Using pip with a Virtual Environment

1. **Create and activate a virtual environment:**

   ```shell
   # Create a virtual environment (replace "chemex_env" with your preferred name)
   python -m venv chemex_env
   source chemex_env/bin/activate  # On Windows, use `chemex_env\Scripts\activate`
   ```

2. **Install ChemEx using pip:**

   ```shell
   pip install chemex
   ```

### Using mamba or conda

For users who prefer `conda` environments, **mamba** is a faster alternative to `conda` for package management but is optional.

1. **Install mamba (optional but recommended):**

   ```shell
   conda install mamba -n base -c conda-forge
   ```

2. **Create and activate a new environment:**

   ```shell
   mamba create -n chemex python=3.11
   conda activate chemex
   ```

3. **Install ChemEx via conda-forge channel:**

   ```shell
   mamba install -c conda-forge chemex
   ```

### Installing from GitHub

To install the latest development version of ChemEx directly from GitHub, use the following command:

```shell
pip install git+https://github.com/gbouvignies/ChemEx.git
```

This method provides access to the latest features and updates that may not yet be included in the pip or conda-forge releases.
