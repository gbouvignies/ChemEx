from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.utils.hooks import collect_submodules


datas = collect_data_files("chemex.experiments.configs")
hiddenimports = collect_submodules("chemex.experiments")
