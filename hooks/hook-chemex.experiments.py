from PyInstaller.utils.hooks import collect_submodules, collect_data_files

datas = collect_data_files("chemex.experiments.configs")
hiddenimports = collect_submodules("chemex.experiments")