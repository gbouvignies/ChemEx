# -*- coding: utf-8 -*-


"""setup.py: setuptools control."""
from ez_setup import use_setuptools
use_setuptools()

use_setuptools()

from setuptools import setup, find_packages
import versioneer

with open("README.md", "rb") as f:
    long_description = f.read().decode("utf-8")


setup(
    name='chemex',

    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),

    description='Program to fit chemical exchange induced shift and relaxation data',
    long_description=long_description,

    author='Guillaume Bouvignies',
    author_email='gbouvignies@gmail.com',

    url='https://github.com/gbouvignies/chemex',

    license='3-Clause BSD',

    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: BSD License',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='nmr protein dynamics chemical exchange cpmg cest relaxation data fitting',

    packages=find_packages(),

    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'lmfit',
    ],

    entry_points={
        'console_scripts': [
            'chemex = chemex.chemex:main',
        ],
    },
)
