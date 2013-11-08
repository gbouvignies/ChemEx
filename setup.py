#!/usr/bin/env python

import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

import chemex

scripts = [
    'bin/chemex_fit.py',
    'bin/chemex_mc.py',
]

setup(
    name='chemex',
    version=chemex.__version__,
    description='Program to fit chemical exchange induced shift and relaxation data.',
    long_description=open('README.md').read(),
    license='GPLv3',
    author='Guillaume Bouvignies',
    author_email='gbouvignies@gmail.com',
    packages=find_packages(),
    scripts=scripts,
    install_requires=[
        "numpy >= 1.6.0",
        "scipy >= 0.9.0",
        "matplotlib >= 1.0.0",
    ],
)
