#!/usr/bin/env python

import sys

import ez_setup
import chemex


if sys.version_info[0] == 2 and sys.version_info[1] < 7:
    print "Sorry, Python < 2.7 is not supported"
    exit()

ez_setup.use_setuptools()

from setuptools import setup, find_packages

scripts = [
    'bin/chemex_fit.py',
    'bin/chemex_mc.py',
    'bin/chemex_bootstrap.py',
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
        'numpy',
        'scipy',
        'matplotlib>=1.1'
    ],
)
