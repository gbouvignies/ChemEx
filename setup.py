# -*- coding: utf-8 -*-


"""setup.py: setuptools control."""

import re

from setuptools import setup, find_packages

version = re.search(
    '^__version__\s*=\s*["|\'](.*)["|\']',
    open('chemex/version.py').read(),
    re.M
).group(1)

with open("README.md", "rb") as f:
    long_description = f.read().decode("utf-8")


setup(
    name='chemex',

    version=version,

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
