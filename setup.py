from ez_setup import use_setuptools

use_setuptools()

import sys

if sys.version_info[0] == 2 and sys.version_info[1] < 7:
    print "Sorry, Python < 2.7 is not supported"
    exit()

from setuptools import setup, find_packages

with open('chemex/version.py') as f:
    exec (f.read())  # read and set the variable __version__

description = ('Program to fit chemical exchange induced shift and '
               'relaxation data.')

install_requires = [
    'numpy',
    'scipy',
    'matplotlib>=1.1',
]

entry_points = {
    'console_scripts': [
        'chemex = chemex.__main__:main',
    ]
}

setup(
    name='chemex',
    version=__version__,
    description=description,
    long_description=open('README.md').read(),
    license='BSD 3-Clause',
    author='Guillaume Bouvignies',
    author_email='gbouvignies@gmail.com',
    packages=find_packages(),
    install_requires=install_requires,
    entry_points=entry_points
)
