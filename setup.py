from setuptools import find_packages, setup
import versioneer

with open("README.md", "rb") as f:
    long_description = f.read().decode("utf-8")

setup(
    name="chemex",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="ChemEx is an analysis program for chemical exchange detected by NMR.",
    long_description=long_description,
    author="Guillaume Bouvignies",
    author_email="gbouvignies@gmail.com",
    url="https://github.com/gbouvignies/chemex",
    license="3-Clause BSD",
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    keywords="nmr protein dynamics chemical exchange cpmg cest relaxation data fitting",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.0",
        "scipy>=1.0",
        "matplotlib>=2.0",
        "lmfit>=0.9.11",
        "asteval>=0.9.11",
    ],
    python_requires=">=3.5",
    entry_points={"console_scripts": ["chemex = chemex.chemex:main"]},
)
