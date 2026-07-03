# setup.py
from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="screening_scalar_tensor_tov_solver",
    version="1.0.0",
    author="Joan Bachs Esteban",
    author_email="joan.bachs@tecnico.ulisboa.pt",
    description="Tolman–Oppenheimer–Volkoff solver for Scalar–Tensor theories with Screening Mechanisms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JoanBachs/Screening_ScalarTensor_TOV_Solver",
    license="GPL-3.0-or-later",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "matplotlib",
        "mpmath",
    ],
    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    keywords="white dwarfs, neutron stars, scalar-tensor gravity, TOV, screening mechanisms, modified gravity",
)