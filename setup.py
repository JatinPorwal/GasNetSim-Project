#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2022.
#     Developed by Yifei Lu
#     Last change on 7/26/22, 1:49 PM
#     Last change by yifei
#    *****************************************************************************

import distutils.command.install as orig
import inspect
import os
import re
import subprocess
import sys

import setuptools
from setuptools import setup


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


def pip_install(package: str):
    """Installs a package manually with pip"""
    return subprocess.check_call(["pip", "install", package])


# def extract_version(version_file_name):
#     """Extracts the version from a python file.
#     The statement setting the __version__ variable must not be indented. Comments after that
#     statement are allowed.
#     """
#     regex = re.compile(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]\s*(#.*)?$")
#     with open(version_file_name, "r") as version_file:
#         lines = version_file.read().splitlines()
#     for line in reversed(lines):
#         version_match = regex.match(line)
#         if version_match:
#             return version_match.group(1)
#     else:
#         raise RuntimeError("Unable to find version string.")


def main():
    """Main-function which is called if this file is executed as script."""
    with open("README.md", "r") as file:
        long_description = file.read()

    version = "0.1.0"

    print(version)

    return setup(name="GasNetSim",
                 version=version,
                 author="IEK-10: Energy Systems Engineering, Forschungszentrum Jülich GmbH",
                 author_email="yifei.lu@fz-juelich.de",
                 description="A tool for gas network steady-state simulation.",
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 python_requires=">=3.7",
                 install_requires=["numpy",
                                   "matplotlib",
                                   "scipy",
                                   "pandas",
                                   "pytest",
                                   "fluids",
                                   "pint",
                                   "setuptools",
                                   "requests",
                                   "pyparsing~=3.0.7",
                                   "cantera~=2.5.1"],
                 classifiers=["Programming Language :: Python :: 3.7",
                              "License :: OSI Approved :: Mozilla Public License 2.0",
                              "Operating System :: OS Independent"],
                 url="https://jugit.fz-juelich.de/iek-10/public/simulation/gasnetsim",
                 project_urls={
                     "Bug Tracker": "https://jugit.fz-juelich.de/iek-10/public/simulation/gasnetsim/-/issues",
                 },
                 package_dir={"": "GasNetSim"},
                 packages=setuptools.find_packages(where="GasNetSim"),
                 )


if __name__ == "__main__":
    main()