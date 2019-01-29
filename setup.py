#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import chain
import os.path
import sys

from setuptools import setup

import versioneer


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

base_reqs = parse_requirements("requirements.txt")

# Add cyordereddict for Python <=3.5 for performance boost
if sys.version_info[:2] < (3, 6):
    pre36_reqs = parse_requirements("requirements/pre36.txt")
else:
    pre36_reqs = []

requirements = chain(base_reqs, pre36_reqs)

test_requirements = parse_requirements("requirements/test.txt")

setup(
    name="vcfpy",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=("Python 3 VCF library with good support for both reading and writing"),
    long_description=readme + "\n\n" + history,
    author="Manuel Holtgrewe",
    author_email="manuel.holtgrewe@bihealth.de",
    url="https://github.com/bihealth/vcfpy",
    packages=["vcfpy"],
    package_dir={"vcfpy": "vcfpy"},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords="vcfpy",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    test_suite="tests",
    tests_require=test_requirements,
)
