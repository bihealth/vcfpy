#!/usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import chain
import sys

from setuptools import setup
import pip
from pip.req import parse_requirements

import versioneer

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

base_reqs = parse_requirements(
    'requirements.txt', session=pip.download.PipSession())

# Add cyordereddict for Python <=3.5 for performance boost
if sys.version_info[:2] < (3, 6):
    pre36_reqs = parse_requirements(
        'requirements/pre36.txt', session=pip.download.PipSession())
else:
    pre36_reqs = []

requirements = [str(ir.req) for ir in chain(base_reqs, pre36_reqs)]

test_requirements = [
    str(ir.req)
    for ir in parse_requirements(
        'requirements/test.txt', session=pip.download.PipSession())
]

setup(
    name='vcfpy',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=(
        'Python 3 VCF library with good support for both reading and writing'),
    long_description=readme + '\n\n' + history,
    author="Manuel Holtgrewe",
    author_email='manuel.holtgrewe@bihealth.de',
    url='https://github.com/bihealth/vcfpy',
    packages=[
        'vcfpy',
    ],
    package_dir={'vcfpy':
                 'vcfpy'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='vcfpy',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
