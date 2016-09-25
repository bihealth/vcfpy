#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'pysam>=0.9',  # for tabix support
]

# Add cyordereddict for Python <=3.5 for performance boost
if sys.version_info[:2] < (3, 5):
    requirements.append('cyordereddict>=1.0.0')

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='vcfpy',
    version='0.7.0',
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
