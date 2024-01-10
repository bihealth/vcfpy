[![pypi](https://img.shields.io/pypi/v/vcfpy.svg)](https://pypi.python.org/pypi/vcfpy)
[![bioconda](https://img.shields.io/conda/dn/bioconda/vcfpy.svg?label=Bioconda)](https://bioconda.github.io/recipes/vcfpy/README.html)
[![CI](https://github.com/bihealth/vcfpy/actions/workflows/main.yml/badge.svg)](https://github.com/bihealth/vcfpy/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/vcfpy/badge/?version=latest)](https://vcfpy.readthedocs.io/en/latest/?badge=latest)
[![Publication in The Journal of Open Source Software](http://joss.theoj.org/papers/edae85d90ea8a49843dbaaa109e47cba/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00085)

# VCFPy

Python 3 VCF library with good support for both reading and writing

- Free software: MIT license
- Documentation: <https://vcfpy.readthedocs.io>.

## Features

- Support for reading and writing VCF v4.3
- Interface to `INFO` and `FORMAT` fields is based on `OrderedDict` allows for easier modification than PyVCF (also I find this more pythonic)
- Read (and jump in) and write BGZF files just using `vcfpy`

## Why another VCF parser for Python!

I've been using PyVCF with quite some success in the past. However, the
main bottleneck of PyVCF is when you want to modify the per-sample
genotype information. There are some issues in the tracker of PyVCF but
none of them can really be considered solved. I tried several hours to
solve these problems within PyVCF but this never got far or towards a
complete rewrite...

For this reason, VCFPy was born and here it is!

## What's the State?

VCFPy is the result of two full days of development plus some
maintenance work later now (right now). I'm using it in several projects
but it is not as battle-tested as PyVCF.

## Why Python 3 Only?

As I'm only using Python 3 code, I see no advantage in carrying around
support for legacy Python 2 and maintaining it. At a later point when
VCFPy is known to be stable, Python 2 support might be added if someone
contributes a pull request.

