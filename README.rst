=====
VCFPy
=====


.. image:: https://img.shields.io/pypi/v/vcfpy.svg
        :target: https://pypi.python.org/pypi/vcfpy

.. image:: https://img.shields.io/travis/bihealth/vcfpy.svg
        :target: https://travis-ci.org/bihealth/vcfpy

.. image:: https://readthedocs.org/projects/vcfpy/badge/?version=latest
        :target: https://vcfpy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://api.codacy.com/project/badge/Grade/cfe741307ec34e8fb90dfe37e84a2519
        :target: https://www.codacy.com/app/manuel-holtgrewe/vcfpy?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/vcfpy&amp;utm_campaign=Badge_Grade
        :alt: Codacy Analysis

.. image:: https://api.codacy.com/project/badge/Coverage/cfe741307ec34e8fb90dfe37e84a2519
        :alt: Codacy Coverage
        :target: https://www.codacy.com/app/manuel-holtgrewe/vcfpy?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/vcfpy&amp;utm_campaign=Badge_Coverage

.. image:: https://landscape.io/github/bihealth/vcfpy/master/landscape.svg?style=flat
        :alt: Landscape Health
        :target: https://landscape.io/github/bihealth/vcfpy

.. image:: http://joss.theoj.org/papers/edae85d90ea8a49843dbaaa109e47cba/status.svg
        :alt: Publication in The Journal of Open Source Software
        :target: http://joss.theoj.org/papers/10.21105/joss.00085

Python 3 VCF library with good support for both reading and writing

* Free software: MIT license
* Documentation: https://vcfpy.readthedocs.io.


Features
--------

- Support for reading and writing VCF v4.3
- Interface to ``INFO`` and ``FORMAT`` fields is based on ``OrderedDict`` allows for easier modification than PyVCF (also I find this more pythonic)
- Read (and jump in) and write BGZF files just using ``vcfpy``

Why another VCF parser for Python!
----------------------------------

I've been using PyVCF with quite some success in the past.
However, the main bottleneck of PyVCF is when you want to modify the per-sample genotype information.
There are some issues in the tracker of PyVCF but none of them can really be considered solved.
I tried several hours to solve these problems within PyVCF but this never got far or towards a complete rewrite...

For this reason, VCFPy was born and here it is!

What's the State?
-----------------

VCFPy is the result of two full days of development plus some maintenance work later now (right now).
I'm using it in several projects but it is not as battle-tested as PyVCF.

Why Python 3 Only?
------------------

As I'm only using Python 3 code, I see no advantage in carrying around support for legacy Python 2 and maintaining it.
At a later point when VCFPy is known to be stable, Python 2 support might be added if someone contributes a pull request.
