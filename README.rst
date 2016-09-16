=====
VCFPy
=====


.. image:: https://img.shields.io/pypi/v/vcfpy.svg
        :target: https://pypi.python.org/pypi/vcfpy

.. image:: https://img.shields.io/travis/holtgrewe/vcfpy.svg
        :target: https://travis-ci.org/holtgrewe/vcfpy

.. image:: https://readthedocs.org/projects/vcfpy/badge/?version=latest
        :target: https://vcfpy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/holtgrewe/vcfpy/shield.svg
     :target: https://pyup.io/repos/github/holtgrewe/vcfpy/
     :alt: Updates


Python 3 VCF parser that allows both reading and writing

* Free software: MIT license
* Documentation: https://vcfpy.readthedocs.io.

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