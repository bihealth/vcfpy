.. _manual-main:

Welcome to VCFPy's documentation!
=================================

VCFPy is a Python 3 library with good support for both reading and writing VCF files.
The documentation is split into three parts (accessible through the navigation on the left):

Installation & Getting Started
    Instructions for the installation of the module and some examples to get you started.

API Documentation
    This section contains the API documentation for the module

Project Info
    More information on the project, including the changelog, list of contributing authors, and contribution instructions.

Quick Example
-------------

.. literalinclude:: ../examples/add_filter/add_filter.py
    :language: python

Features
--------

* Support for reading and writing VCF v4.3
* Interface to ``INFO`` and ``FORMAT`` fields is based on ``OrderedDict`` allows for easier modification than PyVCF (also I find this more pythonic)
* Read (and jump in) and write BGZF files just using ``vcfpy``

Frequently Asked Questions
--------------------------

Why another Python library for VCF?
    I've been using PyVCF with quite some success in the past.
    However, the main bottleneck of PyVCF is when you want to modify the per-sample genotype information.
    There are some issues in the tracker of PyVCF but none of them can really be considered solved.
    I tried several hours to solve these problems within PyVCF but this never got far or towards a complete rewrite...

    For this reason, VCFPy was born and here it is!

Why Python 3 only?
    As I'm only using Python 3 code, I see no advantage in carrying around support for legacy Python 2 and maintaining it.
    At a later point when VCFPy is known to be stable, Python 2 support might be added if someone contributes a pull request.

What's the state?
    VCFPy is the result of two full days of development plus some maintenance work later now (right now).
    I'm using it in several projects but it is not as battle-tested as PyVCF.

What's the difference to PyVCF?
    The main difference is technical.
    Instead of using ``collections.namedtuple`` for storing the call annotation, VCFPy uses ``collections.OrderedDict``.
    This has the advantage that (1) access to optional settings is much more pythonic using ``.get(KEY, DEFAULT)`` instead of ``getattr()``.
    Further, (2) adding call annotations (``FORMAT``) fields is able without any performance penalty where for PyVCF, ``copy.deepcopy`` has to be used at some point which is very slow.
    There has not been any movement in supporting modifying ``FORMAT`` fields in PyVCF and here is a library that does this well.

What's the aim?
    The aim of the project is to provide simple yet efficient read and write access to VCF files.
    Eventually, PySAM will probably be a better choice once it has a Python wrapper for the VCF part of ``htslib``.
    However, as this is still misssing, ``VCFPy`` is a good solution for the time being.

.. toctree::
    :caption: Installation & Getting Started
    :name: getting-started
    :hidden:
    :maxdepth: 1

    installation
    getting_started
    examples
    best_practice

.. toctree::
    :caption: API Reference
    :name: api-reference
    :hidden:
    :maxdepth: 1
    :titlesonly:

    api_header
    api_io
    api_exceptions
    api_record

.. toctree::
    :caption: Project Info
    :name: project-info
    :hidden:
    :maxdepth: 1
    :titlesonly:

    contributing
    authors
    history
    license

.. Generated pages, should not appear

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
