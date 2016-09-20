.. _api_header:

======
Header
======

.. contents::

vcfpy.OrderedDict
-----------------

Convenience export of ``OrderedDict``.
When available, the ``cyordereddict``, a Cython-reimplementation of ``OrderedDict`` is used for Python before 3.5 (from 3.5, Python ships with a fast, C implementation of ``OrderedDict``).

.. autoclass:: vcfpy.OrderedDict
    :members:

vcfpy.Header
------------

.. autoclass:: vcfpy.Header
    :members:

vcfpy.HeaderLine
----------------

.. autoclass:: vcfpy.HeaderLine
    :members:

vcfpy.header_without_lines
--------------------------

.. autofunction:: vcfpy.header_without_lines

vcfpy.SimpleHeaderFile
----------------------

.. autoclass:: vcfpy.SimpleHeaderFile
    :members:

vcfpy.ContigHeaderLine
----------------------

.. autoclass:: vcfpy.ContigHeaderLine
    :members:

vcfpy.FilterHeaderLine
----------------------

.. autoclass:: vcfpy.FilterHeaderLine
    :members:

vcfpy.CompoundHeaderLine
------------------------

.. autoclass:: vcfpy.CompoundHeaderLine
    :members:

vcfpy.InfoHeaderLine
--------------------

.. autoclass:: vcfpy.InfoHeaderLine
    :members:

vcfpy.FormatHeaderLine
----------------------

.. autoclass:: vcfpy.FormatHeaderLine
    :members:

vcfpy.FieldInfo
---------------

.. autoclass:: vcfpy.FieldInfo
    :members:

vcfpy.SamplesInfos
------------------

.. autoclass:: vcfpy.SamplesInfos
    :members:
