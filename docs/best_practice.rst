.. _best_practice:

=============
Best Practice
=============

While not strictly part of the documentation of VCFPy, we include some notes on hints that we consider best practice when building VCF processing applications.

Keep Input Verbatim Where Possible
==================================

Try to keep the input verbatim if there is no strong reason for adjusting it.
Strong reasons include fixing ``Type`` or ``Number`` in header lines describing arrays of strings, for example.

Whenever possible, keep the header order intact.
VCFPy does this automatically for you (in contrast to PyVCF).

Prefer Soft-Filters over Hard-Filters
=====================================

**Soft**-filters mean annotating your VCF records in the ``FILTER`` column whereas **Hard**-filters mean removing records from VCF file.
In many situations, it is useful to keep around all VCF records and just annotate why they are to be dropped.
Then, in the last step, only the interesting ones are kept.

This makes tracing back easier when and why a record was removed.