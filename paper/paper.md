---
title: 'VCFPy: a Python 3 library with good support for both reading and writing VCF'
tags:
  - VCF
  - Python
  - software library
authors:
 - name: Manuel Holtgrewe
   orcid: 0000-0002-3051-1763
   affiliation: 1
 - name: Dieter Beule
   orcid: 0000-0002-3284-0632
   affiliation: 1
affiliations:
 - name: Berlin Institute of Health, Kapelle-Ufer 2, 10117 Berlin
   index: 1
date: 28 September 2016
bibliography: paper.bib
---

# Summary

VCF file format [@Danecek2011] is the standard file format for genetic variants, both small and structural variants.
It has broad adaption in the Bioinformatics community and is used both by most projects, software, and databases these days.

There is a number of Python libraries for processing VCF, but most focus on reading VCF and not allowing for easily creating or augmenting VCF headers and records.
For example, the most popular library PyVCF does not allow for built-in modification of the per-sample `FORMAT/*` records.
PySAM (the wrapper for htslib) does only have very limited support for modifyin VCF records at all.

VCFPy addresses these issues and provides a well-documented, easy to use, and pythonic interface to reading and writing VCF files.
It supports VCF v4.3, reading and writing of both plain-text and bgzip-compressed VCF files, as well as Tabix indices.
Further, the project is well-documented and uses automatic testing as well as static code analysis for enforcing software quality standards.

# References
