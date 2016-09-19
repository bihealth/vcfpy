# -*- coding: utf-8 -*-
"""Test that the public API classes are available from the ``vcfpy`` package
"""

import vcfpy


def test_from_exceptions():
    assert vcfpy.VCFPyException
    assert vcfpy.InvalidHeaderException
    assert vcfpy.InvalidRecordException
    assert vcfpy.IncorrectVCFFormat
    assert vcfpy.HeaderNotFound


def test_from_header():
    assert vcfpy.Header
    assert vcfpy.HeaderLine
    assert vcfpy.SimpleHeaderFile
    assert vcfpy.SimpleHeaderFile
    assert vcfpy.ContigHeaderLine
    assert vcfpy.FilterHeaderLine
    assert vcfpy.CompoundHeaderLine
    assert vcfpy.InfoHeaderLine
    assert vcfpy.FormatHeaderLine
    assert vcfpy.SamplesInfos


def test_from_record():
    assert vcfpy.Record
    assert vcfpy.Call
    assert vcfpy.AltRecord
    assert vcfpy.Substitution
    assert vcfpy.SV
    assert vcfpy.BreakEnd
    assert vcfpy.SymbolicAllele


def test_from_reader():
    assert vcfpy.Reader


def test_from_writer():
    assert vcfpy.Writer
