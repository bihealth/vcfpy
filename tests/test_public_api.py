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
    assert vcfpy.VCFHeader
    assert vcfpy.VCFHeaderLine
    assert vcfpy.VCFSimpleHeaderLine
    assert vcfpy.VCFSimpleHeaderLine
    assert vcfpy.VCFContigHeaderLine
    assert vcfpy.VCFFilterHeaderLine
    assert vcfpy.VCFCompoundHeaderLine
    assert vcfpy.VCFInfoHeaderLine
    assert vcfpy.VCFFormatHeaderLine
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
    assert vcfpy.VCFReader


def test_from_writer():
    assert vcfpy.VCFWriter
