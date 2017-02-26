# -*- coding: utf-8 -*-
"""Reading of VCF header from plain and bgzip-ed file
"""

import os

from vcfpy import reader

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


def test_read_text():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/from_vcf43.vcf')
    r = reader.Reader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    EXPECTED = "HeaderLine('fileformat', 'VCFv4.3')"
    assert str(r.header.lines[0]) == EXPECTED
    EXPECTED = (
        "FormatHeaderLine('FORMAT', '<ID=HQ,Number=2,Type=Integer,"
        "Description=\"Haplotype Quality\">', OrderedDict([('ID', 'HQ'), "
        "('Number', 2), ('Type', 'Integer'), ('Description', "
        "'Haplotype Quality')]))")
    assert str(r.header.lines[-1]) == EXPECTED
    assert r.header.samples
    assert r.header.samples.names == ['NA00001', 'NA00002', 'NA00003']


def test_read_bgzip():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/from_vcf43.vcf.gz')
    r = reader.Reader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    EXPECTED = "HeaderLine('fileformat', 'VCFv4.3')"
    assert str(r.header.lines[0]) == EXPECTED
    EXPECTED = (
        "FormatHeaderLine('FORMAT', '<ID=HQ,Number=2,Type=Integer,"
        "Description=\"Haplotype Quality\">', OrderedDict([('ID', 'HQ'), "
        "('Number', 2), ('Type', 'Integer'), ('Description', "
        "'Haplotype Quality')]))")
    assert str(r.header.lines[-1]) == EXPECTED
    assert r.header.samples
    assert r.header.samples.names == ['NA00001', 'NA00002', 'NA00003']
