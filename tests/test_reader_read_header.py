# -*- coding: utf-8 -*-
"""Reading of VCF header from plain and bgzip-ed file
"""

import os

from vcfpy import reader


def test_read_text():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/from_vcf43.vcf')
    r = reader.VCFReader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    EXPECTED = "VCFHeaderLine('fileformat', 'VCFv4.3')"
    assert str(r.header.lines[0]) == EXPECTED
    EXPECTED = (
        "VCFFormatHeaderLine('FORMAT', '<ID=HQ,Number=2,Type=Integer,"
        "Description=\"Haplotype Quality\">', OrderedDict([('ID', 'HQ'), "
        "('Number', 2), ('Type', 'Integer'), ('Description', "
        "'Haplotype Quality')]))")
    assert str(r.header.lines[-1]) == EXPECTED
    assert r.samples
    assert r.samples.names == ['NA00001', 'NA00002', 'NA00003']


def test_read_bgzip():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/from_vcf43.vcf.gz')
    r = reader.VCFReader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    EXPECTED = "VCFHeaderLine('fileformat', 'VCFv4.3')"
    assert str(r.header.lines[0]) == EXPECTED
    EXPECTED = (
        "VCFFormatHeaderLine('FORMAT', '<ID=HQ,Number=2,Type=Integer,"
        "Description=\"Haplotype Quality\">', OrderedDict([('ID', 'HQ'), "
        "('Number', 2), ('Type', 'Integer'), ('Description', "
        "'Haplotype Quality')]))")
    assert str(r.header.lines[-1]) == EXPECTED
    assert r.samples
    assert r.samples.names == ['NA00001', 'NA00002', 'NA00003']
