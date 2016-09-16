# -*- coding: utf-8 -*-
"""Reading of VCF files from plain and bgzip-ed files
"""

import os

from vcfpy import reader


def test_read_text():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/full_vcf43.vcf')
    r = reader.VCFReader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    assert r.samples
    assert r.samples.names == ['NA00001', 'NA00002', 'NA00003']
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5


def test_read_bgzip():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/full_vcf43.vcf.gz')
    r = reader.VCFReader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    assert r.samples
    assert r.samples.names == ['NA00001', 'NA00002', 'NA00003']
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5
