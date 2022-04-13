# -*- coding: utf-8 -*-
"""Reading of VCF files from plain and bgzip-ed files
"""

import io
import os

from vcfpy import reader

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def test_read_text():
    path = os.path.join(os.path.dirname(__file__), "vcfs/full_vcf43.vcf")
    r = reader.Reader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    assert r.header.samples
    assert r.header.samples.names == ["NA00001", "NA00002", "NA00003"]
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5


def test_read_bgzip():
    path = os.path.join(os.path.dirname(__file__), "vcfs/full_vcf43.vcf.gz")
    r = reader.Reader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    assert r.header.samples
    assert r.header.samples.names == ["NA00001", "NA00002", "NA00003"]
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5


def test_read_text_no_samples():
    path = os.path.join(os.path.dirname(__file__), "vcfs/full_vcf43_no_samples.vcf")
    r = reader.Reader.from_path(path)
    assert r.parser
    assert r.header
    assert len(r.header.lines) == 18
    assert r.header.samples
    assert r.header.samples.names == []
    records = []
    for record in r:
        records.append(record)
    assert len(records) == 5


def test_read_info_flag():
    """Test reading INFO field flag with inconsistent header metadata."""
    # In the INFO field, `MH` is a flag but the header specifies it as a string.
    string_buffer = io.StringIO(
        r"""##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=MH,Number=1,Type=String,Description="Microhomology">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	FOO	BAR
1	4798729	.	CT	C	621	PASS	MH	GT	0/0:25	0/1:64
"""
    )

    # Verify that no TypeError is raised.
    with reader.Reader.from_stream(string_buffer) as r:
        line = next(r)
        assert line.INFO["MH"]
