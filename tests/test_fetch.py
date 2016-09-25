# -*- coding: utf-8 -*-
"""Test fetching within tabix files
"""

import os

from vcfpy import reader

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


def test_fetch_no_records():
    path = os.path.join(os.path.dirname(__file__), 'vcfs',
                        'multi_contig.vcf.gz')
    r = reader.Reader.from_path(path)

    records = [vcf_rec for vcf_rec in r.fetch('20', 1110698, 1230236)]

    assert len(records) == 0


def test_fetch_one_record():
    path = os.path.join(os.path.dirname(__file__), 'vcfs',
                        'multi_contig.vcf.gz')
    r = reader.Reader.from_path(path)

    records = [vcf_rec for vcf_rec in r.fetch('20', 1110695, 1230236)]

    assert len(records) == 1

    assert records[0].CHROM == '20'
    assert records[0].POS == 1110696


def test_fetch_two_records():
    path = os.path.join(os.path.dirname(__file__), 'vcfs',
                        'multi_contig.vcf.gz')
    r = reader.Reader.from_path(path)

    records = [vcf_rec for vcf_rec in r.fetch('20', 1110697, 1234568)]

    assert len(records) == 2

    assert records[0].CHROM == '20'
    assert records[0].POS == 1230237

    assert records[1].CHROM == '20'
    assert records[1].POS == 1234567
