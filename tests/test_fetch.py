# -*- coding: utf-8 -*-
"""Test fetching within tabix files"""

import os

from vcfpy import reader

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


# Test fetch with chrom/begin/end ---------------------------------------------


def test_fetch_no_records_values():
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    records = list(r.fetch("20", 1110698, 1230236))

    assert len(records) == 0


def test_fetch_one_record_values():
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    records = list(r.fetch("20", 1110695, 1230236))

    assert len(records) == 1

    assert records[0].CHROM == "20"
    assert records[0].POS == 1110696


def test_fetch_two_records_values():
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    records = list(r.fetch("20", 1110697, 1234568))

    assert len(records) == 2

    assert records[0].CHROM == "20"
    assert records[0].POS == 1230237

    assert records[1].CHROM == "20"
    assert records[1].POS == 1234567


# Test fetch with samtools region ---------------------------------------------


def test_fetch_no_records_region():
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    records = list(r.fetch("20:1,110,699-1,230,236"))

    assert len(records) == 0


def test_fetch_one_record_region():
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    records = list(r.fetch("20:1,110,696-1,230,236"))

    assert len(records) == 1

    assert records[0].CHROM == "20"
    assert records[0].POS == 1110696


def test_fetch_two_records_region():
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    records = list(r.fetch("20:1,110,698-1,234,568"))

    assert len(records) == 2

    assert records[0].CHROM == "20"
    assert records[0].POS == 1230237

    assert records[1].CHROM == "20"
    assert records[1].POS == 1234567


def test_fetch_multiple_calls_same_reader():
    """Reuse the open index while restarting queries, including unfinished ones."""
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    with reader.Reader.from_path(path) as r:
        records1 = list(r.fetch("20", 1110695, 1230236))
        assert [(record.CHROM, record.POS) for record in records1] == [("20", 1110696)]
        tabix_file = r.tabix_file

        records2 = list(r.fetch("20", 1234566, 1234568))
        assert [(record.CHROM, record.POS) for record in records2] == [("20", 1234567)]
        assert r.tabix_file is tabix_file

        # Rewind using the region-string API, then abandon a partially read query.
        records3 = list(r.fetch("20:1,230,236-1,230,238"))
        assert [(record.CHROM, record.POS) for record in records3] == [("20", 1230237)]
        assert next(r.fetch("20")).POS == 1110696
        assert [(record.CHROM, record.POS) for record in r.fetch("1")] == [("1", 14370)]
        assert list(r.fetch("20:1-10")) == []
        assert [(record.CHROM, record.POS) for record in r.fetch("2")] == [("2", 17330)]
        assert r.tabix_file is tabix_file
        assert not tabix_file.closed

    assert tabix_file.closed
    assert r.stream.closed
