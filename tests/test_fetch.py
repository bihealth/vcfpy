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
    """Test calling fetch multiple times on the same reader to cover tabix_file.close() path"""
    path = os.path.join(os.path.dirname(__file__), "vcfs", "multi_contig.vcf.gz")
    r = reader.Reader.from_path(path)

    # First fetch call - this will open tabix_file
    records1 = list(r.fetch("20", 1110695, 1230236))
    assert len(records1) == 1
    assert records1[0].CHROM == "20"
    assert records1[0].POS == 1110696

    # Second fetch call - this will close the existing tabix_file and open a new one
    # This triggers line 157: self.tabix_file.close()
    records2 = list(r.fetch("20", 1234566, 1234568))
    assert len(records2) == 1
    assert records2[0].CHROM == "20"
    assert records2[0].POS == 1234567

    # Third fetch with region string to ensure it also works
    records3 = list(r.fetch("20:1,230,236-1,230,238"))
    assert len(records3) == 1
    assert records3[0].CHROM == "20"
    assert records3[0].POS == 1230237
