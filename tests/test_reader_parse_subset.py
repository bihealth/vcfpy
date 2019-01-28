# -*- coding: utf-8 -*-
"""Tests for reading with parsing only a subset of samples"""

from vcfpy import Reader, Writer, Call, UnparsedCall


def test_reading_parse_subset(tmpdir, multisample_vcf_file):
    # Perform record-wise copying, saving results in records
    records = []
    out_path = str(tmpdir.mkdir("output").join("output.vcf"))
    with Reader.from_path(multisample_vcf_file, parsed_samples=["NA00001"]) as reader:
        with Writer.from_path(out_path, reader.header) as writer:
            for record in reader:
                records.append(record)
                writer.write_record(record)
    # Check resulting records, checking the first and last records is enough
    assert len(records) == 5
    assert set(records[0].call_for_sample.keys()) == {"NA00001", "NA00002", "NA00003"}
    assert set(records[-1].call_for_sample.keys()) == {"NA00001", "NA00002", "NA00003"}
    assert isinstance(records[0].call_for_sample["NA00001"], Call)
    assert isinstance(records[0].call_for_sample["NA00002"], UnparsedCall)
    assert isinstance(records[0].call_for_sample["NA00003"], UnparsedCall)
    assert isinstance(records[-1].call_for_sample["NA00001"], Call)
    assert isinstance(records[-1].call_for_sample["NA00002"], UnparsedCall)
    assert isinstance(records[-1].call_for_sample["NA00003"], UnparsedCall)
    # Check resulting file
    with open(multisample_vcf_file, "rt") as inf, open(out_path, "rt") as outf:
        assert inf.read() == outf.read()
