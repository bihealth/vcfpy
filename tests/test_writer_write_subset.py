# -*- coding: utf-8 -*-
"""Tests for reading writing in the case of writing only a subset of the
samples
"""


from vcfpy import Reader, Writer, Header, SamplesInfos, Call, UnparsedCall


def test_reading_and_write_subset_of_samples(
        tmpdir, multisample_vcf_file, multisample_vcf_only_NA00002):
    # Perform record-wise copying, saving results in records, not writing
    # out for some samples
    records = []
    out_path = str(tmpdir.mkdir('output').join('output.vcf'))
    with Reader.from_path(multisample_vcf_file) as reader:
        samples = SamplesInfos(['NA00002'])
        header = Header(reader.header.lines, samples)
        with Writer.from_path(out_path, header) as writer:
            for record in reader:
                records.append(record)
                writer.write_record(record)
    # Check resulting file
    with open(out_path, 'rt') as outf:
        assert multisample_vcf_only_NA00002 == outf.read()


def test_reading_and_write_subset_of_samples_parse_different_subset(
        tmpdir, multisample_vcf_file, multisample_vcf_only_NA00002):
    # Perform record-wise copying, saving results in records, not writing
    # out for some samples
    records = []
    out_path = str(tmpdir.mkdir('output').join('output.vcf'))
    with Reader.from_path(
            multisample_vcf_file, parsed_samples=['NA00001']) as reader:
        samples = SamplesInfos(['NA00002'])
        header = Header(reader.header.lines, samples)
        with Writer.from_path(out_path, header) as writer:
            for record in reader:
                records.append(record)
                writer.write_record(record)
    # Check resulting file
    with open(out_path, 'rt') as outf:
        assert multisample_vcf_only_NA00002 == outf.read()
