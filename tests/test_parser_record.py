# -*- coding: utf-8 -*-
"""Test parsing of full VCF record lines"""

import io
import warnings
from io import StringIO

import pytest

import vcfpy
from vcfpy import exceptions, parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


SMALL_HEADER = """
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


MEDIUM_HEADER = """
##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##INFO=<ID=ANNO,Number=.,Type=String,Description="Additional annotation">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FILTER=<ID=FOO,Description="A test filter">
##FILTER=<ID=BAZ,Description="A second test filter">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Call-wise filters">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


def vcf_parser(lines):
    return parser.Parser(io.StringIO(MEDIUM_HEADER + lines), "<builtin>")


def test_parse_minimal_record():
    # Setup parser with stock header and lines to parse
    LINES = "20\t1\t.\tC\tG\t.\t.\t.\tGT\t0/1\t0/2\t.\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 1, [], 'C', [Substitution(type_='SNV', value='G')], None, [], {}, ['GT'], "
        "[Call('NA00001', {'GT': '0/1'}), Call('NA00002', {'GT': '0/2'}), Call('NA00003', {'GT': None})])"
    )
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED


def test_parse_record_with_info():
    # Setup parser with stock header and lines to parse
    LINES = "20\t1\t.\tC\tG\t.\t.\tAA=G\tGT\t0/1\t0/1\t.\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 1, [], 'C', [Substitution(type_='SNV', value='G')], None, [], {'AA': 'G'}, "
        "['GT'], [Call('NA00001', {'GT': '0/1'}), Call('NA00002', {'GT': '0/1'}), Call('NA00003', {'GT': None})])"
    )
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED


def test_parse_record_with_escaping():
    # Setup parser with stock header and lines to parse
    LINES = "20\t100\t.\tC\tG\t.\t.\tANNO=Here%2Care%25some chars,%2525\tGT:FT\t0/1:FOO\t0/0:.\t1/1:.\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 100, [], 'C', [Substitution(type_='SNV', value='G')], None, [], "
        "{'ANNO': ['Here,are%some chars', '%25']}, ['GT', 'FT'], "
        "[Call('NA00001', {'GT': '0/1', 'FT': ['FOO']}), Call('NA00002', {'GT': '0/0', 'FT': []}), "
        "Call('NA00003', {'GT': '1/1', 'FT': []})])"
    )
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED


def test_parse_record_with_filter_warning():
    # Setup parser with stock header and lines to parse
    LINES = "20\t1\t.\tC\tG\t.\tBAZ\t.\tGT:FT\t0/1:.\t0/2:BAR\t.:.\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 1, [], 'C', [Substitution(type_='SNV', value='G')], None, ['BAZ'], {}, "
        "['GT', 'FT'], [Call('NA00001', {'GT': '0/1', 'FT': []}),"
        " Call('NA00002', {'GT': '0/2', 'FT': ['BAR']}), Call('NA00003', {'GT': None, 'FT': []})])"
    )
    with pytest.warns((exceptions.UnknownFilter,)):
        RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED


def test_parse_record_with_filter_no_warning(recwarn: pytest.WarningsRecorder):
    # Setup parser with stock header and lines to parse
    LINES = "20\t1\t.\tC\tG\t.\tFOO\t.\tGT:FT\t0/1:.\t0/2:FOO;BAZ\t.:.\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 1, [], 'C', [Substitution(type_='SNV', value='G')], None, ['FOO'], {}, "
        "['GT', 'FT'], [Call('NA00001', {'GT': '0/1', 'FT': []}),"
        " Call('NA00002', {'GT': '0/2', 'FT': ['FOO', 'BAZ']}), Call('NA00003', {'GT': None, 'FT': []})])"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED


def test_missing_pass(recwarn: pytest.WarningsRecorder):
    """Test parsing VCF file with missing ``PASS`` definition in header."""
    # Setup parser with stock header and lines to parse
    LINES = "20\t1\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t0/2\t.\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 1, [], 'C', [Substitution(type_='SNV', value='G')], None, ['PASS'], {}, ['GT'], "
        "[Call('NA00001', {'GT': '0/1'}), Call('NA00002', {'GT': '0/2'}), Call('NA00003', {'GT': None})])"
    )
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert list(recwarn) == []


def test_parse_line_invalid_number_of_fields():
    """Test parsing VCF file exception message"""
    HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
    LINES = "20\t1\t.\tC\tG\t.\tPASS\t.\tGT\n"
    EXPECTED = "The line contains an invalid number of fields. Was 9 but expected 8"
    p = parser.Parser(io.StringIO(HEADER + LINES), "<builtin>")
    p.parse_header()
    with pytest.raises(exceptions.InvalidRecordException) as record_error:
        p.parse_next_record()

    assert EXPECTED in str(record_error.value)


def test_parse_record_with_gt_data():
    LINES = "20\t1\t.\tC\tG\t.\t.\t.\tGT\t0|1\t1/1\t.\n"
    p = parser.Parser(io.StringIO(MEDIUM_HEADER + LINES))
    p.parse_header()
    record = p.parse_next_record()
    assert record.calls[0].data["GT"] == "0|1"
    assert record.calls[1].data["GT"] == "1/1"
    assert record.calls[2].data["GT"] is None


def test_reader_parse_warnings():
    """Test reader parsing with warnings for missing header lines"""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	30	PASS	DP=10
"""

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        stream = StringIO(vcf_content)
        reader = vcfpy.Reader.from_stream(stream)
        # Just ensure it parses without error
        records = list(reader)
        assert len(records) == 1


def test_parser_quoted_string_splitter_edge_cases():
    """Test QuotedStringSplitter with edge cases"""
    splitter = parser.QuotedStringSplitter()

    # Test with escaped quotes
    result = splitter.run('key="val\\"with\\"quotes"')
    assert len(result) == 1

    # Test with nested brackets
    result = splitter.run("key1=[nested[value]],key2=simple")
    assert len(result) == 2


def test_parse_mapping_correct_format():
    """Test parse_mapping with correct angular bracket format"""
    # Test with proper angular brackets
    result = parser.parse_mapping('<ID=test,Description="Complex value with, commas">')
    assert "ID" in result
    assert "Description" in result
    assert result["ID"] == "test"
