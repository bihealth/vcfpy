# -*- coding: utf-8 -*-
"""Test parsing of VCF header lines from strings
"""

import sys

import pytest

from vcfpy import header
from vcfpy import parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


# parser.StupidHeaderLineParser.parse_key_value() --------------------------


def test_stupid_vcf_header_line_parser_file_format():
    p = parser.StupidHeaderLineParser()
    INPUT = ("fileFormat", "VCFv4.2")
    EXPECTED = "HeaderLine('fileFormat', 'VCFv4.2')"
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


# parser.MappingHeaderLineParser.parse_key_value() -------------------------


def test_mapping_vcf_header_line_parser_parse_key_value_filter():
    p = parser.MappingHeaderLineParser(header.FilterHeaderLine)
    INPUT = ("FILTER", '<ID=q10,Description="Quality below 10">')
    if sys.version_info < (3, 6):
        EXPECTED = (
            "FilterHeaderLine('FILTER', '<ID=q10,Description=\"Quality below 10\">', "
            "OrderedDict([('ID', 'q10'), ('Description', 'Quality below 10')]))"
        )
    else:
        EXPECTED = (
            "FilterHeaderLine('FILTER', '<ID=q10,Description=\"Quality below 10\">', "
            "{'ID': 'q10', 'Description': 'Quality below 10'})"
        )
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_format():
    p = parser.MappingHeaderLineParser(header.FormatHeaderLine)
    INPUT = ("FORMAT", '<ID=GT,Number=1,Type=String,Description="Genotype">')
    if sys.version_info < (3, 6):
        EXPECTED = (
            "FormatHeaderLine('FORMAT', '<ID=GT,Number=1,Type=String,Description=\"Genotype\">', "
            "OrderedDict([('ID', 'GT'), ('Number', 1), ('Type', 'String'), ('Description', 'Genotype')]))"
        )
    else:
        EXPECTED = (
            "FormatHeaderLine('FORMAT', '<ID=GT,Number=1,Type=String,Description=\"Genotype\">', "
            "{'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'})"
        )
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_info():
    p = parser.MappingHeaderLineParser(header.InfoHeaderLine)
    INPUT = ("INFO", "<ID=NS,Number=1,Type=Integer,Description=" '"Number of Samples With Data">')
    if sys.version_info < (3, 6):
        EXPECTED = (
            "InfoHeaderLine('INFO', '<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">', "
            "OrderedDict([('ID', 'NS'), ('Number', 1), ('Type', 'Integer'), "
            "('Description', 'Number of Samples With Data')]))"
        )
    else:
        EXPECTED = (
            "InfoHeaderLine('INFO', '<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">', "
            "{'ID': 'NS', 'Number': 1, 'Type': 'Integer', 'Description': 'Number of Samples With Data'})"
        )
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_contig():
    p = parser.MappingHeaderLineParser(header.ContigHeaderLine)
    INPUT = (
        "contig",
        '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>',
    )
    if sys.version_info < (3, 6):
        EXPECTED = (
            "ContigHeaderLine('contig', '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,"
            "species=\"Homo sapiens\",taxonomy=x>', OrderedDict([('ID', '20'), ('length', '62435964'), "
            "('assembly', 'B36'), ('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'), ('species', 'Homo sapiens'), "
            "('taxonomy', 'x')]))"
        )
    else:
        EXPECTED = (
            "ContigHeaderLine('contig', '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,"
            "species=\"Homo sapiens\",taxonomy=x>', {'ID': '20', 'length': '62435964', 'assembly': 'B36', "
            "'md5': 'f126cdf8a6e0c7f379d618ff66beb2da', 'species': 'Homo sapiens', 'taxonomy': 'x'})"
        )
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


# parser.HeaderParser.parse_line() -----------------------------------------


def test_vcf_header_parser_file_format():
    p = parser.HeaderParser()
    INPUT = "##fileFormat=VCFv4.2\n"
    EXPECTED = "HeaderLine('fileFormat', 'VCFv4.2')"
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_vcf_header_parser_parse_line_filter():
    p = parser.HeaderParser()
    INPUT = '##FILTER=<ID=q10,Description="Quality below 10">\n'
    if sys.version_info < (3, 6):
        EXPECTED = (
            "FilterHeaderLine('FILTER', '<ID=q10,Description=\"Quality below 10\">', OrderedDict([('ID', 'q10'), "
            "('Description', 'Quality below 10')]))"
        )
    else:
        EXPECTED = (
            "FilterHeaderLine('FILTER', '<ID=q10,Description=\"Quality below 10\">', "
            "{'ID': 'q10', 'Description': 'Quality below 10'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_format():
    p = parser.HeaderParser()
    INPUT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    if sys.version_info < (3, 6):
        EXPECTED = (
            "FormatHeaderLine('FORMAT', '<ID=GT,Number=1,Type=String,Description=\"Genotype\">', "
            "OrderedDict([('ID', 'GT'), ('Number', 1), ('Type', 'String'), ('Description', 'Genotype')]))"
        )
    else:
        EXPECTED = (
            "FormatHeaderLine('FORMAT', '<ID=GT,Number=1,Type=String,Description=\"Genotype\">', "
            "{'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_info():
    p = parser.HeaderParser()
    INPUT = "##INFO=" "<ID=NS,Number=1,Type=Integer,Description=" '"Number of Samples With Data">\n'
    if sys.version_info < (3, 6):
        EXPECTED = (
            "InfoHeaderLine('INFO', '<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">', "
            "OrderedDict([('ID', 'NS'), ('Number', 1), ('Type', 'Integer'), ('Description', "
            "'Number of Samples With Data')]))"
        )
    else:
        EXPECTED = (
            "InfoHeaderLine('INFO', '<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">', "
            "{'ID': 'NS', 'Number': 1, 'Type': 'Integer', 'Description': 'Number of Samples With Data'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_contig():
    p = parser.HeaderParser()
    INPUT = (
        "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,"
        'species="Homo sapiens",taxonomy=x>\n'
    )
    if sys.version_info < (3, 6):
        EXPECTED = (
            "ContigHeaderLine('contig', '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,"
            "species=\"Homo sapiens\",taxonomy=x>', OrderedDict([('ID', '20'), ('length', '62435964'), "
            "('assembly', 'B36'), ('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'), ('species', 'Homo sapiens'), "
            "('taxonomy', 'x')]))"
        )
    else:
        EXPECTED = (
            "ContigHeaderLine('contig', '<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,"
            "species=\"Homo sapiens\",taxonomy=x>', {'ID': '20', 'length': '62435964', 'assembly': 'B36', "
            "'md5': 'f126cdf8a6e0c7f379d618ff66beb2da', 'species': 'Homo sapiens', 'taxonomy': 'x'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_alt_allele():
    p = parser.HeaderParser()
    INPUT = "##ALT=" '<ID=R,Description="IUPAC code R = A/G">\n'
    if sys.version_info < (3, 6):
        EXPECTED = (
            "AltAlleleHeaderLine('ALT', "
            "'<ID=R,Description=\"IUPAC code R = A/G\">', "
            "OrderedDict([('ID', 'R'), "
            "('Description', 'IUPAC code R = A/G')]))"
        )
    else:
        EXPECTED = (
            "AltAlleleHeaderLine('ALT', "
            "'<ID=R,Description=\"IUPAC code R = A/G\">', "
            "{'ID': 'R', 'Description': 'IUPAC code R = A/G'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_meta():
    p = parser.HeaderParser()
    INPUT = "##META=" "<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>\n"
    if sys.version_info < (3, 6):
        EXPECTED = (
            "MetaHeaderLine('META', '<ID=Assay,Type=String,Number=.,"
            "Values=[WholeGenome, Exome]>', OrderedDict([('ID', 'Assay'), "
            "('Type', 'String'), ('Number', '.'), ('Values', ['WholeGenome', "
            "'Exome'])]))"
        )
    else:
        EXPECTED = (
            "MetaHeaderLine('META', '<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>', "
            "{'ID': 'Assay', 'Type': 'String', 'Number': '.', 'Values': ['WholeGenome', 'Exome']})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_pedigree():
    p = parser.HeaderParser()
    INPUT = "##PEDIGREE=" "<ID=TumourSample,Original=GermlineID>\n"
    if sys.version_info < (3, 6):
        EXPECTED = (
            "PedigreeHeaderLine('PEDIGREE', "
            "'<ID=TumourSample,Original=GermlineID>',"
            " OrderedDict([('ID', 'TumourSample'), "
            "('Original', 'GermlineID')]))"
        )
    else:
        EXPECTED = (
            "PedigreeHeaderLine('PEDIGREE', '<ID=TumourSample,Original=GermlineID>', "
            "{'ID': 'TumourSample', 'Original': 'GermlineID'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_sample():
    p = parser.HeaderParser()
    INPUT = (
        "##SAMPLE="
        "<ID=Sample1,Assay=WholeGenome,Ethnicity=AFR,Disease=None,"
        'Description="Patient germline genome from unaffected",'
        "DOI=url>\n"
    )
    if sys.version_info < (3, 6):
        EXPECTED = (
            "SampleHeaderLine('SAMPLE', '<ID=Sample1,Assay=WholeGenome,"
            'Ethnicity=AFR,Disease=None,Description="Patient germline genome from '
            "unaffected\",DOI=url>', OrderedDict([('ID', 'Sample1'), ('Assay', "
            "'WholeGenome'), ('Ethnicity', 'AFR'), ('Disease', 'None'), "
            "('Description', 'Patient germline genome from unaffected'), "
            "('DOI', 'url')]))"
        )
    else:
        EXPECTED = (
            "SampleHeaderLine('SAMPLE', '<ID=Sample1,Assay=WholeGenome,"
            'Ethnicity=AFR,Disease=None,Description="Patient germline genome from '
            "unaffected\",DOI=url>', {'ID': 'Sample1', 'Assay': 'WholeGenome', 'Ethnicity': 'AFR', 'Disease': 'None', "
            "'Description': 'Patient germline genome from unaffected', 'DOI': 'url'})"
        )
    assert str(p.parse_line(INPUT)) == EXPECTED
