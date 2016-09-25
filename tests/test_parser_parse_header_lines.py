# -*- coding: utf-8 -*-
"""Test parsing of VCF header lines from strings
"""

import pytest

from vcfpy import header
from vcfpy import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


@pytest.fixture
def warning_helper():
    return parser.WarningHelper()


# parser.StupidHeaderLineParser.parse_key_value() --------------------------


def test_stupid_vcf_header_line_parser_file_format(warning_helper):
    p = parser.StupidHeaderLineParser(warning_helper)
    INPUT = ('fileFormat', 'VCFv4.2')
    EXPECTED = "HeaderLine('fileFormat', 'VCFv4.2')"
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


# parser.MappingHeaderLineParser.parse_key_value() -------------------------


def test_mapping_vcf_header_line_parser_parse_key_value_filter(
        warning_helper):
    p = parser.MappingHeaderLineParser(
        warning_helper, header.FilterHeaderLine)
    INPUT = ('FILTER', '<ID=q10,Description="Quality below 10">')
    EXPECTED = ('FilterHeaderLine(\'FILTER\', \'<ID=q10,Description="'
                'Quality below 10">\', '
                "OrderedDict([('ID', 'q10'), ('Description', "
                "'Quality below 10')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_format(
        warning_helper):
    p = parser.MappingHeaderLineParser(
        warning_helper, header.FormatHeaderLine)
    INPUT = ('FORMAT', '<ID=GT,Number=1,Type=String,Description="Genotype">')
    EXPECTED = ('FormatHeaderLine(\'FORMAT\', \'<ID=GT,Number=1,'
                'Type=String,Description="Genotype">\', '
                "OrderedDict([('ID', 'GT'), ('Number', 1), "
                "('Type', 'String'), ('Description', 'Genotype')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_info(
        warning_helper):
    p = parser.MappingHeaderLineParser(
        warning_helper, header.InfoHeaderLine)
    INPUT = ('INFO',
             '<ID=NS,Number=1,Type=Integer,Description='
             '"Number of Samples With Data">')
    EXPECTED = ('InfoHeaderLine(\'INFO\', \'<ID=NS,Number=1,Type=Integer,'
                'Description="Number of Samples With Data">\', '
                "OrderedDict([('ID', 'NS'), ('Number', 1), "
                "('Type', 'Integer'), ('Description', "
                "'Number of Samples With Data')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_contig(
        warning_helper):
    p = parser.MappingHeaderLineParser(
        warning_helper, header.ContigHeaderLine)
    INPUT = ('contig',
             '<ID=20,length=62435964,assembly=B36,'
             'md5=f126cdf8a6e0c7f379d618ff66beb2da,'
             'species="Homo sapiens",taxonomy=x>')
    EXPECTED = ('ContigHeaderLine(\'contig\', \'<ID=20,length=62435964,'
                'assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species='
                '"Homo sapiens",taxonomy=x>\', '
                "OrderedDict([('ID', '20'), ('length', '62435964'), "
                "('assembly', 'B36'), "
                "('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'), "
                "('species', 'Homo sapiens'), ('taxonomy', 'x')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


# parser.HeaderParser.parse_line() -----------------------------------------


def test_vcf_header_parser_file_format(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = '##fileFormat=VCFv4.2\n'
    EXPECTED = "HeaderLine('fileFormat', 'VCFv4.2')"
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_vcf_header_parser_parse_line_filter(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = '##FILTER=<ID=q10,Description="Quality below 10">\n'
    EXPECTED = ('FilterHeaderLine(\'FILTER\', \'<ID=q10,Description="'
                'Quality below 10">\', '
                "OrderedDict([('ID', 'q10'), ('Description', "
                "'Quality below 10')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_format(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    EXPECTED = ('FormatHeaderLine(\'FORMAT\', \'<ID=GT,Number=1,'
                'Type=String,Description="Genotype">\', '
                "OrderedDict([('ID', 'GT'), ('Number', 1), "
                "('Type', 'String'), ('Description', 'Genotype')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_info(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = ('##INFO='
             '<ID=NS,Number=1,Type=Integer,Description='
             '"Number of Samples With Data">\n')
    EXPECTED = ('InfoHeaderLine(\'INFO\', \'<ID=NS,Number=1,Type=Integer,'
                'Description="Number of Samples With Data">\', '
                "OrderedDict([('ID', 'NS'), ('Number', 1), "
                "('Type', 'Integer'), ('Description', "
                "'Number of Samples With Data')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_contig(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = ('##contig='
             '<ID=20,length=62435964,assembly=B36,'
             'md5=f126cdf8a6e0c7f379d618ff66beb2da,'
             'species="Homo sapiens",taxonomy=x>\n')
    EXPECTED = ('ContigHeaderLine(\'contig\', \'<ID=20,length=62435964,'
                'assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species='
                '"Homo sapiens",taxonomy=x>\', '
                "OrderedDict([('ID', '20'), ('length', '62435964'), "
                "('assembly', 'B36'), "
                "('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'), "
                "('species', 'Homo sapiens'), ('taxonomy', 'x')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_alt_allele(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = ('##ALT='
             '<ID=R,Description="IUPAC code R = A/G">\n')
    EXPECTED = ("AltAlleleHeaderLine('ALT', "
                "'<ID=R,Description=\"IUPAC code R = A/G\">', "
                "OrderedDict([('ID', 'R'), "
                "('Description', 'IUPAC code R = A/G')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_meta(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = ('##META='
             '<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>\n')
    EXPECTED = (
        "MetaHeaderLine('META', '<ID=Assay,Type=String,Number=.,"
        "Values=[WholeGenome, Exome]>', OrderedDict([('ID', 'Assay'), "
        "('Type', 'String'), ('Number', '.'), ('Values', ['WholeGenome', "
        "'Exome'])]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_pedigree(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = ('##PEDIGREE='
             '<ID=TumourSample,Original=GermlineID>\n')
    EXPECTED = ("PedigreeHeaderLine('PEDIGREE', "
                "'<ID=TumourSample,Original=GermlineID>',"
                " OrderedDict([('ID', 'TumourSample'), "
                "('Original', 'GermlineID')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_sample(warning_helper):
    p = parser.HeaderParser(parser.build_header_parsers(warning_helper))
    INPUT = ('##SAMPLE='
             '<ID=Sample1,Assay=WholeGenome,Ethnicity=AFR,Disease=None,'
             'Description="Patient germline genome from unaffected",'
             'DOI=url>\n')
    EXPECTED = (
        "SampleHeaderLine('SAMPLE', '<ID=Sample1,Assay=WholeGenome,"
        'Ethnicity=AFR,Disease=None,Description="Patient germline genome from '
        "unaffected\",DOI=url>', OrderedDict([('ID', 'Sample1'), ('Assay', "
        "'WholeGenome'), ('Ethnicity', 'AFR'), ('Disease', 'None'), "
        "('Description', 'Patient germline genome from unaffected'), "
        "('DOI', 'url')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED
