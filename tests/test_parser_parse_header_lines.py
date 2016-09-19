# -*- coding: utf-8 -*-
"""Test parsing of VCF header lines from strings
"""

from vcfpy import header
from vcfpy import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# parser.StupidVCFHeaderLineParser.parse_key_value() --------------------------


def test_stupid_vcf_header_line_parser_file_format():
    p = parser.StupidVCFHeaderLineParser()
    INPUT = ('fileFormat', 'VCFv4.2')
    EXPECTED = "VCFHeaderLine('fileFormat', 'VCFv4.2')"
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


# parser.MappingVCFHeaderLineParser.parse_key_value() -------------------------


def test_mapping_vcf_header_line_parser_parse_key_value_filter():
    p = parser.MappingVCFHeaderLineParser(header.VCFFilterHeaderLine)
    INPUT = ('FILTER', '<ID=q10,Description="Quality below 10">')
    EXPECTED = ('VCFFilterHeaderLine(\'FILTER\', \'<ID=q10,Description="'
                'Quality below 10">\', '
                "OrderedDict([('ID', 'q10'), ('Description', "
                "'Quality below 10')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_format():
    p = parser.MappingVCFHeaderLineParser(header.VCFFormatHeaderLine)
    INPUT = ('FORMAT', '<ID=GT,Number=1,Type=String,Description="Genotype">')
    EXPECTED = ('VCFFormatHeaderLine(\'FORMAT\', \'<ID=GT,Number=1,'
                'Type=String,Description="Genotype">\', '
                "OrderedDict([('ID', 'GT'), ('Number', 1), "
                "('Type', 'String'), ('Description', 'Genotype')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_info():
    p = parser.MappingVCFHeaderLineParser(header.VCFInfoHeaderLine)
    INPUT = ('INFO',
             '<ID=NS,Number=1,Type=Integer,Description='
             '"Number of Samples With Data">')
    EXPECTED = ('VCFInfoHeaderLine(\'INFO\', \'<ID=NS,Number=1,Type=Integer,'
                'Description="Number of Samples With Data">\', '
                "OrderedDict([('ID', 'NS'), ('Number', 1), "
                "('Type', 'Integer'), ('Description', "
                "'Number of Samples With Data')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


def test_mapping_vcf_header_line_parser_parse_key_value_contig():
    p = parser.MappingVCFHeaderLineParser(header.VCFContigHeaderLine)
    INPUT = ('contig',
             '<ID=20,length=62435964,assembly=B36,'
             'md5=f126cdf8a6e0c7f379d618ff66beb2da,'
             'species="Homo sapiens",taxonomy=x>')
    EXPECTED = ('VCFContigHeaderLine(\'contig\', \'<ID=20,length=62435964,'
                'assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species='
                '"Homo sapiens",taxonomy=x>\', '
                "OrderedDict([('ID', '20'), ('length', '62435964'), "
                "('assembly', 'B36'), "
                "('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'), "
                "('species', 'Homo sapiens'), ('taxonomy', 'x')]))")
    assert str(p.parse_key_value(*INPUT)) == EXPECTED


# parser.VCFHeaderParser.parse_line() -----------------------------------------

def test_vcf_header_parser_file_format():
    p = parser.VCFHeaderParser(parser.HEADER_PARSERS)
    INPUT = '##fileFormat=VCFv4.2\n'
    EXPECTED = "VCFHeaderLine('fileFormat', 'VCFv4.2')"
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_vcf_header_parser_parse_line_filter():
    p = parser.VCFHeaderParser(parser.HEADER_PARSERS)
    INPUT = '##FILTER=<ID=q10,Description="Quality below 10">\n'
    EXPECTED = ('VCFFilterHeaderLine(\'FILTER\', \'<ID=q10,Description="'
                'Quality below 10">\', '
                "OrderedDict([('ID', 'q10'), ('Description', "
                "'Quality below 10')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_format():
    p = parser.VCFHeaderParser(parser.HEADER_PARSERS)
    INPUT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    EXPECTED = ('VCFFormatHeaderLine(\'FORMAT\', \'<ID=GT,Number=1,'
                'Type=String,Description="Genotype">\', '
                "OrderedDict([('ID', 'GT'), ('Number', 1), "
                "('Type', 'String'), ('Description', 'Genotype')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_info():
    p = parser.VCFHeaderParser(parser.HEADER_PARSERS)
    INPUT = ('##INFO='
             '<ID=NS,Number=1,Type=Integer,Description='
             '"Number of Samples With Data">\n')
    EXPECTED = ('VCFInfoHeaderLine(\'INFO\', \'<ID=NS,Number=1,Type=Integer,'
                'Description="Number of Samples With Data">\', '
                "OrderedDict([('ID', 'NS'), ('Number', 1), "
                "('Type', 'Integer'), ('Description', "
                "'Number of Samples With Data')]))")
    assert str(p.parse_line(INPUT)) == EXPECTED


def test_mapping_vcf_header_parser_parse_line_contig():
    p = parser.VCFHeaderParser(parser.HEADER_PARSERS)
    INPUT = ('##contig='
             '<ID=20,length=62435964,assembly=B36,'
             'md5=f126cdf8a6e0c7f379d618ff66beb2da,'
             'species="Homo sapiens",taxonomy=x>\n')
    EXPECTED = ('VCFContigHeaderLine(\'contig\', \'<ID=20,length=62435964,'
                'assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species='
                '"Homo sapiens",taxonomy=x>\', '
                "OrderedDict([('ID', '20'), ('length', '62435964'), "
                "('assembly', 'B36'), "
                "('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'), "
                "('species', 'Homo sapiens'), ('taxonomy', 'x')]))")
    record = p.parse_line(INPUT)
    assert str(p.parse_line(INPUT)) == EXPECTED
