# -*- coding: utf-8 -*-
"""Test for the helper routines in the vcfpy.parser module"""

from vcfpy import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# parser.split_quoted_string() ------------------------------------------------


def test_split_quoted_string_one_noquote():
    INPUT = 'foo=bar'
    EXPECTED = ['foo=bar']
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_two_noquote():
    INPUT = 'foo=bar,bar=baz'
    EXPECTED = ['foo=bar', 'bar=baz']
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_one_quote():
    INPUT = 'foo="bar,,"asdf'
    EXPECTED = ['foo="bar,,"asdf']
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_two_quote():
    INPUT = 'foo="bar,,"asdf,"bar"=baz'
    EXPECTED = ['foo="bar,,"asdf', '"bar"=baz']
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_one_escape():
    INPUT = 'foo="bar\\"asdf\\"asdf"'
    EXPECTED = ['foo="bar\\"asdf\\"asdf"']
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_two_escape():
    INPUT = 'foo="bar\\",,asdf","bar"=baz'
    EXPECTED = ['foo="bar\\",,asdf"', '"bar"=baz']
    assert EXPECTED == parser.split_quoted_string(INPUT)


# parser.VCFheaderLineParser.parse_mapping() ----------------------------------


def test_vcf_header_line_parser_parse_mapping_simple():
    INPUT = r'<key=value,key2=value2>'
    EXPECTED = (('key', 'value'),
                ('key2', 'value2'))
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_vcf_header_line_parser_parse_mapping_flag():
    INPUT = r'<key=value,key2=value,yay>'
    EXPECTED = (('key', 'value'),
                ('key2', 'value'),
                ('yay', True))
    p = parser.MappingHeaderLineParser(None)
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_vcf_header_line_parser_parse_mapping_quoted():
    INPUT = r'<key=value,key2="value,value">'
    EXPECTED = (('key', 'value'),
                ('key2', 'value,value'))
    p = parser.MappingHeaderLineParser(None)
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_vcf_header_line_parser_parse_mapping_escaped():
    INPUT = r'<key=value,key2="value,value=\"asdf">'
    EXPECTED = (('key', 'value'),
                ('key2', 'value,value="asdf'))
    p = parser.MappingHeaderLineParser(None)
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())
