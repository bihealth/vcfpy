# -*- coding: utf-8 -*-
"""Test for the helper routines in the vcfpy.parser module"""

import pytest

from vcfpy import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


@pytest.fixture
def warning_helper():
    return parser.WarningHelper()


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


def test_split_quoted_string_array_syntax_simple():
    INPUT = 'foo=[1, 2, 3],bar=[baz, gnaa]'
    EXPECTED = ['foo=[1, 2, 3]', 'bar=[baz, gnaa]']
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_array_syntax_recursion():
    """Support for only one level, failure is OK here"""
    INPUT = 'foo=[1, 2, 3],bar=[baz, [], gnaa]'
    EXPECTED = ['foo=[1, 2, 3]', 'bar=[baz, []', ' gnaa]']
    assert EXPECTED == parser.split_quoted_string(INPUT)


# parser.VCFheaderLineParser.parse_mapping() ----------------------------------


def test_vcf_header_line_parser_parse_mapping_simple(warning_helper):
    INPUT = r'<key=value,key2=value2>'
    EXPECTED = (('key', 'value'),
                ('key2', 'value2'))
    assert EXPECTED == tuple(parser.parse_mapping(
        INPUT, warning_helper).items())


def test_vcf_header_line_parser_parse_mapping_flag(warning_helper):
    INPUT = r'<key=value,key2=value,yay>'
    EXPECTED = (('key', 'value'),
                ('key2', 'value'),
                ('yay', True))
    parser.MappingHeaderLineParser(warning_helper, None)
    assert EXPECTED == tuple(parser.parse_mapping(
        INPUT, warning_helper).items())


def test_vcf_header_line_parser_parse_mapping_quoted(warning_helper):
    INPUT = r'<key=value,key2="value,value">'
    EXPECTED = (('key', 'value'),
                ('key2', 'value,value'))
    parser.MappingHeaderLineParser(warning_helper, None)
    assert EXPECTED == tuple(parser.parse_mapping(
        INPUT, warning_helper).items())


def test_vcf_header_line_parser_parse_mapping_escaped(warning_helper):
    INPUT = r'<key=value,key2="value,value=\"asdf">'
    EXPECTED = (('key', 'value'),
                ('key2', 'value,value="asdf'))
    parser.MappingHeaderLineParser(warning_helper, None)
    assert EXPECTED == tuple(parser.parse_mapping(
        INPUT, warning_helper).items())
