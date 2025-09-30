# -*- coding: utf-8 -*-
"""Test for the helper routines in the vcfpy.parser module"""

from vcfpy import parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


# parser.split_quoted_string() ------------------------------------------------


def test_split_quoted_string_one_noquote():
    INPUT = "foo=bar"
    EXPECTED = ["foo=bar"]
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_two_noquote():
    INPUT = "foo=bar,bar=baz"
    EXPECTED = ["foo=bar", "bar=baz"]
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
    INPUT = "foo=[1, 2, 3],bar=[baz, gnaa]"
    EXPECTED = ["foo=[1, 2, 3]", "bar=[baz, gnaa]"]
    assert EXPECTED == parser.split_quoted_string(INPUT)


def test_split_quoted_string_array_syntax_recursion():
    """Support for only one level, failure is OK here"""
    INPUT = "foo=[1, 2, 3],bar=[baz, [], gnaa]"
    EXPECTED = ["foo=[1, 2, 3]", "bar=[baz, []", " gnaa]"]
    assert EXPECTED == parser.split_quoted_string(INPUT)


# parser.VCFheaderLineParser.parse_mapping() ----------------------------------


def test_vcf_header_line_parser_parse_mapping_simple():
    INPUT = r"<key=value,key2=value2>"
    EXPECTED = (("key", "value"), ("key2", "value2"))
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_vcf_header_line_parser_parse_mapping_flag():
    INPUT = r"<key=value,key2=value,yay>"
    EXPECTED = (("key", "value"), ("key2", "value"), ("yay", True))
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_vcf_header_line_parser_parse_mapping_quoted():
    INPUT = r'<key=value,key2="value,value">'
    EXPECTED = (("key", "value"), ("key2", "value,value"))
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_vcf_header_line_parser_parse_mapping_escaped():
    INPUT = r'<key=value,key2="value,value=\"asdf">'
    EXPECTED = (("key", "value"), ("key2", 'value,value="asdf'))
    assert EXPECTED == tuple(parser.parse_mapping(INPUT).items())


def test_parse_mapping_invalid_brackets():
    """Test parse_mapping with invalid angular brackets"""
    import pytest

    from vcfpy import exceptions

    # Missing opening bracket
    with pytest.raises(
        exceptions.InvalidHeaderException, match="Header mapping value was not wrapped in angular brackets"
    ):
        parser.parse_mapping("key=value>")

    # Missing closing bracket
    with pytest.raises(
        exceptions.InvalidHeaderException, match="Header mapping value was not wrapped in angular brackets"
    ):
        parser.parse_mapping("<key=value")

    # No brackets at all
    with pytest.raises(
        exceptions.InvalidHeaderException, match="Header mapping value was not wrapped in angular brackets"
    ):
        parser.parse_mapping("key=value")


def test_parse_mapping_with_array_values():
    """Test parse_mapping with array-style values"""
    INPUT = "<ID=DP,Number=A,Type=Float,Values=[1,2,3]>"
    result = parser.parse_mapping(INPUT)
    assert result["Values"] == ["1", "2", "3"]


def test_parse_mapping_with_flags():
    """Test parse_mapping with flag (no value)"""
    INPUT = "<ID=DP,flag>"
    result = parser.parse_mapping(INPUT)
    assert result["ID"] == "DP"
    assert result["flag"] is True


def test_split_mapping_edge_cases():
    """Test split_mapping function with edge cases"""
    # Test basic case
    result = parser.split_mapping("key=value")
    assert result == ("key", "value")

    # Test with equals in value (should only split on first equals)
    result = parser.split_mapping("key=value=with=equals")
    assert result == ("key", "value=with=equals")


def test_split_quoted_string_edge_cases():
    """Test split_quoted_string with additional edge cases"""
    # Test with different delimiters
    result = parser.split_quoted_string("a;b;c", delim=";")
    assert result == ["a", "b", "c"]

    # Test with different quote characters
    result = parser.split_quoted_string("a='b,c',d", quote="'")
    assert result == ["a='b,c'", "d"]

    # Test empty string
    result = parser.split_quoted_string("")
    assert result == [""]

    # Test single value
    result = parser.split_quoted_string("single")
    assert result == ["single"]
