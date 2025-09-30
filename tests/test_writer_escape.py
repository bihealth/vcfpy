# -*- coding: utf-8 -*-
"""Test escaping on writing VCF records"""

from vcfpy import header, writer

# vcfpy.writer.format_atomic() ------------------------------------------------


def test_format_atomic_with_escape_info():
    expected = "%3B%3D%25%2C%0D%0A%09"
    result = writer.format_atomic(";=%,\r\n\t", "INFO")
    assert expected == result


def test_format_atomic_with_escape_format():
    expected = "%3A%3D%25%2C%0D%0A%09"
    result = writer.format_atomic(":=%,\r\n\t", "FORMAT")
    assert expected == result


def test_format_atomic_without_escape_info():
    expected = "This is a legal string:"
    result = writer.format_atomic("This is a legal string:", "INFO")
    assert expected == result


def test_format_atomic_without_escape_format():
    expected = "This is a legal string;"
    result = writer.format_atomic("This is a legal string;", "FORMAT")
    assert expected == result


# vcfpy.writer.format_value() -------------------------------------------------


def test_format_value_with_escape():
    expected = "%3A%3B%3D%25%2C%0D%0A%09,%25"
    result = writer.format_value(header.FieldInfo("String", 2), [":;=%,\r\n\t", "%"], "INFO")
    assert expected == result


def test_format_value_without_escape():
    expected = "This is a legal string,me too"
    result = writer.format_value(header.FieldInfo("String", 2), ["This is a legal string", "me too"], "INFO")
    assert expected == result


def test_format_value_empty_list():
    """Test format_value with an empty list for a field with number != 1"""
    expected = "."
    result = writer.format_value(header.FieldInfo("String", 2), [], "INFO")
    assert expected == result
