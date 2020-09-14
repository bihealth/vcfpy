# -*- coding: utf-8 -*-
"""Test escaping on writing VCF records
"""

from vcfpy import writer
from vcfpy import header

# vcfpy.writer.format_atomic() ------------------------------------------------


def test_format_atomic_with_escape_info():
    EXPECTED = "%3B%3D%25%2C%0D%0A%09"
    RESULT = writer.format_atomic(";=%,\r\n\t", "INFO")
    assert EXPECTED == RESULT


def test_format_atomic_with_escape_format():
    EXPECTED = "%3A%3D%25%2C%0D%0A%09"
    RESULT = writer.format_atomic(":=%,\r\n\t", "FORMAT")
    assert EXPECTED == RESULT


def test_format_atomic_without_escape_info():
    EXPECTED = "This is a legal string:"
    RESULT = writer.format_atomic("This is a legal string:", "INFO")
    assert EXPECTED == RESULT


def test_format_atomic_without_escape_format():
    EXPECTED = "This is a legal string;"
    RESULT = writer.format_atomic("This is a legal string;", "FORMAT")
    assert EXPECTED == RESULT


# vcfpy.writer.format_value() -------------------------------------------------


def test_format_value_with_escape():
    EXPECTED = "%3A%3B%3D%25%2C%0D%0A%09,%25"
    RESULT = writer.format_value(header.FieldInfo("String", 2), (":;=%,\r\n\t", "%"), "INFO")
    assert EXPECTED == RESULT


def test_format_value_without_escape():
    EXPECTED = "This is a legal string,me too"
    RESULT = writer.format_value(
        header.FieldInfo("String", 2), ("This is a legal string", "me too"), "INFO"
    )
    assert EXPECTED == RESULT
