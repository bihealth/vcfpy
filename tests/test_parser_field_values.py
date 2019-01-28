# -*- coding: utf-8 -*-
"""Parsing and converting of field values
"""

import pytest

from vcfpy import header
from vcfpy import parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# parser.convert_field_value() ------------------------------------------------


def test_convert_field_values_integer_value():
    EXPECTED = 42
    RESULT = parser.convert_field_value("Integer", "42")
    assert EXPECTED == RESULT


def test_convert_field_values_integer_none():
    EXPECTED = None
    RESULT = parser.convert_field_value("Integer", ".")
    assert EXPECTED == RESULT


def test_convert_field_values_float_value():
    EXPECTED = 42.0
    RESULT = parser.convert_field_value("Float", "42")
    assert EXPECTED == RESULT


def test_convert_field_values_float_none():
    EXPECTED = None
    RESULT = parser.convert_field_value("Float", ".")
    assert EXPECTED == RESULT


def test_convert_field_values_flag():
    EXPECTED = True
    RESULT = parser.convert_field_value("Flag", True)
    assert EXPECTED == RESULT


def test_convert_field_values_character_value():
    EXPECTED = "X"
    RESULT = parser.convert_field_value("Character", "X")
    assert EXPECTED == RESULT


def test_convert_field_values_character_none():
    EXPECTED = None
    RESULT = parser.convert_field_value("Character", ".")
    assert EXPECTED == RESULT


def test_convert_field_values_string_value():
    EXPECTED = "Value"
    RESULT = parser.convert_field_value("String", "Value")
    assert EXPECTED == RESULT


def test_convert_field_values_string_none():
    EXPECTED = None
    RESULT = parser.convert_field_value("String", ".")
    assert EXPECTED == RESULT


# parser.parse_field_value() --------------------------------------------------


def test_parse_field_value_integer_one():
    EXPECTED = 42
    RESULT = parser.parse_field_value(header.FieldInfo("Integer", 1), "42")
    assert EXPECTED == RESULT


@pytest.mark.parametrize("number", [2, "A", "R", "G", "."])
def test_parse_field_value_integer_more(number):
    EXPECTED = [42, 43]
    RESULT = parser.parse_field_value(header.FieldInfo("Integer", number), "42,43")
    assert EXPECTED == RESULT


@pytest.mark.parametrize("number", [2, "A", "R", "G", "."])
def test_parse_field_value_integer_more_empty(number):
    EXPECTED = []
    RESULT = parser.parse_field_value(header.FieldInfo("Integer", number), ".")
    assert EXPECTED == RESULT


def test_parse_field_value_string_one():
    EXPECTED = "42"
    RESULT = parser.parse_field_value(header.FieldInfo("String", 1), "42")
    assert EXPECTED == RESULT


@pytest.mark.parametrize("number", [2, "A", "R", "G", "."])
def test_parse_field_value_string_more(number):
    EXPECTED = ["42", "43"]
    RESULT = parser.parse_field_value(header.FieldInfo("String", number), "42,43")
    assert EXPECTED == RESULT


@pytest.mark.parametrize("number", [2, "A", "R", "G", "."])
def test_parse_field_value_string_more_empty(number):
    EXPECTED = []
    RESULT = parser.parse_field_value(header.FieldInfo("String", number), ".")
    assert EXPECTED == RESULT
