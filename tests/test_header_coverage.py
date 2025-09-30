# -*- coding: utf-8 -*-
"""Additional tests to achieve 100% coverage in header.py"""

import warnings

import pytest

from vcfpy import exceptions, header

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def test_header_filter_ids():
    """Test Header.filter_ids() method - line 397"""
    lines = [
        header.FilterHeaderLine("FILTER", "Description=Test filter 1", {"ID": "q10", "Description": "Test filter 1"}),
        header.FilterHeaderLine("FILTER", "Description=Test filter 2", {"ID": "s50", "Description": "Test filter 2"}),
    ]
    hdr = header.Header(lines)
    filter_ids = hdr.filter_ids()
    assert set(filter_ids) == {"q10", "s50"}


def test_header_info_ids():
    """Test Header.info_ids() method - line 405"""
    lines = [
        header.InfoHeaderLine(
            "INFO",
            "Description=Test info 1",
            {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Test info 1"},
        ),
        header.InfoHeaderLine(
            "INFO",
            "Description=Test info 2",
            {"ID": "AF", "Number": "A", "Type": "Float", "Description": "Test info 2"},
        ),
    ]
    hdr = header.Header(lines)
    info_ids = hdr.info_ids()
    assert set(info_ids) == {"DP", "AF"}


def test_header_get_lines():
    """Test Header.get_lines() method - lines 409-412"""
    lines = [
        header.InfoHeaderLine(
            "INFO", "Description=Test info", {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Test info"}
        ),
        header.FilterHeaderLine("FILTER", "Description=Test filter", {"ID": "q10", "Description": "Test filter"}),
    ]
    hdr = header.Header(lines)

    # Test getting INFO lines
    info_lines = list(hdr.get_lines("INFO"))
    assert len(info_lines) == 1
    assert info_lines[0].mapping["ID"] == "DP"

    # Test getting FILTER lines
    filter_lines = list(hdr.get_lines("FILTER"))
    assert len(filter_lines) == 1
    assert filter_lines[0].mapping["ID"] == "q10"

    # Test getting non-existent type
    other_lines = list(hdr.get_lines("NONEXISTENT"))
    assert len(other_lines) == 0


def test_header_has_header_line_false():
    """Test Header.has_header_line() returning False - line 425"""
    hdr = header.Header([])

    # Test with non-existent key
    assert not hdr.has_header_line("NONEXISTENT", "any_id")

    # Test with existing key but non-existent ID
    lines = [
        header.InfoHeaderLine(
            "INFO", "Description=Test", {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Test"}
        )
    ]
    hdr = header.Header(lines)
    assert not hdr.has_header_line("INFO", "NONEXISTENT_ID")


def test_header_field_info_warning():
    """Test _get_field_info warning path - line 470"""
    hdr = header.Header([])

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = hdr._get_field_info("INFO", "UNKNOWN_KEY", {})

        assert len(w) == 1
        assert "UNKNOWN_KEY not found" in str(w[0].message)
        assert result.type == "String"
        assert result.number == "."


def test_header_repr():
    """Test Header.__repr__() method - line 495"""
    lines = [header.HeaderLine("test", "value")]
    samples = header.SamplesInfos(["sample1"])
    hdr = header.Header(lines, samples)

    repr_str = repr(hdr)
    assert "Header(" in repr_str
    assert "lines=" in repr_str
    assert "samples=" in repr_str


def test_simple_header_line_missing_id():
    """Test SimpleHeaderLine (via AltAlleleHeaderLine) with missing ID - line 564"""
    with pytest.raises(exceptions.InvalidHeaderException) as exc_info:
        # Create mapping without ID - use AltAlleleHeaderLine which inherits from SimpleHeaderLine
        mapping = {"Description": "Test alt allele"}
        header.AltAlleleHeaderLine("ALT", "value", mapping)

    assert 'Missing key "ID"' in str(exc_info.value)


def test_simple_header_line_eq_false():
    """Test SimpleHeaderLine.__eq__() returning False - line 581"""
    line1 = header.SimpleHeaderLine("test", "value1", {"ID": "test1"})
    line2 = header.SimpleHeaderLine("test", "value2", {"ID": "test2"})
    line3 = header.SimpleHeaderLine("different", "value1", {"ID": "test1"})

    assert line1 != line2
    assert line1 != line3
    assert line1 != "not a header line"


def test_compound_header_line_invalid_number_warning():
    """Test CompoundHeaderLine with invalid Number field - lines 770-775"""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        mapping = {"ID": "test", "Number": "invalid", "Type": "Integer", "Description": "Test"}
        line = header.CompoundHeaderLine("INFO", "value", mapping)

        assert len(w) == 1
        assert "invalid number invalid" in str(w[0].message)
        assert line.mapping["Number"] == "."


def test_compound_header_line_parse_number_valid():
    """Test CompoundHeaderLine._parse_number with valid special numbers - line 794"""
    line = header.CompoundHeaderLine(
        "INFO", "value", {"ID": "test", "Number": "A", "Type": "Integer", "Description": "Test"}
    )

    # Test that A, R, G are preserved
    assert line._parse_number("A") == "A"
    assert line._parse_number("R") == "R"
    assert line._parse_number("G") == "G"
    assert line._parse_number(".") == "."


def test_compound_header_line_str():
    """Test CompoundHeaderLine.__str__() method - line 804"""
    mapping = {"ID": "test", "Number": "1", "Type": "Integer", "Description": "Test"}
    line = header.CompoundHeaderLine("INFO", "value", mapping)

    str_repr = str(line)
    assert "CompoundHeaderLine(" in str_repr
    assert "'INFO'" in str_repr
    assert "test" in str_repr


def test_info_header_line_missing_type_warning():
    """Test InfoHeaderLine with missing Type field - lines 828-832"""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        mapping = {"ID": "test", "Number": "1", "Description": "Test"}
        line = header.InfoHeaderLine("INFO", "value", mapping)

        assert len(w) == 1
        assert 'Field "Type" not found' in str(w[0].message)
        assert line.type == "String"


def test_info_header_line_invalid_type_warning():
    """Test InfoHeaderLine with invalid Type field - lines 834-840"""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        mapping = {"ID": "test", "Number": "1", "Type": "InvalidType", "Description": "Test"}
        line = header.InfoHeaderLine("INFO", "value", mapping)

        assert len(w) == 1
        assert "Invalid INFO value type InvalidType" in str(w[0].message)
        assert line.type == "String"


def test_info_header_line_missing_description_warning():
    """Test InfoHeaderLine with missing Description field - line 845"""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        mapping = {"ID": "test", "Number": "1", "Type": "Integer"}
        line = header.InfoHeaderLine("INFO", "value", mapping)

        assert len(w) == 1
        assert 'Field "Description" not found' in str(w[0].message)
        assert line.description is None


def test_header_format_ids():
    """Test Header.format_ids() method - line 397"""
    lines = [
        header.FormatHeaderLine(
            "FORMAT",
            "Description=Test format 1",
            {"ID": "GT", "Number": "1", "Type": "String", "Description": "Test format 1"},
        ),
        header.FormatHeaderLine(
            "FORMAT",
            "Description=Test format 2",
            {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Test format 2"},
        ),
    ]
    hdr = header.Header(lines)
    format_ids = hdr.format_ids()
    assert set(format_ids) == {"GT", "DP"}


def test_simple_header_line_str():
    """Test SimpleHeaderLine.__str__() method - line 581"""
    mapping = {"ID": "test", "Description": "Test description"}
    line = header.SimpleHeaderLine("TEST", "value", mapping)

    str_repr = str(line)
    assert "SimpleHeaderLine(" in str_repr
    assert "'TEST'" in str_repr
    assert "test" in str_repr
    assert "Test description" in str_repr
