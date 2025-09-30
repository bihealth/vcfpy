# -*- coding: utf-8 -*-
"""Additional tests to achieve 100% coverage in parser.py"""

import warnings

from vcfpy import header, parser, record

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def test_process_sub_shrink_edge_cases():
    """Test process_sub_shrink function edge cases - line 330"""
    # Test case where ref has one character and equals alt_str[0] (INS)
    result = parser.process_sub_shrink("A", "AT")
    assert isinstance(result, record.Substitution)
    assert result.type == record.INS
    assert result.value == "AT"

    # Test case where ref has one character but doesn't equal alt_str[0] (INDEL)
    result = parser.process_sub_shrink("A", "CT")
    assert isinstance(result, record.Substitution)
    assert result.type == record.INDEL
    assert result.value == "CT"

    # Test case where ref has multiple characters (INDEL)
    result = parser.process_sub_shrink("ATG", "C")
    assert isinstance(result, record.Substitution)
    assert result.type == record.INDEL
    assert result.value == "C"


def test_record_parser_format_checkers():
    """Test RecordParser with FORMAT checking enabled - lines 418, 423"""
    # Create a header with FORMAT fields
    hdr = header.Header(
        [
            header.HeaderLine("fileformat", "VCFv4.3"),
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "GT", "Number": "1", "Type": "String", "Description": "Genotype"}
            ),
        ]
    )

    # Create SamplesInfos
    samples = header.SamplesInfos(["sample1"])

    # Create RecordParser with FORMAT checking enabled
    parser_obj = parser.RecordParser(hdr, samples, record_checks=["FORMAT"])

    # Verify FORMAT checker is set up
    assert isinstance(parser_obj._format_checker, parser.FormatChecker)
    assert not isinstance(parser_obj._format_checker, parser.NoopFormatChecker)


def test_record_parser_info_checkers():
    """Test RecordParser with INFO checking enabled - line 418"""
    # Create a header with INFO fields
    hdr = header.Header(
        [
            header.HeaderLine("fileformat", "VCFv4.3"),
            header.InfoHeaderLine(
                "INFO", "", {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Total depth"}
            ),
        ]
    )

    # Create SamplesInfos
    samples = header.SamplesInfos(["sample1"])

    # Create RecordParser with INFO checking enabled
    parser_obj = parser.RecordParser(hdr, samples, record_checks=["INFO"])

    # Verify INFO checker is set up
    assert isinstance(parser_obj._info_checker, parser.InfoChecker)
    assert not isinstance(parser_obj._info_checker, parser.NoopInfoChecker)


def test_binomial_edge_cases():
    """Test binomial function edge cases - lines 602-606"""
    # Test normal cases
    assert parser.binomial(5, 2) == 10
    assert parser.binomial(4, 0) == 1
    assert parser.binomial(4, 4) == 1

    # Test edge case that triggers ValueError (negative values)
    assert parser.binomial(-1, 2) == 0
    assert parser.binomial(2, -1) == 0
    assert parser.binomial(-1, -1) == 0


def test_info_checker_run():
    """Test InfoChecker.run method - line 629"""
    # Create header with INFO field
    hdr = header.Header(
        [
            header.InfoHeaderLine(
                "INFO", "", {"ID": "AF", "Number": "A", "Type": "Float", "Description": "Allele frequency"}
            ),
        ]
    )

    checker = parser.InfoChecker(hdr)

    # Test with non-list value (should return early)
    checker.run("AF", "0.5", 1)  # This should not raise any warnings

    # Test with list value and correct count
    checker.run("AF", [0.5], 1)  # 1 alt allele, expecting 1 value - should be OK


def test_info_checker_incorrect_count():
    """Test InfoChecker count checking - lines 640-652"""
    # Create header with INFO fields
    hdr = header.Header(
        [
            header.InfoHeaderLine(
                "INFO", "", {"ID": "AF", "Number": "A", "Type": "Float", "Description": "Allele frequency"}
            ),
            header.InfoHeaderLine(
                "INFO", "", {"ID": "AC", "Number": "R", "Type": "Integer", "Description": "Allele count"}
            ),
            header.InfoHeaderLine(
                "INFO", "", {"ID": "AN", "Number": "G", "Type": "Integer", "Description": "Allele number"}
            ),
            header.InfoHeaderLine(
                "INFO", "", {"ID": "CUSTOM", "Number": "2", "Type": "Integer", "Description": "Custom field"}
            ),
        ]
    )

    checker = parser.InfoChecker(hdr)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        # Test A number: 2 alt alleles but 1 value given (should warn)
        checker.run("AF", [0.5], 2)

        # Test R number: 2 alt alleles (expecting 3 values) but 2 values given (should warn)
        checker.run("AC", [10, 5], 2)

        # Test G number: 2 alt alleles (expecting 3 values for diploid) but 2 values given (should warn)
        checker.run("AN", [20, 10], 2)

        # Test fixed number: expecting 2 values but 1 given (should warn)
        checker.run("CUSTOM", [100], 2)

        assert len(w) == 4
        for warning in w:
            assert "Number of elements for INFO field" in str(warning.message)


def test_format_checker_init():
    """Test FormatChecker initialization - line 674"""
    hdr = header.Header(
        [
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "GT", "Number": "1", "Type": "String", "Description": "Genotype"}
            ),
        ]
    )

    checker = parser.FormatChecker(hdr)
    assert checker.header == hdr


def test_format_checker_run():
    """Test FormatChecker.run method - lines 681-682"""
    hdr = header.Header(
        [
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "GT", "Number": "1", "Type": "String", "Description": "Genotype"}
            ),
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Read depth"}
            ),
        ]
    )

    checker = parser.FormatChecker(hdr)

    # Create a call with values that won't trigger warnings
    call = record.Call("sample1", {"GT": "0", "DP": "30"})  # Single values for Number=1
    call.gt_alleles = [0, 1]  # Set gt_alleles

    # This should run without issues for string values
    checker.run(call, 1)


def test_format_checker_count_validation():
    """Test FormatChecker._check_count method - lines 685-698"""
    hdr = header.Header(
        [
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "AD", "Number": "R", "Type": "Integer", "Description": "Allelic depths"}
            ),
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "PL", "Number": "G", "Type": "Integer", "Description": "Phred-scaled likelihoods"}
            ),
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "CUSTOM", "Number": "2", "Type": "Integer", "Description": "Custom field"}
            ),
        ]
    )

    checker = parser.FormatChecker(hdr)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        # Test case that should trigger the warning path
        call1 = record.Call("sample1", {"CUSTOM": "1"})  # Fixed number: expecting 2 but 1 given
        call1.gt_alleles = [0, 1]
        checker.run(call1, 1)

        # Test another case with different field
        call2 = record.Call("sample2", {"AD": "1"})  # R number: 1 alt expecting 2 but 1 given
        call2.gt_alleles = [0, 1]
        checker.run(call2, 1)

        # Check that we got at least one warning
        assert len(w) >= 1
        for warning in w:
            assert "Number of elements for FORMAT field" in str(warning.message)


def test_format_checker_list_value_skip():
    """Test FormatChecker with list values that should be skipped - line 687"""
    hdr = header.Header(
        [
            header.FormatHeaderLine(
                "FORMAT", "", {"ID": "AD", "Number": "R", "Type": "Integer", "Description": "Allelic depths"}
            ),
        ]
    )

    checker = parser.FormatChecker(hdr)

    # Create a call with list value (should be skipped in _check_count)
    call = record.Call("sample1", {"AD": [10, 5]})  # This is a list
    call.gt_alleles = [0, 1]

    # This should not generate warnings since list values are skipped
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        checker.run(call, 1)
        # Should have no warnings since list values are returned early
        warning_messages = [str(warning.message) for warning in w]
        format_warnings = [msg for msg in warning_messages if "Number of elements for FORMAT field" in msg]
        assert len(format_warnings) == 0
