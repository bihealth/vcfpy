# -*- coding: utf-8 -*-
"""Test that the public API classes are available from the ``vcfpy`` package
"""

import io
import textwrap

import pytest

from vcfpy import exceptions
from vcfpy import header
from vcfpy import parser


# WarningHelper ---------------------------------------------------------------


def test_warning_helper(capsys):
    stream = io.StringIO()
    helper = parser.WarningHelper(stream=stream)
    helper.warn_once('foo')
    helper.warn_once('bar')
    helper.warn_once('foo')

    EXPECTED = textwrap.dedent(r"""
    [vcfpy] foo
    (Subsequent identical messages will not be printed)
    [vcfpy] bar
    (Subsequent identical messages will not be printed)
    """).lstrip()
    assert stream.getvalue() == EXPECTED

    helper.print_summary()

    EXPECTED = textwrap.dedent("""
    [vcfpy] foo
    (Subsequent identical messages will not be printed)
    [vcfpy] bar
    (Subsequent identical messages will not be printed)
    WARNINGS

         2\tfoo
         1\tbar
    """).lstrip()
    assert stream.getvalue() == EXPECTED

# HeaderChecker ---------------------------------------------------------------


def test_header_checker_missing_fileformat():
    # construct Header
    lines = [
        header.HeaderLine('key', 'value'),
    ]
    hdr = header.Header(lines=lines, samples=header.SamplesInfos(['NA001']))
    # setup HeaderChecker
    stream = io.StringIO()
    helper = parser.WarningHelper(stream=stream)
    checker = parser.HeaderChecker(helper)
    # execute
    with pytest.raises(exceptions.InvalidHeaderException):
        checker.run(hdr)


def test_header_checker_unknown_vcf_version():
    # construct Header
    lines = [
        header.HeaderLine('fileformat', 'v4.3'),
    ]
    hdr = header.Header(lines=lines, samples=header.SamplesInfos(['NA001']))
    # setup HeaderChecker
    stream = io.StringIO()
    helper = parser.WarningHelper(stream=stream)
    checker = parser.HeaderChecker(helper)
    # execute
    checker.run(hdr)
    # check result
    EXPECTED = textwrap.dedent(r"""
    [vcfpy] WARNING: The VCF version v4.3 is not known, going on regardlessly
    (Subsequent identical messages will not be printed)
    """).lstrip()
    assert stream.getvalue() == EXPECTED


def test_header_checker_known_vcf_version():
    # construct Header
    lines = [
        header.HeaderLine('fileformat', 'VCFv4.3'),
    ]
    hdr = header.Header(lines=lines, samples=header.SamplesInfos(['NA001']))
    # setup HeaderChecker
    stream = io.StringIO()
    helper = parser.WarningHelper(stream=stream)
    checker = parser.HeaderChecker(helper)
    # execute
    checker.run(hdr)
    # check result
    EXPECTED = ''
    assert stream.getvalue() == EXPECTED
