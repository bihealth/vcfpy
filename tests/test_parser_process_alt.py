# -*- coding: utf-8 -*-
"""Test procesing of ALT column
"""

import pytest

from vcfpy import header
from vcfpy import parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# TODO: add more tests


@pytest.fixture
def vcf_header():
    return header.Header()


# result type Substitution ----------------------------------------------------


def test_substitution_snv(vcf_header):
    EXPECTED = "Substitution(type_='SNV', value='T')"
    RESULT = parser.process_alt(vcf_header, "C", "T")
    assert str(RESULT) == EXPECTED


def test_substitution_mnv(vcf_header):
    EXPECTED = "Substitution(type_='MNV', value='TG')"
    RESULT = parser.process_alt(vcf_header, "GC", "TG")
    assert str(RESULT) == EXPECTED


def test_substitution_ins(vcf_header):
    EXPECTED = "Substitution(type_='INS', value='CT')"
    RESULT = parser.process_alt(vcf_header, "C", "CT")
    assert str(RESULT) == EXPECTED


def test_substitution_del(vcf_header):
    EXPECTED = "Substitution(type_='DEL', value='C')"
    RESULT = parser.process_alt(vcf_header, "CT", "C")
    assert str(RESULT) == EXPECTED


def test_substitution_indel(vcf_header):
    EXPECTED = "Substitution(type_='INDEL', value='TG')"
    RESULT = parser.process_alt(vcf_header, "AAAC", "TG")
    assert str(RESULT) == EXPECTED
