# -*- coding: utf-8 -*-
"""Tests for vcfpy.header
"""

import sys

import vcfpy
from vcfpy import header

import pytest


def test_header_field_info():
    """Test the builtin functions of the FieldInfo class"""
    info1 = header.FieldInfo("Integer", 1, "Some description")
    info2 = header.FieldInfo("Integer", 1, "Some description")
    info3 = header.FieldInfo("Integer", ".", "Some description")
    assert info1 == info2
    assert info1 != info3
    assert hash(info1) == hash(info2)
    assert str(info1) == "FieldInfo('Integer', 1, 'Some description', None)"
    assert repr(info1) == "FieldInfo('Integer', 1, 'Some description', None)"


def test_sample_infos():
    info1 = header.SamplesInfos(["one", "two", "three"])
    info2 = header.SamplesInfos(["one", "two", "three"])
    info3 = header.SamplesInfos(["one", "two", "four"])
    assert info1 == info2
    assert info1 != info3
    with pytest.raises(TypeError):
        assert hash(info1)
    assert (
        str(info1)
        == "SamplesInfos(names=['one', 'two', 'three'], name_to_idx={'one': 0, 'three': 2, 'two': 1})"
    )
    assert (
        repr(info1)
        == "SamplesInfos(names=['one', 'two', 'three'], name_to_idx={'one': 0, 'three': 2, 'two': 1})"
    )


def test_header_header():
    lines1 = [header.HeaderLine("foo", "bar"), header.HeaderLine("foo2", "bar2")]
    samples1 = header.SamplesInfos(["one", "two", "three"])
    hdr1 = header.Header(lines1, samples1)

    lines2 = [header.HeaderLine("foo", "bar"), header.HeaderLine("foo2", "bar2")]
    samples2 = header.SamplesInfos(["one", "two", "three"])
    hdr2 = header.Header(lines2, samples2)

    lines3 = [header.HeaderLine("foo3", "bar"), header.HeaderLine("foo2", "bar2")]
    samples3 = header.SamplesInfos(["one", "two", "three"])
    hdr3 = header.Header(lines3, samples3)

    assert hdr1 == hdr2
    assert hdr1 != hdr3
    EXPECTED = (
        "Header(lines=[HeaderLine('foo', 'bar'), HeaderLine('foo2', 'bar2')], "
        "samples=SamplesInfos(names=['one', 'two', 'three'], "
        "name_to_idx={'one': 0, 'three': 2, 'two': 1}))"
    )
    assert str(hdr1) == EXPECTED

    with pytest.raises(TypeError):
        hash(hdr1)


def test_header_without_lines():
    lines = [header.HeaderLine("foo", "bar"), header.HeaderLine("foo2", "bar2")]
    samples = header.SamplesInfos(["one", "two", "three"])
    hdr = header.Header(lines, samples)
    hdr.add_filter_line(vcfpy.OrderedDict([("ID", "PASS")]))
    hdr.add_filter_line(vcfpy.OrderedDict([("ID", "q30")]))
    assert len(hdr.lines) == 4

    hdr2 = header.header_without_lines(hdr, [("foo", "bar"), ("FILTER", "q30")])
    assert len(hdr2.lines) == 2
    assert hdr2.samples == hdr.samples


def test_header_header_line():
    line1 = header.HeaderLine("key", "value")
    line2 = header.HeaderLine("key", "value")
    line3 = header.HeaderLine("key2", "value")
    assert line1 == line2
    assert line1 != line3
    assert str(line1) == "HeaderLine('key', 'value')"
    assert repr(line1) == "HeaderLine('key', 'value')"
    assert line1.value == "value"
    assert line1.serialize() == "##key=value"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_alt_allele_header_line():
    line1 = header.AltAlleleHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "DEL"), ("Description", "deletion")])
    )
    line2 = header.AltAlleleHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "DEL"), ("Description", "deletion")])
    )
    line3 = header.AltAlleleHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "DUP"), ("Description", "duplication")])
    )
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            """AltAlleleHeaderLine('ALT', '<ID=DEL,Description="deletion">', """
            """OrderedDict([('ID', 'DEL'), ('Description', 'deletion')]))"""
        )
        assert repr(line1) == (
            """AltAlleleHeaderLine('ALT', '<ID=DEL,Description="deletion">', """
            """OrderedDict([('ID', 'DEL'), ('Description', 'deletion')]))"""
        )
    else:
        assert str(line1) == (
            "AltAlleleHeaderLine('ALT', '<ID=DEL,Description=\"deletion\">', "
            "{'ID': 'DEL', 'Description': 'deletion'})"
        )
        assert repr(line1) == (
            "AltAlleleHeaderLine('ALT', '<ID=DEL,Description=\"deletion\">', "
            "{'ID': 'DEL', 'Description': 'deletion'})"
        )
    assert line1.value == '<ID=DEL,Description="deletion">'
    assert line1.serialize() == '##ALT=<ID=DEL,Description="deletion">'
    with pytest.raises(TypeError):
        hash(line1)


def test_header_contig_header_line():
    line1 = header.ContigHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "1"), ("length", 234)]))
    line2 = header.ContigHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "1"), ("length", 234)]))
    line3 = header.ContigHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "2"), ("length", 123)]))
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            "ContigHeaderLine('contig', '<ID=1,length=234>', OrderedDict([('ID', '1'), ('length', 234)]))"
        )
        assert repr(line1) == (
            "ContigHeaderLine('contig', '<ID=1,length=234>', OrderedDict([('ID', '1'), ('length', 234)]))"
        )
    else:
        assert str(line1) == (
            "ContigHeaderLine('contig', '<ID=1,length=234>', {'ID': '1', 'length': 234})"
        )
        assert repr(line1) == (
            "ContigHeaderLine('contig', '<ID=1,length=234>', {'ID': '1', 'length': 234})"
        )
    assert line1.value == "<ID=1,length=234>"
    assert line1.serialize() == "##contig=<ID=1,length=234>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_filter_header_line():
    line1 = header.FilterHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "PASS"), ("Description", "All filters passed")])
    )
    line2 = header.FilterHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "PASS"), ("Description", "All filters passed")])
    )
    line3 = header.FilterHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "q30"), ("Description", "Phred score <30")])
    )
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            "FilterHeaderLine('FILTER', '<ID=PASS,Description=\"All filters passed\">', "
            "OrderedDict([('ID', 'PASS'), ('Description', 'All filters passed')]))"
        )
        assert repr(line1) == (
            "FilterHeaderLine('FILTER', '<ID=PASS,Description=\"All filters passed\">', "
            "OrderedDict([('ID', 'PASS'), ('Description', 'All filters passed')]))"
        )
    else:
        assert str(line1) == (
            "FilterHeaderLine('FILTER', '<ID=PASS,Description=\"All filters passed\">', "
            "{'ID': 'PASS', 'Description': 'All filters passed'})"
        )
        assert repr(line1) == (
            "FilterHeaderLine('FILTER', '<ID=PASS,Description=\"All filters passed\">', "
            "{'ID': 'PASS', 'Description': 'All filters passed'})"
        )
    assert line1.value == '<ID=PASS,Description="All filters passed">'
    assert line1.serialize() == '##FILTER=<ID=PASS,Description="All filters passed">'
    with pytest.raises(TypeError):
        hash(line1)


def test_header_pedigree_header_line():
    line1 = header.PedigreeHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "child"), ("Father", "father")])
    )
    line2 = header.PedigreeHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "child"), ("Father", "father")])
    )
    line3 = header.PedigreeHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "father")]))
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            "PedigreeHeaderLine('PEDIGREE', '<ID=child,Father=father>', "
            "OrderedDict([('ID', 'child'), ('Father', 'father')]))"
        )
        assert repr(line1) == (
            "PedigreeHeaderLine('PEDIGREE', '<ID=child,Father=father>', "
            "OrderedDict([('ID', 'child'), ('Father', 'father')]))"
        )
    else:
        assert str(line1) == (
            "PedigreeHeaderLine('PEDIGREE', '<ID=child,Father=father>', {'ID': 'child', 'Father': 'father'})"
        )
        assert repr(line1) == (
            "PedigreeHeaderLine('PEDIGREE', '<ID=child,Father=father>', {'ID': 'child', 'Father': 'father'})"
        )
    assert line1.value == "<ID=child,Father=father>"
    assert line1.serialize() == "##PEDIGREE=<ID=child,Father=father>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_sample_header_line():
    line1 = header.SampleHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "sample1")]))
    line2 = header.SampleHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "sample1")]))
    line3 = header.SampleHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "sample2")]))
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            "SampleHeaderLine('SAMPLE', '<ID=sample1>', OrderedDict([('ID', 'sample1')]))"
        )
        assert repr(line1) == (
            "SampleHeaderLine('SAMPLE', '<ID=sample1>', OrderedDict([('ID', 'sample1')]))"
        )
    else:
        assert str(line1) == ("SampleHeaderLine('SAMPLE', '<ID=sample1>', {'ID': 'sample1'})")
        assert repr(line1) == ("SampleHeaderLine('SAMPLE', '<ID=sample1>', {'ID': 'sample1'})")
    assert line1.value == "<ID=sample1>"
    assert line1.serialize() == "##SAMPLE=<ID=sample1>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_info_header_line():
    line1 = header.InfoHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "SVTYPE"), ("Number", 1), ("Type", "String")])
    )
    line2 = header.InfoHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "SVTYPE"), ("Number", 1), ("Type", "String")])
    )
    line3 = header.InfoHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "END"), ("Number", 1), ("Type", "Integer")])
    )
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            "InfoHeaderLine('INFO', '<ID=SVTYPE,Number=1,Type=String>', "
            "OrderedDict([('ID', 'SVTYPE'), ('Number', 1), ('Type', 'String')]))"
        )
        assert repr(line1) == (
            "InfoHeaderLine('INFO', '<ID=SVTYPE,Number=1,Type=String>', "
            "OrderedDict([('ID', 'SVTYPE'), ('Number', 1), ('Type', 'String')]))"
        )
    else:
        assert str(line1) == (
            "InfoHeaderLine('INFO', '<ID=SVTYPE,Number=1,Type=String>', "
            "{'ID': 'SVTYPE', 'Number': 1, 'Type': 'String'})"
        )
        assert repr(line1) == (
            "InfoHeaderLine('INFO', '<ID=SVTYPE,Number=1,Type=String>', "
            "{'ID': 'SVTYPE', 'Number': 1, 'Type': 'String'})"
        )
    assert line1.value == "<ID=SVTYPE,Number=1,Type=String>"
    assert line1.serialize() == "##INFO=<ID=SVTYPE,Number=1,Type=String>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_format_header_line():
    line1 = header.FormatHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "AD"), ("Number", "R"), ("Type", "Integer")])
    )
    line2 = header.FormatHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "AD"), ("Number", "R"), ("Type", "Integer")])
    )
    line3 = header.FormatHeaderLine.from_mapping(
        vcfpy.OrderedDict([("ID", "DP"), ("Number", 1), ("Type", "Integer")])
    )
    assert line1 == line2
    assert line1 != line3
    if sys.version_info < (3, 6):
        assert str(line1) == (
            "FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer>', "
            "OrderedDict([('ID', 'AD'), ('Number', 'R'), ('Type', 'Integer')]))"
        )
        assert repr(line1) == (
            "FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer>', "
            "OrderedDict([('ID', 'AD'), ('Number', 'R'), ('Type', 'Integer')]))"
        )
    else:
        assert str(line1) == (
            "FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer>', "
            "{'ID': 'AD', 'Number': 'R', 'Type': 'Integer'})"
        )
        assert repr(line1) == (
            "FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer>', "
            "{'ID': 'AD', 'Number': 'R', 'Type': 'Integer'})"
        )
    assert line1.value == "<ID=AD,Number=R,Type=Integer>"
    assert line1.serialize() == "##FORMAT=<ID=AD,Number=R,Type=Integer>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_has_header_line_positive():
    lines = [
        header.FormatHeaderLine.from_mapping(
            vcfpy.OrderedDict([("ID", "DP"), ("Number", "R"), ("Type", "Integer")])
        ),
        header.InfoHeaderLine.from_mapping(
            vcfpy.OrderedDict([("ID", "AD"), ("Number", "R"), ("Type", "Integer")])
        ),
        header.FilterHeaderLine.from_mapping(
            vcfpy.OrderedDict([("ID", "PASS"), ("Description", "All filters passed")])
        ),
        header.ContigHeaderLine.from_mapping(vcfpy.OrderedDict([("ID", "1"), ("length", 234)])),
    ]
    samples = header.SamplesInfos(["one", "two", "three"])
    hdr = header.Header(lines, samples)

    assert hdr.has_header_line("FORMAT", "DP")
    assert hdr.has_header_line("INFO", "AD")
    assert hdr.has_header_line("FILTER", "PASS")
    assert hdr.has_header_line("contig", "1")


def test_header_has_header_line_positive_no_samples():
    lines = []
    samples = header.SamplesInfos(["one", "two", "three"])
    hdr = header.Header(lines, samples)

    assert not hdr.has_header_line("FORMAT", "DP")
    assert not hdr.has_header_line("INFO", "AD")
    assert not hdr.has_header_line("FILTER", "PASS")
    assert not hdr.has_header_line("contig", "1")
