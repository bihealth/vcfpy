# -*- coding: utf-8 -*-
"""Tests for vcfpy.header"""

import warnings

import pytest

from vcfpy import exceptions, header
from vcfpy.header import (
    AltAlleleHeaderLine,
    CompoundHeaderLine,
    ContigHeaderLine,
    FieldInfo,
    FilterHeaderLine,
    FormatHeaderLine,
    Header,
    HeaderInvalidType,
    HeaderLine,
    HeaderMissingDescription,
    MetaHeaderLine,
    PedigreeHeaderLine,
    SampleHeaderLine,
    SamplesInfos,
    header_without_lines,
    mapping_to_str,
    serialize_for_header,
)


def test_header_remove_lines_edge_cases():
    """Test header removal with edge cases"""
    # Create header with various line types
    lines = [
        HeaderLine("fileformat", "VCFv4.2"),
        CompoundHeaderLine(
            "INFO",
            '<ID=DP,Number=1,Type=Integer,Description="Depth">',
            {"ID": "DP", "Number": "1", "Type": "Integer", "Description": "Depth"},
        ),
        HeaderLine("reference", "test.fa"),
    ]

    test_header = header.Header(lines)

    # Remove lines that don't have mapping
    result = header_without_lines(test_header, [("reference", "test.fa")])
    assert len(result.lines) == 2


def test_parser_field_info_edge_cases():
    """Test FieldInfo with various edge cases"""
    # Test with missing number field
    field_info = FieldInfo("TEST", None, "String", "Test")
    assert field_info.number is None

    # FieldInfo should be hashable
    hash(field_info)  # Should not raise


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
    assert str(info1) == "SamplesInfos(names=['one', 'two', 'three'], name_to_idx={'one': 0, 'three': 2, 'two': 1})"
    assert repr(info1) == "SamplesInfos(names=['one', 'two', 'three'], name_to_idx={'one': 0, 'three': 2, 'two': 1})"


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
    hdr.add_filter_line({"ID": "PASS", "Description": "All filters passed"})
    hdr.add_filter_line({"ID": "q30", "Description": "Phred score <30"})
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
    line1 = header.AltAlleleHeaderLine.from_mapping({"ID": "DEL", "Description": "deletion"})
    line2 = header.AltAlleleHeaderLine.from_mapping({"ID": "DEL", "Description": "deletion"})
    line3 = header.AltAlleleHeaderLine.from_mapping({"ID": "DUP", "Description": "duplication"})
    assert line1 == line2
    assert line1 != line3
    assert (
        str(line1) == "AltAlleleHeaderLine('ALT', '<ID=DEL,Description=\"deletion\">', "
        "{'ID': 'DEL', 'Description': 'deletion'})"
    )
    assert (
        repr(line1) == "AltAlleleHeaderLine('ALT', '<ID=DEL,Description=\"deletion\">', "
        "{'ID': 'DEL', 'Description': 'deletion'})"
    )
    assert line1.value == '<ID=DEL,Description="deletion">'
    assert line1.serialize() == '##ALT=<ID=DEL,Description="deletion">'
    with pytest.raises(TypeError):
        hash(line1)


def test_header_contig_header_line():
    line1 = header.ContigHeaderLine.from_mapping({"ID": "1", "length": 234})
    line2 = header.ContigHeaderLine.from_mapping({"ID": "1", "length": 234})
    line3 = header.ContigHeaderLine.from_mapping({"ID": "2", "length": 123})
    assert line1 == line2
    assert line1 != line3
    assert str(line1) == "ContigHeaderLine('contig', '<ID=1,length=234>', {'ID': '1', 'length': 234})"
    assert repr(line1) == "ContigHeaderLine('contig', '<ID=1,length=234>', {'ID': '1', 'length': 234})"
    assert line1.value == "<ID=1,length=234>"
    assert line1.serialize() == "##contig=<ID=1,length=234>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_filter_header_line():
    line1 = header.FilterHeaderLine.from_mapping({"ID": "PASS", "Description": "All filters passed"})
    line2 = header.FilterHeaderLine.from_mapping({"ID": "PASS", "Description": "All filters passed"})
    line3 = header.FilterHeaderLine.from_mapping({"ID": "q30", "Description": "Phred score <30"})
    assert line1 == line2
    assert line1 != line3
    assert (
        str(line1) == "FilterHeaderLine('FILTER', '<ID=PASS,Description=\"All filters passed\">', "
        "{'ID': 'PASS', 'Description': 'All filters passed'})"
    )
    assert (
        repr(line1) == "FilterHeaderLine('FILTER', '<ID=PASS,Description=\"All filters passed\">', "
        "{'ID': 'PASS', 'Description': 'All filters passed'})"
    )
    assert line1.value == '<ID=PASS,Description="All filters passed">'
    assert line1.serialize() == '##FILTER=<ID=PASS,Description="All filters passed">'
    with pytest.raises(TypeError):
        hash(line1)


def test_header_pedigree_header_line():
    line1 = header.PedigreeHeaderLine.from_mapping({"ID": "child", "Father": "father"})
    line2 = header.PedigreeHeaderLine.from_mapping({"ID": "child", "Father": "father"})
    line3 = header.PedigreeHeaderLine.from_mapping({"ID": "father"})
    assert line1 == line2
    assert line1 != line3
    assert (
        str(line1) == "PedigreeHeaderLine('PEDIGREE', '<ID=child,Father=father>', {'ID': 'child', 'Father': 'father'})"
    )
    assert (
        repr(line1) == "PedigreeHeaderLine('PEDIGREE', '<ID=child,Father=father>', {'ID': 'child', 'Father': 'father'})"
    )
    assert line1.value == "<ID=child,Father=father>"
    assert line1.serialize() == "##PEDIGREE=<ID=child,Father=father>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_sample_header_line():
    line1 = header.SampleHeaderLine.from_mapping({"ID": "sample1"})
    line2 = header.SampleHeaderLine.from_mapping({"ID": "sample1"})
    line3 = header.SampleHeaderLine.from_mapping({"ID": "sample2"})
    assert line1 == line2
    assert line1 != line3
    assert str(line1) == "SampleHeaderLine('SAMPLE', '<ID=sample1>', {'ID': 'sample1'})"
    assert repr(line1) == "SampleHeaderLine('SAMPLE', '<ID=sample1>', {'ID': 'sample1'})"
    assert line1.value == "<ID=sample1>"
    assert line1.serialize() == "##SAMPLE=<ID=sample1>"
    with pytest.raises(TypeError):
        hash(line1)


def test_header_info_header_line():
    line1 = header.InfoHeaderLine.from_mapping(
        {"ID": "SVTYPE", "Number": 1, "Type": "String", "Description": "Type of structural variant"}
    )
    line2 = header.InfoHeaderLine.from_mapping(
        {"ID": "SVTYPE", "Number": 1, "Type": "String", "Description": "Type of structural variant"}
    )
    line3 = header.InfoHeaderLine.from_mapping(
        {"ID": "END", "Number": 1, "Type": "Integer", "Description": "End position"}
    )
    assert line1 == line2
    assert line1 != line3
    assert (
        str(line1)
        == "InfoHeaderLine('INFO', '<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">', "
        "{'ID': 'SVTYPE', 'Number': 1, 'Type': 'String', 'Description': 'Type of structural variant'})"
    )
    assert (
        repr(line1)
        == "InfoHeaderLine('INFO', '<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">', "
        "{'ID': 'SVTYPE', 'Number': 1, 'Type': 'String', 'Description': 'Type of structural variant'})"
    )
    assert line1.value == '<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
    assert line1.serialize() == '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
    with pytest.raises(TypeError):
        hash(line1)


def test_header_format_header_line():
    line1 = header.FormatHeaderLine.from_mapping(
        {"ID": "AD", "Number": "R", "Type": "Integer", "Description": "Allelic depths"}
    )
    line2 = header.FormatHeaderLine.from_mapping(
        {"ID": "AD", "Number": "R", "Type": "Integer", "Description": "Allelic depths"}
    )
    line3 = header.FormatHeaderLine.from_mapping(
        {"ID": "DP", "Number": 1, "Type": "Integer", "Description": "Read depth"}
    )
    assert line1 == line2
    assert line1 != line3
    assert (
        str(line1) == "FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">', "
        "{'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths'})"
    )
    assert (
        repr(line1) == "FormatHeaderLine('FORMAT', '<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">', "
        "{'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths'})"
    )
    assert line1.value == '<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'
    assert line1.serialize() == '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'
    with pytest.raises(TypeError):
        hash(line1)


def test_header_has_header_line_positive():
    lines = [
        header.FormatHeaderLine.from_mapping(
            {"ID": "DP", "Number": "R", "Type": "Integer", "Description": "Depth of coverage"}
        ),
        header.InfoHeaderLine.from_mapping(
            {"ID": "AD", "Number": "R", "Type": "Integer", "Description": "Allelic depths"}
        ),
        header.FilterHeaderLine.from_mapping({"ID": "PASS", "Description": "All filters passed"}),
        header.ContigHeaderLine.from_mapping({"ID": "1", "length": 234}),
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


def test_header_get_format_field_info():
    lines = []
    samples = header.SamplesInfos(["one", "two", "three"])
    hdr = header.Header(lines, samples)
    with pytest.warns(exceptions.FieldInfoNotFound):
        gt_field_info = hdr.get_format_field_info("GT")

    expected = header.RESERVED_FORMAT["GT"]

    assert gt_field_info is expected


def test_header_get_info_format_field_info():
    lines = []
    samples = header.SamplesInfos(["one", "two", "three"])
    hdr = header.Header(lines, samples)
    with pytest.warns(exceptions.FieldInfoNotFound):
        gt_field_info = hdr.get_info_field_info("AA")

    expected = header.RESERVED_INFO["AA"]

    assert gt_field_info is expected


def test_header_line_equality_and_hash():
    """Test HeaderLine equality and hash functionality"""
    line1 = header.HeaderLine("INFO", "DP")
    line2 = header.HeaderLine("INFO", "DP")
    line3 = header.HeaderLine("FORMAT", "DP")

    # Test equality
    assert line1 == line2
    assert line1 != line3
    assert line1 != "not a header line"

    # Test hash - HeaderLine is intentionally unhashable
    with pytest.raises(TypeError, match="Unhashable type: HeaderLine"):
        hash(line1)


def test_field_info_basic():
    """Test FieldInfo basic functionality"""
    field_info = header.FieldInfo("Integer", 1, "Some description")
    assert field_info.type == "Integer"
    assert field_info.number == 1
    assert field_info.description == "Some description"


def test_samples_infos_hash_fails():
    """Test that SamplesInfos hash raises TypeError (contains unhashable list)"""
    samples = header.SamplesInfos(["one", "two", "three"])
    with pytest.raises(TypeError):
        hash(samples)


def test_samples_infos_inequality():
    """Test SamplesInfos inequality method"""
    samples1 = header.SamplesInfos(["one", "two", "three"])
    samples2 = header.SamplesInfos(["one", "two", "three"])
    samples3 = header.SamplesInfos(["one", "two", "four"])

    # Test __ne__ method
    assert not (samples1 != samples2)  # Should be False
    assert samples1 != samples3  # Should be True
    assert samples1 != "not samples"  # Should be True


def test_header_get_lines_by_key():
    """Test header get_lines functionality for specific keys"""
    lines = [
        header.InfoHeaderLine.from_mapping({"ID": "DP", "Number": 1, "Type": "Integer", "Description": "Total depth"}),
        header.FormatHeaderLine.from_mapping({"ID": "GT", "Number": 1, "Type": "String", "Description": "Genotype"}),
        header.InfoHeaderLine.from_mapping(
            {"ID": "AC", "Number": "A", "Type": "Integer", "Description": "Allele count"}
        ),
    ]
    samples = header.SamplesInfos(["sample1"])
    hdr = header.Header(lines, samples)

    # Get INFO lines
    info_lines = [line for line in hdr.lines if line.key == "INFO"]
    assert len(info_lines) == 2

    # Get FORMAT lines
    format_lines = [line for line in hdr.lines if line.key == "FORMAT"]
    assert len(format_lines) == 1


def test_header_line_hash_raises_error():
    """Test that HeaderLine.__hash__ raises TypeError"""
    line = HeaderLine("test", "value")
    with pytest.raises(TypeError, match="Unhashable type: HeaderLine"):
        hash(line)


def test_header_line_str():
    """Test HeaderLine string representation"""
    line = HeaderLine("test", "value")
    assert str(line) == "HeaderLine('test', 'value')"


def test_header_line_repr():
    """Test HeaderLine repr representation"""
    line = HeaderLine("test", "value")
    assert repr(line) == "HeaderLine('test', 'value')"


def test_meta_header_line_hash_raises_error():
    """Test that MetaHeaderLine.__hash__ raises TypeError"""
    mapping = {"ID": "test", "Description": "Test desc"}
    line = MetaHeaderLine("META", mapping_to_str(mapping), mapping)
    with pytest.raises(TypeError, match="Unhashable type: MetaHeaderLine"):
        hash(line)


def test_meta_header_line_str():
    """Test MetaHeaderLine string representation"""
    mapping = {"ID": "test", "Description": "Test desc"}
    line = MetaHeaderLine("META", mapping_to_str(mapping), mapping)
    expected = f"MetaHeaderLine('META', '{mapping_to_str(mapping)}', {mapping})"
    assert str(line) == expected


def test_meta_header_line_repr():
    """Test MetaHeaderLine repr representation"""
    mapping = {"ID": "test", "Description": "Test desc"}
    line = MetaHeaderLine("META", mapping_to_str(mapping), mapping)
    expected = f"MetaHeaderLine('META', '{mapping_to_str(mapping)}', {mapping})"
    assert repr(line) == expected


def test_samples_infos_hash_raises_error():
    """Test that SamplesInfos.__hash__ raises TypeError"""
    samples = SamplesInfos(["sample1", "sample2"])
    with pytest.raises(TypeError, match="Unhashable type: SamplesInfos"):
        hash(samples)


def test_samples_infos_str():
    """Test SamplesInfos string representation"""
    samples = SamplesInfos(["sample1", "sample2"])
    expected = "SamplesInfos(names=['sample1', 'sample2'], name_to_idx={'sample1': 0, 'sample2': 1})"
    assert str(samples) == expected


def test_samples_infos_repr():
    """Test SamplesInfos repr representation"""
    samples = SamplesInfos(["sample1", "sample2"])
    expected = "SamplesInfos(names=['sample1', 'sample2'], name_to_idx={'sample1': 0, 'sample2': 1})"
    assert repr(samples) == expected


def test_header_invalid_type_exception():
    """Test HeaderInvalidType exception"""
    exc = HeaderInvalidType("Test error message")
    assert str(exc) == "Test error message"
    assert isinstance(exc, Exception)


def test_header_missing_description_exception():
    """Test HeaderMissingDescription exception"""
    exc = HeaderMissingDescription("Missing description")
    assert str(exc) == "Missing description"
    assert isinstance(exc, Exception)


def test_serialize_for_header_function():
    """Test serialize_for_header function with various inputs"""
    # Test with string containing quotes
    result = serialize_for_header("Description", 'Test "quoted" value')
    assert result == '"Test \\"quoted\\" value"'

    # Test with special characters
    result = serialize_for_header("Description", "Test,<>=value")
    assert result == '"Test,<>=value"'


def test_mapping_to_str_function():
    """Test mapping_to_str function"""
    mapping = {"ID": "TEST", "Number": "1", "Type": "String", "Description": "Test field"}
    result = mapping_to_str(mapping)

    # Result should be in the format <key=value,key=value,...>
    assert result.startswith("<") and result.endswith(">")
    assert "ID=TEST" in result
    assert "Number=1" in result
    assert "Type=String" in result
    assert 'Description="Test field"' in result


def test_samples_infos_not_equal():
    """Test SamplesInfos inequality comparison"""
    samples1 = SamplesInfos(["sample1", "sample2"])
    samples2 = SamplesInfos(["sample1", "sample3"])
    result = samples1.__ne__(samples2)
    assert result is True


def test_samples_infos_not_equal_other_type():
    """Test SamplesInfos inequality with different type returns NotImplemented"""
    samples = SamplesInfos(["sample1"])
    result = samples.__ne__("not_samples")
    assert result is NotImplemented


def test_header_build_indices_invalid_line_type():
    """Test Header._build_indices with invalid line type"""

    # Create a mock header line that has mapping but wrong type
    class BadHeaderLine(HeaderLine):
        def __init__(self):
            super().__init__("FORMAT", "test_value")
            self.mapping = {"ID": "test"}

    bad_line = BadHeaderLine()
    header_obj = Header()
    header_obj.lines = [bad_line]

    with pytest.raises(HeaderInvalidType, match="Header line must be of type SimpleHeaderLine or CompoundHeaderLine"):
        header_obj._build_indices()


def test_format_header_line_missing_type():
    """Test FormatHeaderLine with missing Type field"""
    mapping = {"ID": "GT", "Number": "1", "Description": "Genotype"}

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        line = FormatHeaderLine("FORMAT", mapping_to_str(mapping), mapping)

        assert len(w) == 1
        assert issubclass(w[0].category, HeaderInvalidType)
        assert 'Field "Type" not found in header line' in str(w[0].message)
        assert line.type == "String"


def test_format_header_line_invalid_type():
    """Test FormatHeaderLine with invalid Type field"""
    mapping = {"ID": "GT", "Number": "1", "Type": "InvalidType", "Description": "Genotype"}

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        line = FormatHeaderLine("FORMAT", mapping_to_str(mapping), mapping)

        assert len(w) == 1
        assert issubclass(w[0].category, HeaderInvalidType)
        assert "Invalid FORMAT value type InvalidType" in str(w[0].message)
        assert line.type == "String"


def test_format_header_line_missing_description():
    """Test FormatHeaderLine with missing Description field"""
    mapping = {"ID": "GT", "Number": "1", "Type": "String"}

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        FormatHeaderLine("FORMAT", mapping_to_str(mapping), mapping)

        assert len(w) == 1
        assert issubclass(w[0].category, HeaderMissingDescription)
        assert 'Field "Description" not found' in str(w[0].message)


def test_mapping_to_str():
    """Test mapping_to_str function with multiple values"""
    mapping = {"ID": "test", "Number": "1", "Type": "String", "Description": "Test field"}
    result = mapping_to_str(mapping)
    # Order may vary in Python dict iteration, so check components
    assert result.startswith("<")
    assert result.endswith(">")
    assert "ID=test" in result
    assert "Number=1" in result
    assert "Type=String" in result
    assert 'Description="Test field"' in result


def test_serialize_for_header_with_quotes():
    """Test serialize_for_header with strings containing quotes"""
    result = serialize_for_header("Description", 'Test "quoted" value')
    assert result == '"Test \\"quoted\\" value"'


def test_serialize_for_header_with_special_chars():
    """Test serialize_for_header with special characters"""
    result = serialize_for_header("Description", "Test,<>=value")
    assert result == '"Test,<>=value"'


def test_meta_header_line_from_mapping():
    """Test MetaHeaderLine.from_mapping class method"""
    mapping = {"ID": "test_meta", "Description": "Test metadata"}
    line = MetaHeaderLine.from_mapping(mapping)

    assert line.key == "META"
    assert line.mapping == mapping
    assert line.id == "test_meta"
    assert line.value == mapping_to_str(mapping)


def test_header_line_types_coverage():
    """Test various header line types for coverage"""
    # Test PedigreeHeaderLine
    mapping = {"ID": "sample1", "Father": "father1", "Mother": "mother1"}
    pedigree_line = PedigreeHeaderLine("PEDIGREE", mapping_to_str(mapping), mapping)
    assert pedigree_line.id == "sample1"

    # Test SampleHeaderLine
    sample_mapping = {"ID": "sample1", "Genomes": "Germline", "Mixture": "1.0"}
    sample_line = SampleHeaderLine("SAMPLE", mapping_to_str(sample_mapping), sample_mapping)
    assert sample_line.id == "sample1"

    # Test ContigHeaderLine
    contig_mapping = {"ID": "chr1", "length": "249250621"}
    contig_line = ContigHeaderLine("contig", mapping_to_str(contig_mapping), contig_mapping)
    assert contig_line.id == "chr1"

    # Test AltAlleleHeaderLine
    alt_mapping = {"ID": "DEL", "Description": "Deletion"}
    alt_line = AltAlleleHeaderLine("ALT", mapping_to_str(alt_mapping), alt_mapping)
    assert alt_line.id == "DEL"

    # Test FilterHeaderLine
    filter_mapping = {"ID": "LowQual", "Description": "Low quality"}
    filter_line = FilterHeaderLine("FILTER", mapping_to_str(filter_mapping), filter_mapping)
    assert filter_line.id == "LowQual"


def test_samples_infos_edge_cases():
    """Test SamplesInfos edge cases"""
    # Test empty samples
    samples = SamplesInfos([])
    assert len(samples.names) == 0

    # Test with duplicate names - should store last occurrence
    samples = SamplesInfos(["sample1", "sample1", "sample2"])
    assert samples.name_to_idx["sample1"] == 1  # Last occurrence
    assert samples.name_to_idx["sample2"] == 2


def test_field_info_type_validation():
    """Test FieldInfo with edge cases"""
    # Test with missing number field - should work fine
    field_info = FieldInfo("String", None, "Test")
    assert field_info.number is None

    # FieldInfo should be hashable actually
    hash(field_info)  # Should not raise


def test_simple_header_lines():
    """Test SimpleHeaderLine without ID requirement"""
    # Regular HeaderLine should be fine
    line = HeaderLine("fileformat", "VCFv4.2")
    assert line.key == "fileformat"
    assert line.value == "VCFv4.2"


def test_header_line_types_with_mapping_to_str():
    """Test various header line types for coverage"""
    # Test PedigreeHeaderLine
    mapping = {"ID": "sample1", "Father": "father1", "Mother": "mother1"}
    pedigree_line = PedigreeHeaderLine("PEDIGREE", mapping_to_str(mapping), mapping)
    assert pedigree_line.id == "sample1"

    # Test SampleHeaderLine
    sample_mapping = {"ID": "sample1", "Genomes": "Germline", "Mixture": "1.0"}
    sample_line = SampleHeaderLine("SAMPLE", mapping_to_str(sample_mapping), sample_mapping)
    assert sample_line.id == "sample1"

    # Test ContigHeaderLine
    contig_mapping = {"ID": "chr1", "length": "249250621"}
    contig_line = ContigHeaderLine("contig", mapping_to_str(contig_mapping), contig_mapping)
    assert contig_line.id == "chr1"

    # Test AltAlleleHeaderLine
    alt_mapping = {"ID": "DEL", "Description": "Deletion"}
    alt_line = AltAlleleHeaderLine("ALT", mapping_to_str(alt_mapping), alt_mapping)
    assert alt_line.id == "DEL"

    # Test FilterHeaderLine
    filter_mapping = {"ID": "LowQual", "Description": "Low quality"}
    filter_line = FilterHeaderLine("FILTER", mapping_to_str(filter_mapping), filter_mapping)
    assert filter_line.id == "LowQual"


def test_samples_infos_index_behavior():
    """Test SamplesInfos indexing behavior"""
    # Test empty samples
    samples = SamplesInfos([])
    assert len(samples.names) == 0

    # Test with duplicate names - should store last occurrence
    samples = SamplesInfos(["sample1", "sample1", "sample2"])
    assert samples.name_to_idx["sample1"] == 1  # Last occurrence
    assert samples.name_to_idx["sample2"] == 2
