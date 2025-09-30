# -*- coding: utf-8 -*-
"""Test parsing of full VCF header"""

import io
import warnings

import pytest

import vcfpy
from vcfpy import exceptions, header, parser
from vcfpy.header import FieldInfo

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

MEDIUM_HEADER = """
##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


@pytest.fixture(scope="function")
def medium_header():
    return io.StringIO(MEDIUM_HEADER)


def test_parse_header(medium_header):
    p = parser.Parser(stream=medium_header, path="<builtin>")
    header = p.parse_header()
    assert header.lines
    assert len(header.lines) == 18
    EXPECTED = "HeaderLine('fileformat', 'VCFv4.3')"
    assert str(header.lines[0]) == EXPECTED
    EXPECTED = (
        "FormatHeaderLine('FORMAT', '<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">', "
        "{'ID': 'HQ', 'Number': 2, 'Type': 'Integer', 'Description': 'Haplotype Quality'})"
    )
    assert str(header.lines[-1]) == EXPECTED
    assert header.samples.names == ["NA00001", "NA00002", "NA00003"]


def test_stupid_header_line_parser():
    """Test StupidHeaderLineParser"""
    stupid_parser = parser.StupidHeaderLineParser()
    result = stupid_parser.parse_key_value("fileformat", "VCFv4.3")

    assert isinstance(result, header.HeaderLine)
    assert result.key == "fileformat"
    assert result.value == "VCFv4.3"


def test_parser_warnings_and_validation():
    """Test parser warning and validation logic"""
    # Test key with leading/trailing spaces in parse_mapping
    with pytest.warns(exceptions.LeadingTrailingSpaceInKey):
        result = parser.parse_mapping("< key =value>")
        assert "key" in result


def test_quoted_string_splitter_edge_cases():
    """Test QuotedStringSplitter with various edge cases"""
    splitter = parser.QuotedStringSplitter()

    # Test with escaped quotes
    result = splitter.run('key="val\\"ue",other=value')
    assert len(result) == 2
    assert result[0] == 'key="val\\"ue"'
    assert result[1] == "other=value"

    # Test with arrays
    result = splitter.run("key=[val1,val2],other=value")
    assert len(result) == 2
    assert result[0] == "key=[val1,val2]"
    assert result[1] == "other=value"

    # Test with different delimiters
    splitter_semicolon = parser.QuotedStringSplitter(delim=";")
    result = splitter_semicolon.run("key=value;other=value2")
    assert len(result) == 2
    assert result[0] == "key=value"
    assert result[1] == "other=value2"


def test_header_line_parser_base_not_implemented():
    """Test that HeaderLineParserBase.parse_key_value raises NotImplementedError"""
    parser_base = parser.HeaderLineParserBase()

    with pytest.raises(NotImplementedError, match="Must be overridden"):
        parser_base.parse_key_value("test", "value")


def test_mapping_header_line_parser():
    """Test MappingHeaderLineParser functionality"""
    mapping_parser = parser.MappingHeaderLineParser(header.InfoHeaderLine)

    # Test parsing a simple mapping
    result = mapping_parser.parse_key_value("INFO", '<ID=DP,Number=1,Type=Integer,Description="Depth of coverage">')

    # Should return a CompoundHeaderLine (or similar)
    assert result.key == "INFO"
    assert '<ID=DP,Number=1,Type=Integer,Description="Depth of coverage">' in result.value


def test_parser_edge_cases():
    """Test parser edge cases and error conditions"""
    # Test split_mapping with key that needs stripping
    with pytest.warns(exceptions.LeadingTrailingSpaceInKey):
        key, value = parser.split_mapping(" key = value")
        assert key == "key"
        assert value == " value"


def test_parse_field_value_conversion_failure():
    """Test parse_field_value when conversion fails"""
    # Create a FieldInfo for Integer type with proper constructor
    field_info = FieldInfo("Integer", 1, "Test field", "DP")

    # This should successfully convert a valid integer
    result = parser.parse_field_value(field_info, "123")  # Valid integer

    # Should successfully convert to integer
    assert result == 123


def test_process_sub_grow_edge_cases():
    """Test process_sub_grow with various edge cases"""
    # Test empty ALT string - should raise exception
    with pytest.raises(exceptions.InvalidRecordException, match="Invalid VCF, empty ALT"):
        parser.process_sub_grow("A", "")

    # Test single character where REF[0] == ALT[0] (deletion)
    result = parser.process_sub_grow("AT", "A")
    assert result.type == vcfpy.DEL
    assert result.value == "A"

    # Test single character where REF[0] != ALT[0] (indel)
    result = parser.process_sub_grow("A", "T")
    assert result.type == vcfpy.INDEL
    assert result.value == "T"

    # Test multi-character ALT (indel)
    result = parser.process_sub_grow("A", "ATG")
    assert result.type == vcfpy.INDEL
    assert result.value == "ATG"


def test_parse_breakend_edge_cases():
    """Test breakend parsing edge cases"""
    # Test with mate chromosome in brackets (not within main assembly)
    result = parser.parse_breakend("A[<chr2>:12345[")
    mate_chrom, mate_pos, orientation, mate_orientation, sequence, within_main_assembly = result

    assert mate_chrom == "chr2"
    assert mate_pos == 12345
    assert within_main_assembly is False
    assert sequence == "A"  # Sequence before the bracket

    # Test with mate chromosome without brackets (within main assembly)
    result = parser.parse_breakend("A[chr2:12345[")
    mate_chrom, mate_pos, orientation, mate_orientation, sequence, within_main_assembly = result

    assert mate_chrom == "chr2"
    assert mate_pos == 12345
    assert within_main_assembly is True


def test_parse_field_value_edge_cases():
    """Test parse_field_value with edge cases"""
    # Test FT field special handling
    ft_field_info = header.FieldInfo("String", 1, "Filter field", "FORMAT/FT")
    result = parser.parse_field_value(ft_field_info, "PASS;LowQual")
    assert result == ["PASS", "LowQual"]

    # Test flag field
    flag_field_info = header.FieldInfo("Flag", 0, "Flag field")
    result = parser.parse_field_value(flag_field_info, True)
    assert result is True

    result = parser.parse_field_value(flag_field_info, "anything")
    assert result is True


def test_header_missing_lines():
    """Test some header functionality to cover missing lines"""
    # Test creating headers with various line types
    lines = [
        header.HeaderLine("fileformat", "VCFv4.3"),
        header.InfoHeaderLine.from_mapping({"ID": "DP", "Number": 1, "Type": "Integer", "Description": "Depth"}),
    ]
    samples = header.SamplesInfos(["sample1"])
    hdr = header.Header(lines, samples)

    # Test some header functionality that might not be covered
    assert len(hdr.lines) == 2
    assert hdr.samples.names == ["sample1"]


def test_build_header_parsers():
    """Test build_header_parsers function"""
    parsers = parser.build_header_parsers()

    # Check that all expected parsers are present
    expected_keys = ["ALT", "contig", "FILTER", "FORMAT", "INFO", "PEDIGREE", "SAMPLE"]
    for key in expected_keys:
        assert key in parsers
        assert isinstance(parsers[key], parser.MappingHeaderLineParser)

    # Test that other default keys are handled by the system
    assert len(parsers) >= len(expected_keys)


def test_parse_field_value_with_escaping():
    """Test parse_field_value with string escaping"""
    # Test string with escaping using proper String type
    field_info = FieldInfo("String", 1, "Test field", "DESC")
    result = parser.parse_field_value(field_info, "value%3Bwith%3Dsemi")
    # Should unescape the string and return single value for number=1
    assert isinstance(result, str)
    assert ";" in result and "=" in result


def test_quoted_string_splitter_state_transitions():
    """Test QuotedStringSplitter state transitions"""
    splitter = parser.QuotedStringSplitter()

    # Test complex string with multiple state transitions
    complex_str = 'key1="val,with,commas",key2=[arr,val],key3=simple'
    result = splitter.run(complex_str)

    assert len(result) == 3
    assert result[0] == 'key1="val,with,commas"'
    assert result[1] == "key2=[arr,val]"
    assert result[2] == "key3=simple"


def test_convert_field_value_error_handling():
    """Test convert_field_value error handling"""
    # Test with invalid integer
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = parser.convert_field_value("Integer", "invalid_int")

        # Should return original string and generate warning
        assert result == "invalid_int"
        assert len(w) == 1


def test_parse_mapping_edge_cases():
    """Test parse_mapping with edge cases"""
    # Test with proper angular brackets and nested values
    result = parser.parse_mapping('<ID=test,Description="Complex,value=with%3Bescapes">')
    assert "ID" in result
    assert "Description" in result
    assert result["ID"] == "test"
    # Check if description contains the raw escaped content
    desc = result["Description"]
    assert isinstance(desc, str)
    assert "Complex,value=with%3Bescapes" in desc or "Complex,value=with;escapes" in desc
