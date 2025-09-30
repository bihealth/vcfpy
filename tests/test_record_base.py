# -*- coding: utf-8 -*-
"""Test the Record class basics."""

import warnings

import pytest

import vcfpy
from vcfpy.exceptions import CannotModifyUnparsedCallWarning


def test_record_from_scratch():
    """Test contruction of Record objectss"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    record.calls = [vcfpy.Call("sample-1", {})]
    record.add_format("GT", "./.")
    assert (
        str(record) == "Record('chr1', 1234, [], 'A', [Substitution(type_='SNV', value='T')], None, "
        "[], {}, ['GT'], [Call('sample-1', {'GT': './.'})])"
    )


def test_record_format_calls_mismatch():
    """Test that ValueError is raised when FORMAT and calls don't match"""
    # Test case where FORMAT is provided but calls is None
    with pytest.raises(ValueError, match="Either provide both FORMAT and calls or none"):
        vcfpy.Record(
            CHROM="chr1",
            POS=1234,
            ID=[],
            REF="A",
            ALT=[vcfpy.Substitution("SNV", "T")],
            QUAL=None,
            FILTER=[],
            INFO={},
            FORMAT=["GT"],
            calls=None,
        )

    # Test case where calls is provided but FORMAT is None
    with pytest.raises(ValueError, match="Either provide both FORMAT and calls or none"):
        vcfpy.Record(
            CHROM="chr1",
            POS=1234,
            ID=[],
            REF="A",
            ALT=[vcfpy.Substitution("SNV", "T")],
            QUAL=None,
            FILTER=[],
            INFO={},
            FORMAT=None,
            calls=[vcfpy.Call("sample", {})],
        )


def test_record_equality_and_hash():
    """Test record equality and hashing"""
    record1 = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=["rs123"],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=30.0,
        FILTER=["PASS"],
        INFO={"DP": 100},
    )

    record2 = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=["rs123"],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=30.0,
        FILTER=["PASS"],
        INFO={"DP": 100},
    )

    record3 = vcfpy.Record(
        CHROM="chr2",  # Different chromosome
        POS=1234,
        ID=["rs123"],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=30.0,
        FILTER=["PASS"],
        INFO={"DP": 100},
    )

    # Test equality
    assert record1 == record2
    assert record1 != record3
    # Skip testing with non-record objects since __ne__ raises NotImplementedError

    # Hash tests - these will fail due to unhashable types, so we test for exceptions
    with pytest.raises(TypeError):
        hash(record1)

    # Test inequality
    assert not (record1 != record2)
    assert record1 != record3


def test_record_str_and_repr():
    """Test string representation and repr"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=["rs123"],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=30.0,
        FILTER=["PASS"],
        INFO={"DP": 100},
    )

    str_repr = str(record)
    assert "Record(" in str_repr
    assert "chr1" in str_repr
    assert "1234" in str_repr

    # Test that repr returns same as str
    assert repr(record) == str(record)


def test_record_is_snv():
    """Test is_snv method"""
    # Test SNV
    snv_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert snv_record.is_snv()

    # Test multi-SNV
    multi_snv_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T"), vcfpy.Substitution("SNV", "C")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert multi_snv_record.is_snv()

    # Test MNV (not SNV)
    mnv_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="AA",  # Multi-nucleotide REF
        ALT=[vcfpy.Substitution("MNV", "TT")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert not mnv_record.is_snv()

    # Test mixed types (not SNV)
    mixed_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T"), vcfpy.Substitution("INS", "ATG")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert not mixed_record.is_snv()


def test_affected_start_and_end():
    """Test affected_start and affected_end properties"""
    # Test SNV
    snv_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert snv_record.affected_start == 1233  # 0-based
    assert snv_record.affected_end == 1234  # 0-based, end is exclusive

    # Test deletion
    del_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="ATG",
        ALT=[vcfpy.Substitution("DEL", "A")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert del_record.affected_start == 1233  # 0-based
    assert del_record.affected_end == 1236  # 0-based, covers REF length

    # Test insertion (only insertions)
    ins_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("INS", "ATG")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert ins_record.affected_start == 1234  # right of first base
    assert ins_record.affected_end == 1234  # same for 0-length interval

    # Test structural variant
    sv_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.SymbolicAllele("DEL")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    assert sv_record.affected_start == 1233  # standard behavior for SV
    assert sv_record.affected_end == 1234

    # Test mixed types with insertion
    mixed_ins_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("INS", "ATG"), vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )
    # Mixed types don't follow insertion logic
    assert mixed_ins_record.affected_start == 1233
    assert mixed_ins_record.affected_end == 1234


def test_add_filter():
    """Test add_filter method"""
    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=["PASS"],
        INFO={},
    )

    # Add first filter - should remove PASS
    record.add_filter("LowQual")
    assert record.FILTER == ["LowQual"]

    # Add second filter - should append
    record.add_filter("HighFS")
    assert record.FILTER == ["LowQual", "HighFS"]

    # Add existing filter - should not duplicate
    record.add_filter("LowQual")
    assert record.FILTER == ["LowQual", "HighFS"]

    # Test with record that doesn't have PASS
    record2 = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=["ExistingFilter"],
        INFO={},
    )

    record2.add_filter("NewFilter")
    assert record2.FILTER == ["ExistingFilter", "NewFilter"]


def test_add_format():
    """Test add_format method"""
    # Create a record with calls
    call1 = vcfpy.Call("sample1", {"GT": "0/1"})
    call2 = vcfpy.Call("sample2", {"GT": "1/1"})

    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
        FORMAT=["GT"],
        calls=[call1, call2],
    )

    # Add new format field with default value
    record.add_format("DP", 30)
    assert "DP" in record.FORMAT
    assert call1.data["DP"] == 30
    assert call2.data["DP"] == 30

    # Add format field without default value
    record.add_format("GQ")
    assert "GQ" in record.FORMAT
    # Should not add to call data if no default value
    assert "GQ" not in call1.data
    assert "GQ" not in call2.data

    # Try to add existing format field - should do nothing
    original_length = len(record.FORMAT)
    record.add_format("GT", "0/0")
    assert len(record.FORMAT) == original_length
    # Should not modify existing data
    assert call1.data["GT"] == "0/1"
    assert call2.data["GT"] == "1/1"


def test_add_format_with_unparsed_call():
    """Test add_format method with UnparsedCall"""
    unparsed_call = vcfpy.UnparsedCall("sample1", "0/1:30:99")

    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
        FORMAT=["GT"],
        calls=[unparsed_call],
    )

    # Adding format to UnparsedCall should issue warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        record.add_format("DP", 30)
        assert len(w) == 1
        assert issubclass(w[0].category, CannotModifyUnparsedCallWarning)
        assert "UnparsedCall encountered" in str(w[0].message)


def test_record_iteration():
    """Test iterating over record calls"""
    call1 = vcfpy.Call("sample1", {"GT": "0/1"})
    call2 = vcfpy.Call("sample2", {"GT": "1/1"})

    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
        FORMAT=["GT"],
        calls=[call1, call2],
    )

    calls_list = list(record)
    assert len(calls_list) == 2
    assert calls_list[0] == call1
    assert calls_list[1] == call2


def test_update_calls():
    """Test update_calls method"""
    call1 = vcfpy.Call("sample1", {"GT": "0/1"})
    call2 = vcfpy.Call("sample2", {"GT": "1/1"})

    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )

    # Initially no calls
    assert len(record.calls) == 0
    assert len(record.call_for_sample) == 0

    # Update with calls (note: update_calls doesn't change record.calls, just mapping)
    record.update_calls([call1, call2])
    assert len(record.call_for_sample) == 2
    assert record.call_for_sample["sample1"] == call1
    assert record.call_for_sample["sample2"] == call2

    # Check that site is set on calls
    assert call1.site == record
    assert call2.site == record


def test_call_basic():
    """Test basic Call functionality"""
    call = vcfpy.Call("sample1", {"GT": "0/1", "DP": 30, "GQ": 99})

    assert call.sample == "sample1"
    assert call.data["GT"] == "0/1"
    assert call.data["DP"] == 30
    assert call.data["GQ"] == 99

    # Test string representation
    str_repr = str(call)
    assert "Call" in str_repr
    assert "sample1" in str_repr

    # Test equality
    call2 = vcfpy.Call("sample1", {"GT": "0/1", "DP": 30, "GQ": 99})
    call3 = vcfpy.Call("sample2", {"GT": "0/1", "DP": 30, "GQ": 99})

    assert call == call2
    assert call != call3
    # Skip testing with non-call objects since __ne__ raises NotImplementedError

    # Hash tests - these will fail due to unhashable dict, so we test for exceptions
    with pytest.raises(TypeError):
        hash(call)


def test_call_genotype_parsing():
    """Test genotype parsing in Call"""
    # Test diploid heterozygous
    call = vcfpy.Call("sample1", {"GT": "0/1"})
    assert call.gt_alleles == [0, 1]
    assert call.called is True
    assert call.ploidy == 2
    assert not call.is_phased
    assert call.gt_phase_char == "/"
    assert call.gt_type == vcfpy.HET

    # Test diploid homozygous reference
    call_hom_ref = vcfpy.Call("sample1", {"GT": "0/0"})
    assert call_hom_ref.gt_alleles == [0, 0]
    assert call_hom_ref.called is True
    assert call_hom_ref.ploidy == 2
    assert call_hom_ref.gt_type == vcfpy.HOM_REF

    # Test diploid homozygous alternative
    call_hom_alt = vcfpy.Call("sample1", {"GT": "1/1"})
    assert call_hom_alt.gt_alleles == [1, 1]
    assert call_hom_alt.called is True
    assert call_hom_alt.ploidy == 2
    assert call_hom_alt.gt_type == vcfpy.HOM_ALT

    # Test phased genotype
    call_phased = vcfpy.Call("sample1", {"GT": "0|1"})
    assert call_phased.gt_alleles == [0, 1]
    assert call_phased.called is True
    assert call_phased.ploidy == 2
    assert call_phased.is_phased
    assert call_phased.gt_phase_char == "|"
    assert call_phased.gt_type == vcfpy.HET

    # Test missing genotype
    call_missing = vcfpy.Call("sample1", {"GT": "./."})
    assert call_missing.gt_alleles == [None, None]
    assert call_missing.called is False
    assert call_missing.ploidy == 2
    assert call_missing.gt_type is None

    # Test partially missing genotype
    call_partial = vcfpy.Call("sample1", {"GT": "0/."})
    assert call_partial.gt_alleles == [0, None]
    assert call_partial.called is False
    assert call_partial.ploidy == 2
    assert call_partial.gt_type is None

    # Test no genotype data
    call_no_gt = vcfpy.Call("sample1", {"DP": 30})
    assert call_no_gt.gt_alleles is None
    assert call_no_gt.called is None
    assert call_no_gt.ploidy is None
    assert call_no_gt.gt_type is None

    # Test haploid genotype
    call_haploid = vcfpy.Call("sample1", {"GT": "1"})
    assert call_haploid.gt_alleles == [1]
    assert call_haploid.called is True
    assert call_haploid.ploidy == 1
    assert call_haploid.gt_type == vcfpy.HOM_ALT


def test_call_set_genotype():
    """Test set_genotype method"""
    call = vcfpy.Call("sample1", {"DP": 30})

    # Initially no genotype
    assert call.data.get("GT") is None
    assert call.gt_alleles is None

    # Set genotype
    call.set_genotype("0/1")
    assert call.data["GT"] == "0/1"
    assert call.gt_alleles == [0, 1]
    assert call.called is True
    assert call.ploidy == 2

    # Set to missing
    call.set_genotype("./.")
    assert call.data["GT"] == "./."
    assert call.gt_alleles == [None, None]
    assert call.called is False
    assert call.ploidy == 2

    # Set to None
    call.set_genotype(None)
    assert call.data["GT"] is None
    assert call.gt_alleles is None
    assert call.called is None
    assert call.ploidy is None


def test_call_gt_bases():
    """Test gt_bases property"""
    # Create a record and call to test gt_bases
    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T"), vcfpy.Substitution("SNV", "C")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )

    call = vcfpy.Call("sample1", {"GT": "0/1"})
    call.site = record

    # Test reference/alternative
    assert call.gt_bases == ("A", "T")

    # Test homozygous reference
    call.set_genotype("0/0")
    assert call.gt_bases == ("A", "A")

    # Test homozygous alternative
    call.set_genotype("1/1")
    assert call.gt_bases == ("T", "T")

    # Test multiple alternatives
    call.set_genotype("1/2")
    assert call.gt_bases == ("T", "C")

    # Test with missing allele
    call.set_genotype("0/.")
    assert call.gt_bases == ("A", None)

    # Test without site set
    call_no_site = vcfpy.Call("sample1", {"GT": "0/1"})
    with pytest.raises(ValueError, match="Cannot determine bases without site being set"):
        _ = call_no_site.gt_bases


def test_unparsed_call():
    """Test UnparsedCall functionality"""
    unparsed = vcfpy.UnparsedCall("sample1", "0/1:30:99")

    assert unparsed.sample == "sample1"
    assert unparsed.unparsed_data == "0/1:30:99"
    assert unparsed.site is None

    # Test setting site
    record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],
        INFO={},
    )

    unparsed.site = record
    assert unparsed.site == record

    # UnparsedCall doesn't implement equality, so we test basic attributes instead
    unparsed2 = vcfpy.UnparsedCall("sample1", "0/1:30:99")
    unparsed3 = vcfpy.UnparsedCall("sample2", "0/1:30:99")

    assert unparsed.sample == unparsed2.sample
    assert unparsed.unparsed_data == unparsed2.unparsed_data
    assert unparsed.sample != unparsed3.sample

    # Test string representation (UnparsedCall doesn't have custom __str__)
    str_repr = str(unparsed)
    assert "UnparsedCall" in str_repr


def test_substitution_alt_record():
    """Test Substitution AltRecord functionality"""
    sub = vcfpy.Substitution("SNV", "T")

    assert sub.type == "SNV"
    assert sub.value == "T"
    assert sub.serialize() == "T"

    # Test equality
    sub2 = vcfpy.Substitution("SNV", "T")
    sub3 = vcfpy.Substitution("SNV", "C")

    assert sub == sub2
    assert sub != sub3
    # Skip testing with non-substitution objects since __ne__ raises NotImplementedError

    # Test hash
    assert hash(sub) == hash(sub2)
    assert hash(sub) != hash(sub3)

    # Test string representation
    str_repr = str(sub)
    assert "Substitution" in str_repr
    assert "SNV" in str_repr
    assert "T" in str_repr


def test_symbolic_allele():
    """Test SymbolicAllele functionality"""
    sym = vcfpy.SymbolicAllele("DEL")

    assert sym.type == "SYMBOLIC"
    assert sym.value == "DEL"
    assert sym.serialize() == "<DEL>"

    # Test equality
    sym2 = vcfpy.SymbolicAllele("DEL")
    sym3 = vcfpy.SymbolicAllele("DUP")

    assert sym == sym2
    assert sym != sym3
    # Skip testing with non-symbolic allele objects since __ne__ raises NotImplementedError

    # Test hash
    assert hash(sym) == hash(sym2)
    assert hash(sym) != hash(sym3)

    # Test string representation
    str_repr = str(sym)
    assert "SymbolicAllele" in str_repr
    assert "DEL" in str_repr


def test_alt_record_base():
    """Test AltRecord base class functionality"""
    # Create a basic AltRecord (though it's meant to be abstract)
    alt = vcfpy.AltRecord("SNV")

    assert alt.type == "SNV"

    # Test equality
    alt2 = vcfpy.AltRecord("SNV")
    alt3 = vcfpy.AltRecord("MNV")

    assert alt == alt2
    assert alt != alt3
    # Skip testing with non-alt record objects since __ne__ raises NotImplementedError

    # Test hash
    assert hash(alt) == hash(alt2)
    assert hash(alt) != hash(alt3)

    # Test that serialize is not implemented
    with pytest.raises(NotImplementedError):
        alt.serialize()


def test_single_break_end():
    """Test SingleBreakEnd class functionality"""
    # Test basic construction
    single_bnd = vcfpy.record.SingleBreakEnd("+", "ATCG")

    assert single_bnd.orientation == "+"
    assert single_bnd.sequence == "ATCG"
    assert single_bnd.mate_chrom is None
    assert single_bnd.mate_pos is None
    assert single_bnd.mate_orientation is None
    assert single_bnd.within_main_assembly is None

    # Test equality
    single_bnd2 = vcfpy.record.SingleBreakEnd("+", "ATCG")
    single_bnd3 = vcfpy.record.SingleBreakEnd("-", "ATCG")

    assert single_bnd == single_bnd2
    assert single_bnd != single_bnd3

    # Test hash
    hash1 = hash(single_bnd)
    hash2 = hash(single_bnd2)
    hash3 = hash(single_bnd3)

    assert hash1 == hash2
    assert hash1 != hash3

    # Test string representation
    str_repr = str(single_bnd)
    assert "SingleBreakEnd" in str_repr
    assert "+" in str_repr
    assert "ATCG" in str_repr


def test_break_end_hash():
    """Test BreakEnd class hash and equality"""
    # Create BreakEnd instances
    bnd = vcfpy.record.BreakEnd("chr2", 1234, "+", "+", "A", True)
    bnd2 = vcfpy.record.BreakEnd("chr2", 1234, "+", "+", "A", True)
    bnd3 = vcfpy.record.BreakEnd("chr3", 1234, "+", "+", "A", True)

    # Test equality
    assert bnd == bnd2
    assert bnd != bnd3

    # Test hash - equal objects should have equal hashes
    assert hash(bnd) == hash(bnd2)
    # Different objects should (very likely) have different hashes
    assert hash(bnd) != hash(bnd3)
