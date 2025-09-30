# -*- coding: utf-8 -*-
"""Test the structural variant Record classes"""

import pytest

import vcfpy
from vcfpy import record

# BreakEnd ------------------------------------------------------------------


def test_breakend_builtins():
    """Test the builtin functions of the BreakEnd class"""
    rec1 = record.BreakEnd("chr2", 1234, record.FORWARD, record.FORWARD, "A", True)
    rec2 = record.BreakEnd("chr2", 1234, record.FORWARD, record.FORWARD, "A", True)
    rec3 = record.BreakEnd("chr1", 1234, record.FORWARD, record.FORWARD, "A", True)
    assert rec1 == rec2
    assert hash(rec1) == hash(rec2)
    assert rec1 != rec3
    assert str(rec1) == "BreakEnd('chr2', 1234, '+', '+', 'A', True)"
    assert repr(rec1) == "BreakEnd('chr2', 1234, '+', '+', 'A', True)"


def test_breakend_fwd_fwd_true():
    rec = record.BreakEnd("chr2", 1234, record.FORWARD, record.FORWARD, "A", True)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '+', '+', 'A', True)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "[chr2:1234[A"
    assert EXPECTED == rec.serialize()


def test_breakend_fwd_fwd_false():
    rec = record.BreakEnd("chr2", 1234, record.FORWARD, record.FORWARD, "A", False)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '+', '+', 'A', False)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "[<chr2>:1234[A"
    assert EXPECTED == rec.serialize()


def test_breakend_fwd_rev_true():
    rec = record.BreakEnd("chr2", 1234, record.FORWARD, record.REVERSE, "A", True)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '+', '-', 'A', True)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "]chr2:1234]A"
    assert EXPECTED == rec.serialize()


def test_breakend_fwd_rev_false():
    rec = record.BreakEnd("chr2", 1234, record.FORWARD, record.REVERSE, "A", False)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '+', '-', 'A', False)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "]<chr2>:1234]A"
    assert EXPECTED == rec.serialize()


def test_breakend_rev_fwd_true():
    rec = record.BreakEnd("chr2", 1234, record.REVERSE, record.FORWARD, "A", True)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '-', '+', 'A', True)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "A[chr2:1234["
    assert EXPECTED == rec.serialize()


def test_breakend_rev_fwd_false():
    rec = record.BreakEnd("chr2", 1234, record.REVERSE, record.FORWARD, "A", False)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '-', '+', 'A', False)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "A[<chr2>:1234["
    assert EXPECTED == rec.serialize()


def test_breakend_rev_rev_true():
    rec = record.BreakEnd("chr2", 1234, record.REVERSE, record.REVERSE, "A", True)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '-', '-', 'A', True)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "A]chr2:1234]"
    assert EXPECTED == rec.serialize()


def test_breakend_rev_rev_false():
    rec = record.BreakEnd("chr2", 1234, record.REVERSE, record.REVERSE, "A", False)
    # test initialization
    EXPECTED = "BreakEnd('chr2', 1234, '-', '-', 'A', False)"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "A]<chr2>:1234]"
    assert EXPECTED == rec.serialize()


# SingleBreakEnd ------------------------------------------------------------


def test_single_breakend_fdw():
    rec = record.SingleBreakEnd(record.FORWARD, "A")
    # test initialization
    EXPECTED = "SingleBreakEnd('+', 'A')"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = ".A"
    assert EXPECTED == rec.serialize()


def test_single_breakend_rev():
    rec = record.SingleBreakEnd(record.REVERSE, "A")
    # test initialization
    EXPECTED = "SingleBreakEnd('-', 'A')"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "A."
    assert EXPECTED == rec.serialize()


# SingleBreakEnd ------------------------------------------------------------


def test_symbolic_allele():
    rec = record.SymbolicAllele("DUP")
    # test initialization
    EXPECTED = "SymbolicAllele('DUP')"
    assert EXPECTED == str(rec)
    # test serialize()
    EXPECTED = "<DUP>"
    assert EXPECTED == rec.serialize()


def test_single_break_end():
    """Test SingleBreakEnd functionality"""
    single_bnd = record.SingleBreakEnd("+", "ATCG")

    assert single_bnd.orientation == "+"
    assert single_bnd.sequence == "ATCG"
    assert single_bnd.mate_chrom is None
    assert single_bnd.mate_pos is None
    assert single_bnd.mate_orientation is None
    assert single_bnd.within_main_assembly is None

    # Test equality
    single_bnd2 = record.SingleBreakEnd("+", "ATCG")
    single_bnd3 = record.SingleBreakEnd("-", "ATCG")

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
    """Test BreakEnd hash functionality"""
    bnd = record.BreakEnd("chr2", 1234, "+", "+", "A", True)
    bnd2 = record.BreakEnd("chr2", 1234, "+", "+", "A", True)
    bnd3 = record.BreakEnd("chr3", 1234, "+", "+", "A", True)

    # Test hash
    hash1 = hash(bnd)
    hash2 = hash(bnd2)
    hash3 = hash(bnd3)

    assert hash1 == hash2
    assert hash1 != hash3


def test_record_with_breakends():
    """Test record with breakend alternatives for affected positions"""
    # Test breakend records
    bnd_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[record.BreakEnd("chr2", 5678, "+", "+", "A", True)],
        QUAL=None,
        FILTER=[],
        INFO={},
    )

    # Breakends should use standard affected positions
    assert bnd_record.affected_start == 1233  # 0-based
    assert bnd_record.affected_end == 1234


def test_symbolic_allele_hash():
    """Test SymbolicAllele hash functionality"""
    sym1 = record.SymbolicAllele("DEL")
    sym2 = record.SymbolicAllele("DEL")
    sym3 = record.SymbolicAllele("DUP")

    hash1 = hash(sym1)
    hash2 = hash(sym2)
    hash3 = hash(sym3)

    assert hash1 == hash2
    assert hash1 != hash3


def test_substitution_hash():
    """Test Substitution hash functionality"""
    sub1 = record.Substitution("SNV", "T")
    sub2 = record.Substitution("SNV", "T")
    sub3 = record.Substitution("SNV", "C")

    hash1 = hash(sub1)
    hash2 = hash(sub2)
    hash3 = hash(sub3)

    assert hash1 == hash2
    assert hash1 != hash3


def test_call_gt_type_edge_cases():
    """Test Call gt_type with various edge cases"""
    # Test triploid genotype
    call_triploid = vcfpy.Call("sample1", {"GT": "0/1/2"})
    assert call_triploid.gt_alleles == [0, 1, 2]
    assert call_triploid.called is True
    assert call_triploid.ploidy == 3
    assert call_triploid.gt_type == vcfpy.HET  # Mixed alleles = HET

    # Test homozygous alternative triploid
    call_hom_alt_3 = vcfpy.Call("sample1", {"GT": "1/1/1"})
    assert call_hom_alt_3.gt_alleles == [1, 1, 1]
    assert call_hom_alt_3.called is True
    assert call_hom_alt_3.ploidy == 3
    assert call_hom_alt_3.gt_type == vcfpy.HOM_ALT  # All same non-ref = HOM_ALT


def test_call_hash():
    """Test Call hash functionality"""
    call1 = vcfpy.Call("sample1", {"GT": "0/1"})
    vcfpy.Call("sample1", {"GT": "0/1"})
    vcfpy.Call("sample2", {"GT": "0/1"})

    # Hash should fail due to unhashable dict
    with pytest.raises(TypeError):
        hash(call1)


def test_record_mixed_insertion_types():
    """Test record with mixed insertion and other types"""
    # Mixed types with insertion should not use insertion logic
    mixed_record = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[
            record.Substitution("INS", "ATG"),
            record.SymbolicAllele("DEL"),
        ],
        QUAL=None,
        FILTER=[],
        INFO={},
    )

    # Should not use insertion-specific logic due to mixed types
    assert mixed_record.affected_start == 1233
    assert mixed_record.affected_end == 1234


def test_add_filter_without_pass():
    """Test add_filter when no PASS filter exists"""
    record_no_pass = vcfpy.Record(
        CHROM="chr1",
        POS=1234,
        ID=[],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=[],  # No PASS filter
        INFO={},
    )

    record_no_pass.add_filter("LowQual")
    assert record_no_pass.FILTER == ["LowQual"]

    # Add another filter
    record_no_pass.add_filter("HighFS")
    assert record_no_pass.FILTER == ["LowQual", "HighFS"]


def test_call_gt_bases_error_cases():
    """Test Call gt_bases error conditions"""
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

    # Test with allele index beyond available ALT alleles
    call = vcfpy.Call("sample1", {"GT": "0/3"})  # Only 1 ALT allele available (index 1)
    call.site = record

    # This should raise an IndexError for out-of-range allele index
    with pytest.raises(IndexError):
        _ = call.gt_bases
