# -*- coding: utf-8 -*-
"""Test the structural variant Record classes
"""

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
