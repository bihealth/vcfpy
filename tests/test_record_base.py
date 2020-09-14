# -*- coding: utf-8 -*-
"""Test the Record class basics."""

import sys

import vcfpy


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
        INFO=vcfpy.OrderedDict(),
    )
    record.calls = [vcfpy.Call("sample-1", vcfpy.OrderedDict())]
    record.add_format("GT", "./.")
    if sys.version_info < (3, 6):
        assert str(record) == (
            "Record('chr1', 1234, [], 'A', [Substitution(type_='SNV', value='T')], None, "
            "[], OrderedDict(), ['GT'], [Call('sample-1', OrderedDict([('GT', './.')]))])"
        )
    else:
        assert str(record) == (
            "Record('chr1', 1234, [], 'A', [Substitution(type_='SNV', value='T')], None, "
            "[], {}, ['GT'], [Call('sample-1', {'GT': './.'})])"
        )
