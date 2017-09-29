# -*- coding: utf-8 -*-
"""Test parsing of breakend records
"""

import io

from vcfpy import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


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
##INFO=<ID=ANNO,Number=.,Type=String,Description="Additional annotation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Call-wise filters">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


def vcf_parser(lines):
    return parser.Parser(io.StringIO(MEDIUM_HEADER + lines), '<builtin>')


def test_parse_simple_breakend():
    # Setup parser with stock header and lines to parse
    LINES = '2\t321681\tbnd_W\tG\tG]17:198982]\t6\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n'
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        """Record('2', 321681, ['bnd_W'], 'G', """
        """[BreakEnd('17', 198982, '-', '-', 'G', True)], 6, """
        """['PASS'], OrderedDict([('SVTYPE', 'BND')]), ['GT'], """
        """[Call('NA00001', OrderedDict([('GT', '0/1')])), """
        """Call('NA00002', OrderedDict([('GT', '0/0')])), """
        """Call('NA00003', OrderedDict([('GT', '0/0')]))])""")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert RESULT.ALT[0].serialize() == 'G]17:198982]'


def test_parse_breakend_with_seq():
    # Setup parser with stock header and lines to parse
    LINES = '2\t321681\tbnd_V\tT\t]13:123456]AGTNNNNNCAT\t6\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n'
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        """Record('2', 321681, ['bnd_V'], 'T', """
        """[BreakEnd('13', 123456, '+', '-', 'AGTNNNNNCAT', True)], 6, """
        """['PASS'], OrderedDict([('SVTYPE', 'BND')]), ['GT'], """
        """[Call('NA00001', OrderedDict([('GT', '0/1')])), """
        """Call('NA00002', OrderedDict([('GT', '0/0')])), """
        """Call('NA00003', OrderedDict([('GT', '0/0')]))])""")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert RESULT.ALT[0].serialize() == ']13:123456]AGTNNNNNCAT'


def test_parse_breakend_telomere():
    # Setup parser with stock header and lines to parse
    LINES = '2\t321681\tbnd_V\tN\t.[13:123457[\t6\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n'
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        """Record('2', 321681, ['bnd_V'], 'N', """
        """[BreakEnd('13', 123457, '-', '+', '.', True)], 6, """
        """['PASS'], OrderedDict([('SVTYPE', 'BND')]), ['GT'], """
        """[Call('NA00001', OrderedDict([('GT', '0/1')])), """
        """Call('NA00002', OrderedDict([('GT', '0/0')])), """
        """Call('NA00003', OrderedDict([('GT', '0/0')]))])""")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert RESULT.ALT[0].serialize() == '.[13:123457['


def test_parse_breakend_multi_mate():
    # Setup parser with stock header and lines to parse
    LINES = '2\t321681\tbnd_U\tT\tC[2:321682[,C[17:198983\t6\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n'
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        """Record('2', 321681, ['bnd_U'], 'T', """
        """[BreakEnd('2', 321682, '-', '+', 'C', True), """
        """BreakEnd('17', 198983, '-', '+', 'C', True)], 6, ['PASS'], """
        """OrderedDict([('SVTYPE', 'BND')]), ['GT'], """
        """[Call('NA00001', OrderedDict([('GT', '0/1')])), """
        """Call('NA00002', OrderedDict([('GT', '0/0')])), """
        """Call('NA00003', OrderedDict([('GT', '0/0')]))])""")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert RESULT.ALT[0].serialize() == 'C[2:321682['
    assert RESULT.ALT[1].serialize() == 'C[17:198983['


def test_parse_breakend_single_breakend_fwd():
    # Setup parser with stock header and lines to parse
    LINES = '13\t123457\tbnd_X\tA\t.A\t6\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n'
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        """Record('13', 123457, ['bnd_X'], 'A', [SingleBreakEnd('+', 'A')], """
        """6, ['PASS'], OrderedDict([('SVTYPE', 'BND')]), ['GT'], """
        """[Call('NA00001', OrderedDict([('GT', '0/1')])), """
        """Call('NA00002', OrderedDict([('GT', '0/0')])), """
        """Call('NA00003', OrderedDict([('GT', '0/0')]))])""")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert RESULT.ALT[0].serialize() == '.A'


def test_parse_breakend_single_breakend_rev():
    # Setup parser with stock header and lines to parse
    LINES = '13\t123457\tbnd_X\tA\tA.\t6\tPASS\tSVTYPE=BND\tGT\t0/1\t0/0\t0/0\n'
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    EXPECTED = (
        """Record('13', 123457, ['bnd_X'], 'A', [SingleBreakEnd('-', 'A')], """
        """6, ['PASS'], OrderedDict([('SVTYPE', 'BND')]), ['GT'], """
        """[Call('NA00001', OrderedDict([('GT', '0/1')])), """
        """Call('NA00002', OrderedDict([('GT', '0/0')])), """
        """Call('NA00003', OrderedDict([('GT', '0/0')]))])""")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
    assert RESULT.ALT[0].serialize() == 'A.'
