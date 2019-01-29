# -*- coding: utf-8 -*-
"""Test parsing of symbolic records
"""

import io
import sys

from vcfpy import parser

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
##INFO=<ID=ANNO,Number=.,Type=String,Description="Additional annotation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="SV end position">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Call-wise filters">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=R,Description="IUPAC code R = A/G">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


def vcf_parser(lines):
    return parser.Parser(io.StringIO(MEDIUM_HEADER + lines), "<builtin>")


def test_parse_dup():
    # Setup parser with stock header and lines to parse
    LINES = "2\t321681\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=324681;SVLEN=3000\tGT\t0/1\t0/0\t0/0\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    if sys.version_info < (3, 6):
        EXPECTED = (
            "Record('2', 321681, [], 'N', [SymbolicAllele('DUP')], None, ['PASS'], "
            "OrderedDict([('SVTYPE', 'DUP'), ('END', 324681), ('SVLEN', 3000)]), ['GT'], ["
            "Call('NA00001', OrderedDict([('GT', '0/1')])), Call('NA00002', OrderedDict([('GT', '0/0')])), "
            "Call('NA00003', OrderedDict([('GT', '0/0')]))])"
        )
    else:
        EXPECTED = (
            "Record('2', 321681, [], 'N', [SymbolicAllele('DUP')], None, ['PASS'], "
            "{'SVTYPE': 'DUP', 'END': 324681, 'SVLEN': 3000}, ['GT'], ["
            "Call('NA00001', {'GT': '0/1'}), Call('NA00002', {'GT': '0/0'}), Call('NA00003', {'GT': '0/0'})])"
        )
    rec = p.parse_next_record()
    assert str(rec) == EXPECTED
    assert rec.ALT[0].serialize() == "<DUP>"


def test_parse_iupac():
    # Setup parser with stock header and lines to parse
    LINES = "2\t321681\t.\tC\t<R>\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\n"
    p = vcf_parser(LINES)
    p.parse_header()
    # Perform the actual test
    if sys.version_info < (3, 6):
        EXPECTED = (
            "Record('2', 321681, [], 'C', [SymbolicAllele('R')], None, ['PASS'], OrderedDict(), ['GT'], "
            "[Call('NA00001', OrderedDict([('GT', '0/1')])), Call('NA00002', OrderedDict([('GT', '0/0')])), "
            "Call('NA00003', OrderedDict([('GT', '0/0')]))])"
        )
    else:
        EXPECTED = (
            "Record('2', 321681, [], 'C', [SymbolicAllele('R')], None, ['PASS'], {}, ['GT'], "
            "[Call('NA00001', {'GT': '0/1'}), Call('NA00002', {'GT': '0/0'}), Call('NA00003', {'GT': '0/0'})])"
        )
    rec = p.parse_next_record()
    assert str(rec) == EXPECTED
    assert rec.ALT[0].serialize() == "<R>"
