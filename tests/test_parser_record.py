# -*- coding: utf-8 -*-
"""Test parsing of full VCF record lines
"""

import io

import pytest

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
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=FT,Number=.,Type=String,Description="Haplotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


def vcf_parser(lines):
    return parser.VCFParser(io.StringIO(MEDIUM_HEADER + lines), '<builtin>')


def test_parse_minimal_record():
    # Setup parser with stock header and lines to parse
    LINES = '20\t1\t.\tC\tG\t.\t.\t.\tGT\t0/1\t0/2\t.\n'
    import sys; print(LINES, file=sys.stderr)
    p = vcf_parser(LINES)
    hdr = p.parse_header()
    # Perform the actual test
    EXPECTED = (
        "Record('20', 1, [], 'C', [Substitution(type_='SNV', value='G')], "
        "None, [], OrderedDict(), ['GT'], "
        "[Call('NA00001', OrderedDict([('GT', '0/1')])),"
        " Call('NA00002', OrderedDict([('GT', '0/2')])),"
        " Call('NA00003', OrderedDict([('GT', None)]))])")
    RESULT = p.parse_next_record()
    assert str(RESULT) == EXPECTED
