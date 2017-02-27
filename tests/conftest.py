# -*- coding: utf-8 -*-
"""Shared code for tests"""

import textwrap

import pytest


@pytest.fixture
def multisample_vcf():
    """Return string with multi-sample VCF"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\t.\tG\tA\t29\t.\t.\tGT\t0/0\t1/0\t1/1
        20\t17330\t.\tT\tA\t3\t.\t.\tGT\t0/0\t0/1\t0/0
        20\t1110696\t.\tA\tG,T\t67\t.\t.\tGT\t1/2\t2/1\t2/2
        20\t1230237\t.\tT\t.\t47\t.\t.\tGT\t0/0\t0/0\t0/0
        20\t1234567\t.\tGTC\tG,GTCT\t50\t.\t.\tGT\t0/1\t0/2\t1/1
        """).lstrip()


@pytest.fixture
def multisample_vcf_only_NA00002():
    """Return string with multi-sample VCF limited to NA00002"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002
        20\t14370\t.\tG\tA\t29\t.\t.\tGT\t1/0
        20\t17330\t.\tT\tA\t3\t.\t.\tGT\t0/1
        20\t1110696\t.\tA\tG,T\t67\t.\t.\tGT\t2/1
        20\t1230237\t.\tT\t.\t47\t.\t.\tGT\t0/0
        20\t1234567\t.\tGTC\tG,GTCT\t50\t.\t.\tGT\t0/2
        """).lstrip()


@pytest.fixture
def multisample_vcf_reordered():
    """Return string with re-ordered multi-sample VCF"""
    return textwrap.dedent("""
        ##fileformat=VCFv4.3
        ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00002\tNA00003\tNA00001
        20\t14370\t.\tG\tA\t29\t.\t.\tGT\t1/0\t1/1\t0/0
        20\t17330\t.\tT\tA\t3\t.\t.\tGT\t0/1\t0/0\t0/0
        20\t1110696\t.\tA\tG,T\t67\t.\t.\tGT\t2/1\t2/2\t1/2
        20\t1230237\t.\tT\t.\t47\t.\t.\tGT\t0/0\t0/0\t0/0
        20\t1234567\t.\tGTC\tG,GTCT\t50\t.\t.\tGT\t0/2\t1/1\t0/1
        """).lstrip()


@pytest.fixture
def multisample_vcf_file(tmpdir, multisample_vcf):
    """Return path to multi-sample VCF file"""
    p = tmpdir.mkdir('input').join('input.vcf')
    with open(str(p), 'wt') as f:
        f.write(multisample_vcf)
    return str(p)
