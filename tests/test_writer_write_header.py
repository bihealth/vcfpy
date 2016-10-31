# -*- coding: utf-8 -*-
"""Writing of VCF headers
"""

import io

import pytest

from vcfpy import parser
from vcfpy import writer

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
##ALT=<ID=R,Description="IUPAC code R = A/G">
##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
##PEDIGREE=<ID=TumourSample,Original=GermlineID>
##SAMPLE=<ID=Sample1,Assay=WholeGenome,Ethnicity=AFR,Disease=None,Description="Patient germline genome from unaffected",DOI=url>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


@pytest.fixture(scope='function')
def header_samples():
    p = parser.Parser(stream=io.StringIO(MEDIUM_HEADER), path='<builtin>')
    p.parse_header()
    return (p.header, p.samples)


def test_write_header(header_samples, tmpdir_factory):
    path = tmpdir_factory.mktemp('write_header').join('out.vcf')
    header, _ = header_samples
    w = writer.Writer.from_path(path, header)
    w.close()
    RESULT = path.read()
    EXPECTED = MEDIUM_HEADER
    assert RESULT == EXPECTED
