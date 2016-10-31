# -*- coding: utf-8 -*-
"""Tests for the Reader and Writer working as context managers
"""

import io
import os

import pytest

import vcfpy
from vcfpy import reader
from vcfpy import writer
from vcfpy import parser
from vcfpy import record

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


def test_reader():
    path = os.path.join(os.path.dirname(__file__), 'vcfs/full_vcf43.vcf')
    with reader.Reader.from_path(path) as r:
        # should be open now
        assert r.stream
        assert not r.stream.closed
        # read records
        records = [rec for rec in r]
    # should be closed now
    assert r.stream
    assert r.stream.closed
    # check result
    assert len(records) == 5


def test_reader_fetch():
    path = os.path.join(os.path.dirname(__file__), 'vcfs',
                        'multi_contig.vcf.gz')
    with reader.Reader.from_path(path) as r:
        records = [vcf_rec for vcf_rec in r.fetch('20', 1110695, 1230236)]
        assert r.stream
        assert not r.stream.closed
        assert r.tabix_file
        assert not r.tabix_file.closed
    # closed
    assert r.stream
    assert r.stream.closed
    assert r.tabix_file
    assert r.tabix_file.closed
    # check result
    assert len(records) == 1
    assert records[0].CHROM == '20'
    assert records[0].POS == 1110696


def test_writer(header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_header').join('out.vcf')
    header, _ = header_samples
    # construct record to write out from scratch
    r = record.Record(
        '20', 100, [], 'C', [record.Substitution(record.SNV, 'T')],
        None, [], O(), ['GT'],
        [
            record.Call('NA00001', O(GT='0/1')),
            record.Call('NA00002', O(GT='0/0')),
            record.Call('NA00003', O(GT='1/1')),
        ])
    # open writer
    with writer.Writer.from_path(path, header) as w:
        # write out the record
        w.write_record(r)
    # should be closed
    assert w.stream
    assert w.stream.closed
    # compare actual result with expected
    RESULT = path.read()
    LINE = '20\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
    EXPECTED = MEDIUM_HEADER + LINE
    assert EXPECTED == RESULT
