# -*- coding: utf-8 -*-
"""Writing of VCF records
"""

import io

import pytest

import vcfpy
from vcfpy import parser
from vcfpy import writer
from vcfpy import record

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# TODO: cleanup, refactor tests somewhat

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
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=FT,Number=.,Type=String,Description="Call-wise filters">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


@pytest.fixture(scope='function')
def header_samples():
    p = parser.Parser(stream=io.StringIO(MEDIUM_HEADER), path='<builtin>')
    p.parse_header()
    return (p.header, p.samples)


def test_write_minimal_record(header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_header').join('out.vcf')
    header, _ = header_samples
    w = writer.Writer.from_path(path, header)
    # construct record to write out from scratch
    r = record.Record(
        '20', 100, [], 'C', [record.Substitution(record.SNV, 'T')],
        None, [], O(), ['GT'],
        [
            record.Call('NA00001', O(GT='0/1')),
            record.Call('NA00002', O(GT='0/0')),
            record.Call('NA00003', O(GT='1/1')),
        ])
    # write out the record, close file to ensure flushing to disk
    w.write_record(r)
    w.close()
    # compare actual result with expected
    RESULT = path.read()
    LINE = '20\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
    EXPECTED = MEDIUM_HEADER + LINE
    assert EXPECTED == RESULT


def test_write_annotated_record(header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    S = record.Substitution
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_annotated_record').join('out.vcf')
    header, _ = header_samples
    w = writer.Writer.from_path(path, header)
    # construct record to write out from scratch
    r = record.Record(
        '20',
        100,
        ['rs333', 'CSN42'],
        'C',
        [
            record.Substitution(record.SNV, 'T'),
            record.Substitution(record.SNV, 'G'),
        ],
        50,
        ['PASS'],
        O([('DP', 93), ('AF', [0.3, 0.2]), ('DB', True)]),
        ['GT', 'DP', 'GQ', 'HQ'],
        [
            record.Call('NA00001', O(GT='0/1', DP=30, GQ=40, HQ=[1, 2])),
            record.Call('NA00002', O(GT='0/2', DP=31, GQ=41, HQ=[3, 4])),
            record.Call('NA00003', O(GT='1/2', DP=32, GQ=42, HQ=[5, 6])),
        ])
    # write out the record, close file to ensure flushing to disk
    w.write_record(r)
    w.close()
    # compare actual result with expected
    RESULT = path.read()
    LINE = '20\t100\trs333;CSN42\tC\tT,G\t50\tPASS\tDP=93;AF=0.3,0.2;DB\tGT:DP:GQ:HQ\t0/1:30:40:1,2\t0/2:31:41:3,4\t1/2:32:42:5,6\n'
    EXPECTED = MEDIUM_HEADER + LINE
    assert EXPECTED == RESULT


def test_write_record_with_escaping(header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    S = record.Substitution
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_header').join('out.vcf')
    header, _ = header_samples
    w = writer.Writer.from_path(path, header)
    # construct record to write out from scratch
    r = record.Record(
        '20', 100, [], 'C', [record.Substitution(record.SNV, 'T')],
        None, [],
        O([
            ('ANNO', ['Here,are%some chars', '%25'])
        ]),
        ['GT', 'FT'],
        [
            record.Call('NA00001', O(GT='0/1', FT=['%25', 'FOO'])),
            record.Call('NA00002', O(GT='0/0', FT=[])),
            record.Call('NA00003', O(GT='1/1', FT=[])),
        ])
    # write out the record, close file to ensure flushing to disk
    w.write_record(r)
    w.close()
    # compare actual result with expected
    RESULT = path.read()
    LINE = ('20\t100\t.\tC\tT\t.\t.\tANNO=Here%2Care%25some chars,'
            '%2525\tGT:FT\t0/1:%2525;FOO\t0/0:.\t1/1:.\n')
    EXPECTED = MEDIUM_HEADER + LINE
    assert EXPECTED == RESULT
