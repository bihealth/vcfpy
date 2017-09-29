# -*- coding: utf-8 -*-
"""Test writing of BGZF files
"""

import codecs
import gzip
import io

import pytest

import vcfpy
from vcfpy import parser
from vcfpy import writer
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
##INFO=<ID=ANNO,Number=.,Type=String,Description="Additional annotation">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Call-wise filters">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


@pytest.fixture(scope='function')
def header_samples():
    p = parser.Parser(stream=io.StringIO(MEDIUM_HEADER), path='<builtin>')
    p.parse_header()
    return (p.header, p.samples)


def check_file(path, line):
    """Test whether the file is the expected ``.vcf.gz`` file
    """
    raw = path.read(mode='rb')
    assert raw[0] == 0x1f
    assert raw[1] == 0x8b
    # compare actual result with expected
    inflated = gzip.decompress(raw)
    RESULT = codecs.latin_1_decode(inflated)[0]
    LINE = '20\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
    EXPECTED = MEDIUM_HEADER + LINE
    assert EXPECTED == RESULT


def test_write_minimal_record_writer_from_path(header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_header').join('out.vcf.gz')
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
    # check the resulting record
    LINE = '20\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
    check_file(path, LINE)


def test_write_minimal_record_writer_from_stream_path(
        header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_header').join('out.vcf.gz')
    header, _ = header_samples
    with open(str(path), 'wb') as f:
        w = writer.Writer.from_stream(f, header, path=str(path))
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
    # check the resulting record
    LINE = '20\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
    check_file(path, LINE)


def test_write_minimal_record_writer_from_stream_use_bgzf(
        header_samples, tmpdir_factory):
    O = vcfpy.OrderedDict
    # open temporary file and setup the Writer with header
    path = tmpdir_factory.mktemp('write_header').join('out.vcf.gz')
    header, samples = header_samples
    with open(str(path), 'wb') as f:
        w = writer.Writer.from_stream(f, header, samples, use_bgzf=True)
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
    # check the resulting record
    LINE = '20\t100\t.\tC\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n'
    check_file(path, LINE)
