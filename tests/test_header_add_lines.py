# -*- coding: utf-8 -*-
"""Test adding header lines to headers
"""

import io

import pytest

import vcfpy
from vcfpy import parser
from vcfpy import header

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

HEADER = """
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
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
""".lstrip()


@pytest.fixture(scope='function')
def vcf_header():
    p = parser.Parser(io.StringIO(HEADER))
    return p.parse_header()


def test_add_simple_line(vcf_header):
    assert len(vcf_header.lines) == 18
    line = header.HeaderLine('somekey', 'somevalue')
    vcf_header.add_line(line)
    assert len(vcf_header.lines) == 19

    assert vcf_header.lines[-1].key == 'somekey'
    assert vcf_header.lines[-1].value == 'somevalue'


def test_add_contig_line(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18
    assert '20a' not in vcf_header._indices['contig']

    # add header line
    VALUE = (
        '<ID=20a,length=62435964,assembly=B36,'
        'md5=f126cdf8a6e0c7f379d618ff66beb2da,'
        'species="Homo sapiens",taxonomy=x>')
    line = header.ContigHeaderLine(
        'contig', VALUE, vcfpy.OrderedDict([
            ('ID', '20a'),
            ('length', 62435964),
            ('assembly', 'B36'),
            ('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'),
            ('species', 'Homo sapiens'),
            ('taxonomy', 'x'),
        ]))
    vcf_header.add_line(line)

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert '20a' in vcf_header._indices['contig']
    assert vcf_header._indices['contig']['20a'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'contig'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 6
    assert vcf_header.lines[-1].mapping['ID'] == '20a'
    assert vcf_header.lines[-1].mapping['length'] == 62435964
    assert vcf_header.lines[-1].mapping['assembly'] == 'B36'
    assert vcf_header.lines[-1].mapping['md5'] == 'f126cdf8a6e0c7f379d618ff66beb2da'
    assert vcf_header.lines[-1].mapping['species'] == 'Homo sapiens'
    assert vcf_header.lines[-1].mapping['taxonomy'] == 'x'


def test_add_contig_line_shortcut(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18
    assert '20a' not in vcf_header._indices['contig']

    # add header line
    mapping = vcfpy.OrderedDict([
        ('ID', '20a'),
        ('length', 62435964),
        ('assembly', 'B36'),
        ('md5', 'f126cdf8a6e0c7f379d618ff66beb2da'),
        ('species', 'Homo sapiens'),
        ('taxonomy', 'x')])
    vcf_header.add_contig_line(mapping)

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert '20a' in vcf_header._indices['contig']
    assert vcf_header._indices['contig']['20a'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'contig'
    VALUE = (
        '<ID=20a,length=62435964,assembly=B36,'
        'md5=f126cdf8a6e0c7f379d618ff66beb2da,'
        'species="Homo sapiens",taxonomy=x>')
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 6
    assert vcf_header.lines[-1].mapping['ID'] == '20a'
    assert vcf_header.lines[-1].mapping['length'] == 62435964
    assert vcf_header.lines[-1].mapping['assembly'] == 'B36'
    assert vcf_header.lines[-1].mapping['md5'] == 'f126cdf8a6e0c7f379d618ff66beb2da'
    assert vcf_header.lines[-1].mapping['species'] == 'Homo sapiens'
    assert vcf_header.lines[-1].mapping['taxonomy'] == 'x'


def test_add_filter_line(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18

    # add header line
    VALUE = '<ID=q10a,Description="Quality below 10">'
    line = header.FilterHeaderLine(
        'FILTER', VALUE, vcfpy.OrderedDict(
            [('ID', 'q10a'), ('Description', 'Quality below 10')]))
    vcf_header.add_line(line)

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert 'q10a' in vcf_header._indices['FILTER']
    assert vcf_header._indices['FILTER']['q10a'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'FILTER'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 2
    assert vcf_header.lines[-1].mapping['ID'] == 'q10a'
    assert vcf_header.lines[-1].mapping['Description'] == 'Quality below 10'


def test_add_filter_line_shortcut(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18

    # add header line
    mapping = vcfpy.OrderedDict(
        [('ID', 'q10a'), ('Description', 'Quality below 10')])
    vcf_header.add_filter_line(mapping)

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert 'q10a' in vcf_header._indices['FILTER']
    assert vcf_header._indices['FILTER']['q10a'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'FILTER'
    VALUE = '<ID=q10a,Description="Quality below 10">'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 2
    assert vcf_header.lines[-1].mapping['ID'] == 'q10a'
    assert vcf_header.lines[-1].mapping['Description'] == 'Quality below 10'


def test_add_format_line(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18

    # add header line
    VALUE = '<ID=GTa,Number=1,Type=String,Description="Genotype">'
    line = header.FormatHeaderLine(
        'FORMAT', VALUE, vcfpy.OrderedDict(
            [('ID', 'GTa'),
             ('Number', 1),
             ('Type', 'String'),
             ('Description', 'Genotype')]))
    vcf_header.add_line(line)

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert 'GTa' in vcf_header._indices['FORMAT']
    assert vcf_header._indices['FORMAT']['GTa'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'FORMAT'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 4
    assert vcf_header.lines[-1].mapping['ID'] == 'GTa'
    assert vcf_header.lines[-1].mapping['Number'] == 1
    assert vcf_header.lines[-1].mapping['Type'] == 'String'
    assert vcf_header.lines[-1].mapping['Description'] == 'Genotype'


def test_add_format_line_shortcut(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18

    # add header line
    mapping = vcfpy.OrderedDict(
        [('ID', 'GTa'),
         ('Number', 1),
         ('Type', 'String'),
         ('Description', 'Genotype')])
    vcf_header.add_format_line(mapping)

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert 'GTa' in vcf_header._indices['FORMAT']
    assert vcf_header._indices['FORMAT']['GTa'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'FORMAT'
    VALUE = '<ID=GTa,Number=1,Type=String,Description="Genotype">'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 4
    assert vcf_header.lines[-1].mapping['ID'] == 'GTa'
    assert vcf_header.lines[-1].mapping['Number'] == 1
    assert vcf_header.lines[-1].mapping['Type'] == 'String'
    assert vcf_header.lines[-1].mapping['Description'] == 'Genotype'


def test_add_info_line(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18

    # add header line
    VALUE = '<ID=DPa,Number=1,Type=Integer,Description="Total Depth">'
    line = header.FormatHeaderLine(
        'INFO', VALUE, vcfpy.OrderedDict(
            [('ID', 'DPa'),
             ('Number', 1),
             ('Type', 'Integer'),
             ('Description', 'Total Depth')]))
    vcf_header.add_line(line)
    assert len(vcf_header.lines) == 19

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert 'DPa' in vcf_header._indices['INFO']
    assert vcf_header._indices['INFO']['DPa'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'INFO'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 4
    assert vcf_header.lines[-1].mapping['ID'] == 'DPa'
    assert vcf_header.lines[-1].mapping['Number'] == 1
    assert vcf_header.lines[-1].mapping['Type'] == 'Integer'
    assert vcf_header.lines[-1].mapping['Description'] == 'Total Depth'


def test_add_info_line_shortcut(vcf_header):
    # check header before adding
    assert len(vcf_header.lines) == 18

    # add header line
    VALUE = '<ID=DPa,Number=1,Type=Integer,Description="Total Depth">'
    mapping = vcfpy.OrderedDict(
        [('ID', 'DPa'),
         ('Number', 1),
         ('Type', 'Integer'),
         ('Description', 'Total Depth')])
    vcf_header.add_info_line(mapping)
    assert len(vcf_header.lines) == 19

    # check header after adding
    assert len(vcf_header.lines) == 19
    assert 'DPa' in vcf_header._indices['INFO']
    assert vcf_header._indices['INFO']['DPa'] is vcf_header.lines[-1]

    # Check resulting added header line
    assert vcf_header.lines[-1].key == 'INFO'
    assert vcf_header.lines[-1].value == VALUE
    assert len(vcf_header.lines[-1].mapping) == 4
    assert vcf_header.lines[-1].mapping['ID'] == 'DPa'
    assert vcf_header.lines[-1].mapping['Number'] == 1
    assert vcf_header.lines[-1].mapping['Type'] == 'Integer'
    assert vcf_header.lines[-1].mapping['Description'] == 'Total Depth'
