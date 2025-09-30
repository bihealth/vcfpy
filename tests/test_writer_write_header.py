# -*- coding: utf-8 -*-
"""Writing of VCF headers"""

import io
import pathlib
import textwrap
from typing import cast

import pytest

from vcfpy import header, parser, writer

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


@pytest.fixture(scope="function")
def header_samples() -> tuple[header.Header, header.SamplesInfos]:
    cast_stream = cast(io.TextIOWrapper, io.StringIO(MEDIUM_HEADER))
    p = parser.Parser(stream=cast_stream, path="<builtin>")
    p.parse_header()
    assert p.header is not None
    assert p.samples is not None
    return (p.header, p.samples)


def test_write_header(header_samples: tuple[header.Header, header.SamplesInfos], tmpdir: pathlib.Path):
    path = tmpdir / "out.vcf"
    header, _ = header_samples
    w = writer.Writer.from_path(path, header)
    w.close()
    with path.open("rt") as f:
        result = f.read()
    expected = MEDIUM_HEADER
    assert result == expected


def test_write_header_no_samples(tmpdir: pathlib.Path):
    # Create header to write out from scratch
    hdr = header.Header(lines=[header.HeaderLine("fileformat", "VCFv4.0")], samples=header.SamplesInfos([]))
    # Write out header
    path = tmpdir / "out.vcf"
    w = writer.Writer.from_path(path, hdr)
    w.close()
    # Compare result
    with path.open("rt") as f:
        result = f.read()
    expected = textwrap.dedent(
        r"""
    ##fileformat=VCFv4.0
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    """
    ).lstrip()
    assert result == expected
