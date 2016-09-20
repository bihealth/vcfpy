# -*- coding: utf-8 -*-
"""Reading of VCF files "from the wild"

The aim here is to test whether they can be parsed at all or make the parser
crash.
"""

import os

import pytest

from vcfpy import reader

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


FILENAMES = [
    '1kg.sites.vcf',
    '1kg.vcf.gz',
    'bcftools.vcf',
    'contig_idonly.vcf',
    'example-4.0.vcf',
    'example-4.1-info-multiple-values.vcf',
    'example-4.1-ploidy.vcf',
    'example-4.1-sv.vcf',
    'example-4.1.vcf',
    'example-4.2.vcf',
    'freebayes.vcf',
    'gatk_26_meta.vcf',
    'gatk.vcf',
    'gonl.chr20.release4.gtc.vcf',
    'info-type-character.vcf',
    'issue-140-file1.vcf',
    'issue-140-file2.vcf',
    'issue-140-file3.vcf',
    'issue-16.vcf',
    'issue-201.vcf.gz',
    'issue-214.vcf',
    'issue_49.vcf',
    'metadata-whitespace.vcf',
    'mixed-filtering.vcf',
    'null_genotype_mono.vcf',
    'parse-meta-line.vcf',
    'samples-space.vcf',
    'samtools.vcf',
    'strelka.vcf',
    'string_as_flag.vcf',
    'tb.vcf.gz',
    'uncalled_genotypes.vcf',
    'walk_left.vcf',
    'walk_refcall.vcf'
]


@pytest.mark.parametrize('filename', FILENAMES)
def test_read_vcfs_from_the_wild(filename):
    path = os.path.join(os.path.dirname(__file__), 'vcfs_from_the_wild',
                        filename)
    r = reader.Reader.from_path(path)
