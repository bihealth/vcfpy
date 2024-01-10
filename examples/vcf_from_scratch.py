#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create VCF from scratch using vcfpy."""

import sys

import vcfpy


def main():
    if len(sys.argv) != 2:
        print("Usage: vcf_from_scratch.py OUTPUT.vcf", file=sys.stderr)
        return 1

    header = vcfpy.Header(samples=vcfpy.SamplesInfos([]))
    with vcfpy.Writer.from_path(sys.argv[1], header) as writer:
        record = vcfpy.Record(
            CHROM="1", POS=1, ID=[], REF="N", ALT=[], QUAL=None, FILTER=[], INFO={}, FORMAT=[]
        )
        writer.write_record(record)


if __name__ == "__main__":
    sys.exit(main())
