#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple utility that allows for benchmarking the vcfpy module
"""

import argparse
import sys
import time


import vcfpy


def run(args):
    """Main program entry point after parsing arguments"""
    # open VCF reader
    reader = vcfpy.VCFReader.from_path(args.input_vcf)
    # optionally, open VCF writer
    writer = None
    if args.output_vcf:
        writer = vcfpy.VCFReader.from_path(
            reader.header, reader.samples, args.output_vcf)
    # read through input VCF file, optionally also writing out
    start = time.clock()
    for num, r in enumerate(reader):
        if num % 10000 == 0:
            print(num, ''.join(map(str, [r.CHROM, ':', r.POS])), sep='\t',
                  file=sys.stderr)
        if writer:
            writer.write_record(r)
        if args.max_records and num >= args.max_records:
            break
    end = time.clock()
    print('Read {} records in {} seconds'.format(num, (end - start)),
          file=sys.stderr)


def main(argv=None):
    """Main program entry point for parsing command line arguments"""
    parser = argparse.ArgumentParser(description='Benchmark driver')

    parser.add_argument('--max-records', type=int, default=100*1000)
    parser.add_argument('--input-vcf', type=str, required=True,
                        help='Path to VCF file to read')
    parser.add_argument('--output-vcf', type=str, required=False,
                        help='Path to VCF file to write if given')

    args = parser.parse_args(argv)
    run(args)


if __name__ == '__main__':
    sys.exit(main())