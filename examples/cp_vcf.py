#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple utility that allows for benchmarking the vcfpy module
"""

import argparse
import itertools
import sys
import time

import vcf
import vcfpy


class FakeWriter:
    """Used as a dummy object for writing"""

    def write_record(self, record):
        pass  # noop


class BaseRunner:

    def __init__(self, args):
        self.args = args

    def run(self):
        """Main program entry point after parsing arguments"""
        start = time.clock()
        it = iter(self.reader)
        if self.args.max_records:
            it = itertools.islice(it, self.args.max_records)
        num = self.work(it)
        end = time.clock()
        print('Read {} records in {} seconds'.format(num, (end - start)),
              file=sys.stderr)

    def work(self, it):
        for num, r in enumerate(it):
            if num % 10000 == 0:
                print(num, ''.join(map(str, [r.CHROM, ':', r.POS])), sep='\t',
                      file=sys.stderr)
            self.writer.write_record(r)
        return num


class VCFPyRunner(BaseRunner):

    def __init__(self, args):
        super().__init__(args)
        self.reader = vcfpy.Reader.from_path(args.input_vcf)
        self.writer = FakeWriter()
        if args.output_vcf:
            self.writer = vcfpy.VCFWriter.from_path(
                reader.header, reader.samples, args.output_vcf)


class PyVCFRunner(BaseRunner):

    def __init__(self, args):
        super().__init__(args)
        self.reader = vcf.Reader(filename=args.input_vcf)
        self.writer = FakeWriter()


def run_pyvcf(args):
    """Main program entry point after parsing arguments"""
    # open VCF reader
    reader = vcf.Reader(filename=args.input_vcf)
    # optionally, open VCF writer
    writer = None
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
    parser.add_argument('--engine', type=str, choices=('vcfpy', 'pyvcf'),
                        default='vcfpy')
    parser.add_argument('--input-vcf', type=str, required=True,
                        help='Path to VCF file to read')
    parser.add_argument('--output-vcf', type=str, required=False,
                        help='Path to VCF file to write if given')

    args = parser.parse_args(argv)
    if args.engine == 'vcfpy':
        VCFPyRunner(args).run()
    else:
        PyVCFRunner(args).run()


if __name__ == '__main__':
    sys.exit(main())
