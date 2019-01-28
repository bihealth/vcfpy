# -*- coding: utf-8 -*-
"""Test parsing of full VCF record lines
"""

import argparse
import io
import sys
import statistics
import time

from vcfpy import parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


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

LINE = """
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
""".lstrip()


def run(args):
    # Setup parser
    p = parser.VCFParser(io.StringIO(HEADER), "<builtin>")
    # Parse header
    p.parse_header()
    # Parse line several times
    times = []
    for r in range(args.repetitions):
        begin = time.clock()
        for _ in range(args.line_count):
            r = p._record_parser.parse_line(LINE)  # noqa
            if args.debug:
                print(r, file=sys.stderr)
        times.append(time.clock() - begin)
    print(
        "Took {:.3} seconds (stdev {:.3})".format(
            statistics.mean(times), statistics.stdev(times)
        ),
        file=sys.stderr,
    )


def main(argv=None):
    """Main program entry point for parsing command line arguments"""
    parser = argparse.ArgumentParser(description="Parser benchmark")

    parser.add_argument(
        "--debug", default=False, action="store_true", help="Enable debugging"
    )
    parser.add_argument(
        "--repetitions", type=int, default=10, help="Number of repetitions"
    )
    parser.add_argument(
        "--line-count", type=int, default=5000, help="Number of lines to parse"
    )

    args = parser.parse_args(argv)
    run(args)


if __name__ == "__main__":
    sys.exit(main())
