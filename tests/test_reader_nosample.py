# -*- coding: utf-8 -*-
"""Tests for reading with file without samples"""

import textwrap

import vcfpy


def test_reading_parse_nosample_read_only(tmpdir, nosample_vcf_file):
    """Test reading of VCF file without any samples."""
    # Perform record-wise copying, saving results in records
    records = []
    with vcfpy.Reader.from_path(nosample_vcf_file) as reader:
        for record in reader:
            records.append(record)
    # Check resulting records, checking the first and last records is enough
    assert len(records) == 5


def test_reading_parse_nosample_also_write(tmpdir, nosample_vcf_file):
    """Read VCF file without samples, write file with samples."""
    # Perform record-wise copying, saving results in records
    path_out = tmpdir.mkdir("output").join("output.vcf")
    with vcfpy.Reader.from_path(nosample_vcf_file) as reader:
        header = reader.header.copy()
        header.samples = vcfpy.SamplesInfos(["NA00001", "NA00002", "NA00003"])
        with vcfpy.Writer.from_path(str(path_out), header) as writer:
            for record in reader:
                record.update_calls(
                    [vcfpy.Call(sample, {}) for sample in ("NA00001", "NA00002", "NA00003")]
                )
                record.add_format("GT", "./.")
                writer.write_record(record)

    expected = textwrap.dedent(
        """
    ##fileformat=VCFv4.3
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
    20	14370	.	G	A	29	.	.	GT	.	.	.
    20	17330	.	T	A	3	.	.	GT	.	.	.
    20	1110696	.	A	G,T	67	.	.	GT	.	.	.
    20	1230237	.	T	.	47	.	.	GT	.	.	.
    20	1234567	.	GTC	G,GTCT	50	.	.	GT	.	.	.
    """
    ).lstrip()

    assert path_out.open("rt").read() == expected
