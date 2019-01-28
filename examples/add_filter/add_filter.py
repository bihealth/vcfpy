#!/usr/bin/env python
# -*- coding: utf-8 -*-
import vcfpy

# Open input, add FILTER header, and open output file
reader = vcfpy.Reader.from_path("input.vcf")
reader.header.add_filter_line(
    vcfpy.OrderedDict([("ID", "DP10"), ("Description", "total DP < 10")])
)
writer = vcfpy.Writer.from_path("/dev/stdout", reader.header)

# Add "DP10" filter to records having less than 10 reads
for record in reader:
    ad = sum(c.data.get("DP", 0) for c in record.calls)
    if ad < 10:
        record.add_filter("DP10")
    writer.write_record(record)
