"""Additional tests to push coverage to 95%+ by targeting specific missing lines"""

import pytest

import vcfpy
from vcfpy import header


def test_writer_path_error():
    """Test writer when path has issues"""
    header_obj = header.Header()

    # Try to write to invalid path
    with pytest.raises(IOError):
        vcfpy.Writer.from_path("/invalid/path/that/does/not/exist.vcf", header_obj)
