"""Test the BGZF reader functionality."""

import pathlib

import pytest

from vcfpy.bgzf import BgzfReader, make_virtual_offset, split_virtual_offset


def test_bgzf_reader_basic_functionality():
    """Test basic BgzfReader functionality."""
    vcf_path = pathlib.Path(__file__).parent.parent / "tests" / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    reader = BgzfReader(str(vcf_path))

    # Test reading first line
    first_line = reader.readline()
    assert first_line.startswith("##fileformat=VCFv4.3")

    # Test tell position
    pos = reader.tell()
    assert pos > 0

    # Test seeking back to start
    reader.seek(0)
    assert reader.tell() == 0

    # Read first line again to verify seek worked
    first_line_again = reader.readline()
    assert first_line == first_line_again

    reader.close()


def test_bgzf_virtual_offsets():
    """Test virtual offset manipulation functions."""
    # Test make_virtual_offset and split_virtual_offset
    test_cases = [
        (0, 0, 0),
        (0, 1, 1),
        (0, 65535, 65535),
        (1, 0, 65536),
        (1, 1, 65537),
        (100000, 10, 6553600010),
    ]

    for block_start, within_block, expected_voffset in test_cases:
        # Test make_virtual_offset
        voffset = make_virtual_offset(block_start, within_block)
        assert voffset == expected_voffset

        # Test split_virtual_offset
        split_block, split_within = split_virtual_offset(voffset)
        assert split_block == block_start
        assert split_within == within_block


def test_bgzf_reader_seek_functionality():
    """Test BgzfReader seek functionality with virtual offsets."""
    vcf_path = pathlib.Path(__file__).parent.parent / "tests" / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    reader = BgzfReader(str(vcf_path))

    # Record positions and content
    positions_and_content = []
    for _ in range(3):
        pos = reader.tell()
        line = reader.readline()
        positions_and_content.append((pos, line))

    # Test seeking back to each recorded position
    for pos, expected_line in positions_and_content:
        reader.seek(pos)
        actual_line = reader.readline()
        assert actual_line == expected_line

    reader.close()


if __name__ == "__main__":
    test_bgzf_reader_basic_functionality()
    test_bgzf_virtual_offsets()
    test_bgzf_reader_seek_functionality()
    print("All BGZF tests passed!")
