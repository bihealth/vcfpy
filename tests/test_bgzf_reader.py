"""Test the BGZF reader functionality."""

import io
import pathlib
import tempfile

import pytest

from vcfpy.bgzf import BgzfReader, make_virtual_offset, split_virtual_offset


def test_bgzf_reader_basic_functionality():
    """Test basic BgzfReader functionality."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

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
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

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


def test_bgzf_reader_constructor_errors():
    """Test BgzfReader constructor error conditions."""
    # Test providing both filename and fileobj
    with tempfile.NamedTemporaryFile() as tmp:
        with pytest.raises(ValueError, match="Supply either filename or fileobj, not both"):
            BgzfReader(filename=tmp.name, fileobj=tmp)

    # Test invalid modes (covered with pragma: no cover)
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"
    if vcf_path.exists():
        with pytest.raises(ValueError, match="Must use a read mode"):
            BgzfReader(str(vcf_path), mode="w")

    # Test max_cache < 1 (covered with pragma: no cover)
    if vcf_path.exists():
        with pytest.raises(ValueError, match="Use max_cache with a minimum of 1"):
            BgzfReader(str(vcf_path), max_cache=0)

    # Test fileobj not in binary mode (covered with pragma: no cover)
    with io.StringIO("test") as text_file:
        with pytest.raises(ValueError, match="fileobj not opened in binary mode"):
            BgzfReader(fileobj=text_file)

    # Test None filename when no fileobj (covered with pragma: no cover)
    with pytest.raises(ValueError, match="Must provide filename if fileobj is None"):
        BgzfReader(filename=None, fileobj=None)


def test_bgzf_reader_read_method():
    """Test BgzfReader read method with various sizes."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    reader = BgzfReader(str(vcf_path))

    # Test reading negative size (covered with pragma: no cover)
    with pytest.raises(NotImplementedError, match="Don't be greedy"):
        reader.read(-1)

    # Test reading specific sizes
    data1 = reader.read(10)
    assert len(data1) == 10
    assert data1 == "##fileform"

    data2 = reader.read(5)
    assert len(data2) == 5
    assert data2 == "at=VC"

    reader.close()


def test_bgzf_reader_iteration():
    """Test BgzfReader iteration functionality."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    reader = BgzfReader(str(vcf_path))

    # Test __iter__ method (covered with pragma: no cover)
    assert iter(reader) is reader

    # Test __next__ method
    first_line = next(reader)
    assert first_line.startswith("##fileformat=VCFv4.3")

    # Test readlines with hint
    reader.seek(0)
    lines = reader.readlines(hint=3)
    assert len(lines) == 3
    assert all(line.startswith("#") for line in lines)

    # Test readlines without hint (read all)
    reader.seek(0)
    all_lines = reader.readlines()
    assert len(all_lines) > 10  # Should have many lines

    reader.close()


def test_bgzf_reader_unsupported_operations():
    """Test BgzfReader operations that should raise errors."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    reader = BgzfReader(str(vcf_path))

    # Test write operations (covered with pragma: no cover)
    with pytest.raises(OSError, match="not writable"):
        reader.write("test")

    with pytest.raises(OSError, match="not writable"):
        reader.writelines(["test\n"])

    with pytest.raises(OSError, match="not writable"):
        reader.truncate(10)

    # Test seekable (covered with pragma: no cover)
    assert reader.seekable() is True

    # Test flush (covered with pragma: no cover)
    reader.flush()  # Should be a no-op

    reader.close()


def test_bgzf_reader_context_manager():
    """Test BgzfReader as context manager."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    # Test context manager functionality
    with BgzfReader(str(vcf_path)) as reader:
        first_line = reader.readline()
        assert first_line.startswith("##fileformat=VCFv4.3")

        # Test name property
        assert reader.name.endswith("multi_contig.vcf.gz")

    # File should be closed after context manager


def test_bgzf_reader_with_fileobj():
    """Test BgzfReader with file object instead of filename."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    # Test using fileobj parameter
    with open(vcf_path, "rb") as f:
        reader = BgzfReader(fileobj=f)
        first_line = reader.readline()
        assert first_line.startswith("##fileformat=VCFv4.3")
        reader.close()


def test_bgzf_reader_edge_cases():
    """Test BgzfReader edge cases for better coverage."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Test VCF file not found")

    reader = BgzfReader(str(vcf_path))

    # Read some data to get past the first block
    lines = []
    for _ in range(10):
        line = reader.readline()
        lines.append(line)

    # Test readline method (note: size parameter is not implemented in current version)
    reader.seek(0)
    first_line = reader.readline()
    assert first_line.startswith("##fileformat=VCFv4.3")
    assert first_line.endswith("\n")

    # Test reading 0 bytes
    reader.seek(0)
    empty_data = reader.read(0)
    assert empty_data == ""

    reader.close()


if __name__ == "__main__":
    test_bgzf_reader_basic_functionality()
    test_bgzf_virtual_offsets()
    test_bgzf_reader_seek_functionality()
    test_bgzf_reader_constructor_errors()
    test_bgzf_reader_read_method()
    test_bgzf_reader_iteration()
    test_bgzf_reader_unsupported_operations()
    test_bgzf_reader_context_manager()
    test_bgzf_reader_with_fileobj()
    test_bgzf_reader_edge_cases()
    print("All BGZF tests passed!")
