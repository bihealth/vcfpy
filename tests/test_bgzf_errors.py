# -*- coding: utf-8 -*-
"""Test BGZF error conditions and edge cases"""

import io

import pytest

from vcfpy.bgzf import BgzfWriter


def test_bgzf_writer_truncate_not_supported():
    """Test that truncate operation raises OSError"""
    stream = io.BytesIO()
    writer = BgzfWriter(fileobj=stream)

    with pytest.raises(OSError, match="truncate not supported on BGZF files"):
        writer.truncate()

    writer.close()


def test_bgzf_writer_writelines():
    """Test writelines method"""
    stream = io.BytesIO()
    writer = BgzfWriter(fileobj=stream)

    lines = ["line1\n", "line2\n", "line3\n"]
    writer.writelines(lines)

    # Check that content was written before closing
    writer.flush()
    assert len(stream.getvalue()) > 0

    writer.close()


def test_bgzf_writer_context_manager():
    """Test BGZF writer as context manager"""
    stream = io.BytesIO()

    with BgzfWriter(fileobj=stream) as writer:
        writer.write("test content\n")
        writer.flush()
        # Check content before context exits
        assert len(stream.getvalue()) > 0


def test_bgzf_writer_fileno():
    """Test fileno method with actual file"""
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(delete=False) as f:
        temp_path = f.name

    try:
        with BgzfWriter(temp_path) as writer:
            # Test fileno returns an integer
            fileno = writer.fileno()
            assert isinstance(fileno, int)
            assert fileno >= 0
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def test_bgzf_writer_flush():
    """Test flush method"""
    stream = io.BytesIO()
    writer = BgzfWriter(fileobj=stream)

    writer.write("test content")
    writer.flush()  # Should not raise any errors

    writer.close()


def test_bgzf_writer_multiple_writes():
    """Test multiple write operations"""
    stream = io.BytesIO()
    writer = BgzfWriter(fileobj=stream)

    # Write multiple chunks
    for i in range(10):
        writer.write(f"line {i}\n")

    # Check content before closing
    writer.flush()
    assert len(stream.getvalue()) > 0

    writer.close()
