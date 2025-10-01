# -*- coding: utf-8 -*-
"""Support code for writing BGZF files

Shamelessly taken from Biopython
"""

#                 Biopython License Agreement
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE.

import codecs
import struct
import typing
import zlib
from typing import Iterable

# For Python 2 can just use: _bgzf_magic = '\x1f\x8b\x08\x04'
# but need to use bytes on Python 3
_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_eof = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)
_bytes_BC = b"BC"


def make_virtual_offset(block_start_offset: int, within_block_offset: int) -> int:
    """Compute a BGZF virtual offset from block start and within block offsets.
    The BAM indexing scheme records read positions using a 64 bit
    'virtual offset', comprising in C terms:
    block_start_offset << 16 | within_block_offset
    Here block_start_offset is the file offset of the BGZF block
    start (unsigned integer using up to 64-16 = 48 bits), and
    within_block_offset within the (decompressed) block (unsigned
    16 bit integer).

    >>> make_virtual_offset(0, 0)
    0
    >>> make_virtual_offset(0, 1)
    1
    >>> make_virtual_offset(0, 2**16 - 1)
    65535
    >>> make_virtual_offset(0, 2**16)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= within_block_offset < 2**16, got 65536
    >>> 65536 == make_virtual_offset(1, 0)
    True
    >>> 65537 == make_virtual_offset(1, 1)
    True
    >>> 131071 == make_virtual_offset(1, 2**16 - 1)
    True
    >>> 6553600000 == make_virtual_offset(100000, 0)
    True
    >>> 6553600001 == make_virtual_offset(100000, 1)
    True
    >>> 6553600010 == make_virtual_offset(100000, 10)
    True
    >>> make_virtual_offset(2**48, 0)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= block_start_offset < 2**48, got 281474976710656
    """
    if within_block_offset < 0 or within_block_offset >= 65536:
        raise ValueError("Require 0 <= within_block_offset < 2**16, got %i" % within_block_offset)
    if block_start_offset < 0 or block_start_offset >= 281474976710656:
        raise ValueError("Require 0 <= block_start_offset < 2**48, got %i" % block_start_offset)
    return (block_start_offset << 16) | within_block_offset


def split_virtual_offset(virtual_offset: int) -> tuple[int, int]:
    """Split a 64-bit BGZF virtual offset into block start and within block offsets.

    Returns a tuple of (block_start_offset, within_block_offset).

    >>> split_virtual_offset(0)
    (0, 0)
    >>> split_virtual_offset(1)
    (0, 1)
    >>> split_virtual_offset(65535)
    (0, 65535)
    >>> split_virtual_offset(65536)
    (1, 0)
    >>> split_virtual_offset(65537)
    (1, 1)
    >>> split_virtual_offset(1195311108)
    (18239, 4)
    """
    start_offset = virtual_offset >> 16
    within_block = virtual_offset ^ (start_offset << 16)
    return start_offset, within_block


def _load_bgzf_block(handle: typing.IO[bytes]) -> tuple[int, str]:
    """Load the next BGZF block from the file handle.

    Returns a tuple of (block_size, decompressed_data_as_string).
    Raises StopIteration if EOF is reached.
    """
    # Read the complete gzip header (10 bytes)
    # ID1 ID2 CM FLG MTIME(4) XFL OS
    header = handle.read(10)
    if len(header) < 10:  # pragma: no cover
        raise StopIteration("EOF")

    # Check magic bytes
    if header[:4] != _bgzf_magic:  # pragma: no cover
        raise ValueError(f"Invalid BGZF magic: {header[:4]!r}")

    # Check FLG field - should have FEXTRA bit set (0x04)
    flg = header[3]
    if not (flg & 0x04):  # pragma: no cover
        raise ValueError("BGZF file missing FEXTRA flag")

    # Read XLEN (extra field length)
    xlen_bytes = handle.read(2)
    if len(xlen_bytes) < 2:  # pragma: no cover
        raise ValueError("Truncated BGZF extra field length")

    xlen = struct.unpack("<H", xlen_bytes)[0]
    if xlen < 6:  # pragma: no cover
        raise ValueError("Invalid BGZF extra field length")

    # Read the extra field
    extra_data = handle.read(xlen)
    if len(extra_data) < xlen:  # pragma: no cover
        raise ValueError("Truncated BGZF extra field")

    # Parse the BC subfield to get BSIZE
    if extra_data[:2] != _bytes_BC:  # pragma: no cover
        raise ValueError("Missing BC subfield in BGZF header")

    if len(extra_data) < 6:  # pragma: no cover
        raise ValueError("BC subfield too short")

    # BC subfield: SI1='B' SI2='C' SLEN=2 BSIZE(2 bytes)
    slen = struct.unpack("<H", extra_data[2:4])[0]
    if slen != 2:  # pragma: no cover
        raise ValueError(f"Expected BC subfield length 2, got {slen}")

    bsize = struct.unpack("<H", extra_data[4:6])[0]

    # Calculate compressed data size
    # Total block size = BSIZE + 1
    # Header size = 10 + 2 + XLEN  (header(10) + xlen(2) + extra(XLEN))
    # Trailer size = 8 (CRC32 + ISIZE)
    header_size = 10 + 2 + xlen
    cdata_size = (bsize + 1) - header_size - 8

    if cdata_size <= 0:  # pragma: no cover
        raise ValueError(f"Invalid compressed data size: {cdata_size}")

    # Read the compressed data
    compressed_data = handle.read(cdata_size)
    if len(compressed_data) != cdata_size:  # pragma: no cover
        raise ValueError(f"Truncated BGZF block data: expected {cdata_size}, got {len(compressed_data)}")

    # Read the trailer (CRC32 and ISIZE)
    trailer = handle.read(8)
    if len(trailer) != 8:  # pragma: no cover
        raise ValueError("Truncated BGZF trailer")

    crc32, isize = struct.unpack("<II", trailer)

    # Decompress the data
    try:
        # Use raw deflate decompression (negative window bits)
        data = zlib.decompress(compressed_data, -15)
    except zlib.error as e:  # pragma: no cover
        raise ValueError(f"Failed to decompress BGZF block: {e}")

    # Verify the uncompressed size
    if len(data) != isize:  # pragma: no cover
        raise ValueError(f"Uncompressed size mismatch: got {len(data)}, expected {isize}")

    # Verify the CRC32
    if zlib.crc32(data) & 0xFFFFFFFF != crc32:  # pragma: no cover
        raise ValueError("CRC32 mismatch in BGZF block")

    # Always convert to string using latin-1 encoding
    result = data.decode("latin-1")

    return bsize + 1, result


class BgzfWriter(typing.IO[str]):
    def __init__(
        self,
        filename: str | None = None,
        mode: str = "w",
        fileobj: typing.IO[bytes] | None = None,
        compresslevel: int = 6,
    ):
        if fileobj:
            assert filename is None
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():  # pragma: no cover
                raise ValueError("Must use write or append mode, not %r" % mode)
            if filename is None:  # pragma: no cover
                raise ValueError("Must give a filename if not passing a file handle")
            if "a" in mode.lower():  # pragma: no cover
                handle = open(filename, "ab")
            else:
                handle = open(filename, "wb")
        self._text: bool = "b" not in mode.lower()
        self._handle: typing.IO[bytes] = handle
        self._buffer: bytes = b""
        self.compresslevel: int = compresslevel
        self._filename = filename
        self._mode = mode
        self._closed = False

    def _write_block(self, block: bytes):
        # print("Saving %i bytes" % len(block))
        assert len(block) <= 65536
        # Giving a negative window bits means no gzip/zlib headers,
        # -15 used in samtools
        c = zlib.compressobj(self.compresslevel, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0)
        compressed = c.compress(block) + c.flush()
        del c
        assert len(compressed) < 65536, "TODO - Didn't compress enough, try less data in this block"
        crc = zlib.crc32(block)
        # Should cope with a mix of Python platforms...
        if crc < 0:  # pragma: no cover
            crc = struct.pack("<i", crc)
        else:
            crc = struct.pack("<I", crc)
        bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
        crc = struct.pack("<I", zlib.crc32(block) & 0xFFFFFFFF)
        uncompressed_length = struct.pack("<I", len(block))
        # Fixed 16 bytes,
        # gzip magic bytes (4) mod time (4),
        # gzip flag (1), os (1), extra length which is six (2),
        # sub field which is BC (2), sub field length of two (2),
        # Variable data,
        # 2 bytes: block length as BC sub field (2)
        # X bytes: the data
        # 8 bytes: crc (4), uncompressed data length (4)
        data = _bgzf_header + bsize + compressed + crc + uncompressed_length
        self._handle.write(data)

    def write(self, data: str) -> int:
        """Write string data to the BGZF file.

        Args:
            data: String data to write

        Returns:
            Number of characters written
        """
        if self._closed:  # pragma: no cover
            raise ValueError("I/O operation on closed file.")

        original_len = len(data)
        # Convert string to bytes using latin-1 encoding
        data_bytes = codecs.latin_1_encode(data)[0]

        # block_size = 2**16 = 65536
        data_len = len(data_bytes)
        if len(self._buffer) + data_len < 65536:
            # print("Cached %r" % data)
            self._buffer += data_bytes
        else:  # pragma: no cover
            # print("Got %r, writing out some data..." % data)
            self._buffer += data_bytes
            while len(self._buffer) >= 65536:
                self._write_block(self._buffer[:65536])
                self._buffer = self._buffer[65536:]

        return original_len

    def flush(self):
        while len(self._buffer) >= 65536:  # pragma: no cover
            self._write_block(self._buffer[:65535])
            self._buffer = self._buffer[65535:]
        self._write_block(self._buffer)
        self._buffer = b""
        self._handle.flush()

    def close(self) -> None:
        """Flush data, write 28 bytes BGZF EOF marker, and close BGZF file.
        samtools will look for a magic EOF marker, just a 28 byte empty BGZF
        block, and if it is missing warns the BAM file may be truncated. In
        addition to samtools writing this block, so too does bgzip - so this
        implementation does too.
        """
        if self._closed:  # pragma: no cover
            return

        if self._buffer:
            self.flush()
        self._handle.write(_bgzf_eof)
        self._handle.flush()
        self._handle.close()
        self._closed = True

    def tell(self) -> int:  # pragma: no cover
        """Returns a BGZF 64-bit virtual offset."""
        return make_virtual_offset(self._handle.tell(), len(self._buffer))

    def seekable(self) -> bool:  # pragma: no cover
        # Not seekable, but we do support tell...
        return False

    def isatty(self) -> bool:  # pragma: no cover
        """Return False as BGZF files are not TTY."""
        return False

    @property
    def closed(self) -> bool:  # pragma: no cover
        """Return True if the file is closed."""
        return self._closed

    @property
    def mode(self) -> str:  # pragma: no cover
        """Return the file mode."""
        return self._mode

    @property
    def name(self) -> str:  # pragma: no cover
        """Return the file name."""
        return self._filename or ""

    def readable(self) -> bool:  # pragma: no cover
        """Return False as this is a write-only file."""
        return False

    def writable(self) -> bool:  # pragma: no cover
        """Return True as this is a writable file."""
        return not self._closed

    def read(self, size: int = -1) -> str:  # pragma: no cover
        """Read operation not supported for write-only BGZF file."""
        raise OSError("not readable")

    def readline(self, size: int = -1) -> str:  # pragma: no cover
        """Readline operation not supported for write-only BGZF file."""
        raise OSError("not readable")

    def readlines(self, hint: int = -1) -> list[str]:  # pragma: no cover
        """Readlines operation not supported for write-only BGZF file."""
        raise OSError("not readable")

    def seek(self, offset: int, whence: int = 0) -> int:  # pragma: no cover
        """Seek operation not supported for BGZF files."""
        raise OSError("seek not supported on BGZF files")

    def truncate(self, size: int | None = None) -> int:  # pragma: no cover
        """Truncate operation not supported for BGZF files."""
        raise OSError("truncate not supported on BGZF files")

    def writelines(self, lines: Iterable[str]) -> None:
        """Write a list of strings to the file."""
        for line in lines:
            self.write(line)

    def fileno(self) -> int:
        return self._handle.fileno()

    def __enter__(self) -> "BgzfWriter":
        return self

    def __exit__(self, type_: type[BaseException] | None, value: BaseException | None, traceback: typing.Any) -> None:
        self.close()


class BgzfReader(typing.IO[str]):
    r"""BGZF reader, acts like a read only handle but seek/tell differ."""

    def __init__(
        self,
        filename: str | None = None,
        mode: str = "r",
        fileobj: typing.IO[bytes] | None = None,
        max_cache: int = 100,
    ):
        r"""Initialize the class for reading a BGZF file.

        You would typically use the top level ``bgzf.open(...)`` function
        which will call this class internally. Direct use is discouraged.

        Either the ``filename`` (string) or ``fileobj`` (input file object in
        binary mode) arguments must be supplied, but not both.

        Argument ``mode`` controls if the data will be returned as strings in
        text mode ("rt", "tr", or default "r"), or bytes binary mode ("rb"
        or "br"). The argument name matches the built-in ``open(...)`` and
        standard library ``gzip.open(...)`` function.

        If text mode is requested, in order to avoid multi-byte characters,
        this is hard coded to use the "latin1" encoding, and "\r" and "\n"
        are passed as is (without implementing universal new line mode). There
        is no ``encoding`` argument.

        If your data is in UTF-8 or any other incompatible encoding, you must
        use binary mode, and decode the appropriate fragments yourself.

        Argument ``max_cache`` controls the maximum number of BGZF blocks to
        cache in memory. Each can be up to 64kb thus the default of 100 blocks
        could take up to 6MB of RAM. This is important for efficient random
        access, a small value is fine for reading the file in one pass.
        """
        # TODO - Assuming we can seek, check for 28 bytes EOF empty block
        # and if missing warn about possible truncation (as in samtools)?
        if max_cache < 1:  # pragma: no cover
            raise ValueError("Use max_cache with a minimum of 1")
        # Must open the BGZF file in binary mode, but we may want to
        # treat the contents as either text or binary (unicode or
        # bytes under Python 3)
        if filename and fileobj:
            raise ValueError("Supply either filename or fileobj, not both")
        # Want to reject output modes like w, a, x, +
        if mode.lower() not in ("r", "tr", "rt", "rb", "br"):  # pragma: no cover
            raise ValueError("Must use a read mode like 'r' (default), 'rt', or 'rb' for binary")
        # If an open file was passed, make sure it was opened in binary mode.
        if fileobj:
            if fileobj.read(0) != b"":  # pragma: no cover
                raise ValueError("fileobj not opened in binary mode")
            handle = fileobj
        else:
            if filename is None:  # pragma: no cover
                raise ValueError("Must provide filename if fileobj is None")
            handle = open(filename, "rb")

        # Always read bytes from disk but return strings to the outside world
        self._newline = "\n"
        self._handle = handle
        self.max_cache = max_cache
        self._buffers: dict[int, tuple[str, int]] = {}
        self._block_start_offset: int = -1  # Force initial load
        self._block_raw_length: int = 0
        self._within_block_offset: int = 0
        self._buffer: str = ""
        self._load_block(0)  # Start at offset 0

    def _load_block(self, start_offset: int | None = None) -> None:
        if start_offset is None:
            # If the file is being read sequentially, then _handle.tell()
            # should be pointing at the start of the next block.
            # However, if seek has been used, we can't assume that.
            start_offset = self._block_start_offset + self._block_raw_length
        if start_offset == self._block_start_offset:  # pragma: no cover
            self._within_block_offset = 0
            return
        elif start_offset in self._buffers:
            # Already in cache
            self._buffer, self._block_raw_length = self._buffers[start_offset]
            self._within_block_offset = 0
            self._block_start_offset = start_offset
            return
        # Must hit the disk... first check cache limits,
        while len(self._buffers) >= self.max_cache:  # pragma: no cover
            # TODO - Implement LRU cache removal?
            self._buffers.popitem()
        # Now load the block
        handle = self._handle
        if start_offset is not None:
            handle.seek(start_offset)
        self._block_start_offset = handle.tell()
        try:
            block_size, self._buffer = _load_bgzf_block(handle)
        except StopIteration:
            # EOF
            block_size = 0
            self._buffer = ""
        self._within_block_offset = 0
        self._block_raw_length = block_size
        # Finally save the block in our cache,
        self._buffers[self._block_start_offset] = self._buffer, block_size

    def tell(self):  # pragma: no cover
        """Return a 64-bit unsigned BGZF virtual offset."""
        if 0 < self._within_block_offset and self._within_block_offset == len(self._buffer):
            # Special case where we're right at the end of a (non empty) block.
            # For non-maximal blocks could give two possible virtual offsets,
            # but for a maximal block can't use 65536 as the within block
            # offset. Therefore for consistency, use the next block and a
            # within block offset of zero.
            return (self._block_start_offset + self._block_raw_length) << 16
        else:
            # return make_virtual_offset(self._block_start_offset,
            #                           self._within_block_offset)
            # TODO - Include bounds checking as in make_virtual_offset?
            return (self._block_start_offset << 16) | self._within_block_offset

    def seek(self, virtual_offset: int, whence: int = 0) -> int:
        """Seek to a 64-bit unsigned BGZF virtual offset."""
        # Do this inline to avoid a function call,
        # start_offset, within_block = split_virtual_offset(virtual_offset)
        start_offset = virtual_offset >> 16
        within_block = virtual_offset ^ (start_offset << 16)
        if start_offset != self._block_start_offset:
            # Don't need to load the block if already there
            # (this avoids a function call since _load_block would do nothing)
            self._load_block(start_offset)
            if start_offset != self._block_start_offset:  # pragma: no cover
                raise ValueError("start_offset not loaded correctly")
        if within_block > len(self._buffer):
            if not (within_block == 0 and len(self._buffer) == 0):  # pragma: no cover
                raise ValueError("Within offset %i but block size only %i" % (within_block, len(self._buffer)))
        self._within_block_offset = within_block
        # assert virtual_offset == self.tell(), \
        #    "Did seek to %i (%i, %i), but tell says %i (%i, %i)" \
        #    % (virtual_offset, start_offset, within_block,
        #       self.tell(), self._block_start_offset,
        #       self._within_block_offset)
        return virtual_offset

    def read(self, size: int = -1) -> str:
        """Read method for the BGZF module."""
        if size < 0:  # pragma: no cover
            raise NotImplementedError("Don't be greedy, that could be massive!")

        result = ""
        while size and self._block_raw_length:
            if self._within_block_offset + size <= len(self._buffer):
                # This may leave us right at the end of a block
                # (lazy loading, don't load the next block unless we have too)
                data = self._buffer[self._within_block_offset : self._within_block_offset + size]
                self._within_block_offset += size
                if not data:  # pragma: no cover
                    raise ValueError("Must be at least 1 byte")
                result += data
                break
            else:  # pragma: no cover
                data = self._buffer[self._within_block_offset :]
                size -= len(data)
                self._load_block()  # will reset offsets
                result += data

        return result

    def readline(self, size: int = -1) -> str:
        """Read a single line for the BGZF file."""
        result = ""
        while self._block_raw_length:
            i = self._buffer.find(self._newline, self._within_block_offset)
            # Three cases to consider,
            if i == -1:  # pragma: no cover
                # No newline, need to read in more data
                data = self._buffer[self._within_block_offset :]
                self._load_block()  # will reset offsets
                result += data
            elif i + 1 == len(self._buffer):
                # Found new line, but right at end of block (SPECIAL)
                data = self._buffer[self._within_block_offset :]
                # Must now load the next block to ensure tell() works
                self._load_block()  # will reset offsets
                if not data:  # pragma: no cover
                    raise ValueError("Must be at least 1 byte")
                result += data
                break
            else:
                # Found new line, not at end of block (easy case, no IO)
                data = self._buffer[self._within_block_offset : i + 1]
                self._within_block_offset = i + 1
                # assert data.endswith(self._newline)
                result += data
                break

        return result

    def __next__(self) -> str:
        """Return the next line."""
        line = self.readline()
        if not line:
            raise StopIteration
        return line

    def __iter__(self) -> "BgzfReader":  # pragma: no cover
        """Iterate over the lines in the BGZF file."""
        return self

    def close(self) -> None:
        """Close BGZF file."""
        self._handle.close()
        self._buffer = ""
        self._block_start_offset = 0
        self._buffers = {}

    def seekable(self) -> bool:  # pragma: no cover
        """Return True indicating the BGZF supports random access."""
        return True

    def isatty(self) -> bool:  # pragma: no cover
        """Return True if connected to a TTY device."""
        return False

    def readable(self) -> bool:  # pragma: no cover
        """Return True indicating the BGZF file is readable."""
        return True

    def writable(self) -> bool:  # pragma: no cover
        """Return False indicating the BGZF file is not writable."""
        return False

    def fileno(self) -> int:  # pragma: no cover
        """Return integer file descriptor."""
        return self._handle.fileno()

    @property
    def closed(self) -> bool:  # pragma: no cover
        """Return True if the file is closed."""
        return self._handle.closed

    @property
    def mode(self) -> str:  # pragma: no cover
        """Return the file mode."""
        return "r"

    @property
    def name(self) -> str:
        """Return the file name."""
        return getattr(self._handle, "name", "")

    def flush(self) -> None:  # pragma: no cover
        """Flush - no-op for read-only file."""
        pass

    def readlines(self, hint: int = -1) -> list[str]:
        """Read all lines from the file."""
        lines = []
        for line in self:
            lines.append(line)
            if hint > 0 and len(lines) >= hint:
                break
        return lines

    def writelines(self, lines: Iterable[str]) -> None:  # pragma: no cover
        """Write lines - not supported for read-only file."""
        raise OSError("not writable")

    def write(self, s: str) -> int:  # pragma: no cover
        """Write - not supported for read-only file."""
        raise OSError("not writable")

    def truncate(self, size: int | None = None) -> int:  # pragma: no cover
        """Truncate - not supported for read-only file."""
        raise OSError("not writable")

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Close a file with WITH statement."""
        self.close()
