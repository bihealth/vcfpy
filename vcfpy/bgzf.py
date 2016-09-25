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
import zlib


# For Python 2 can just use: _bgzf_magic = '\x1f\x8b\x08\x04'
# but need to use bytes on Python 3
_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42"
                b"\x43\x02\x00")
_bgzf_eof = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
             b"\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00")
_bytes_BC = b"BC"


def make_virtual_offset(block_start_offset, within_block_offset):
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
        raise ValueError("Require 0 <= within_block_offset < 2**16, got %i" %
                         within_block_offset)
    if block_start_offset < 0 or block_start_offset >= 281474976710656:
        raise ValueError("Require 0 <= block_start_offset < 2**48, got %i" %
                         block_start_offset)
    return (block_start_offset << 16) | within_block_offset


class BgzfWriter(object):

    def __init__(self, filename=None, mode="w", fileobj=None, compresslevel=6):
        if fileobj:
            assert filename is None
            handle = fileobj
        else:
            if "w" not in mode.lower() and "a" not in mode.lower():
                raise ValueError(
                    "Must use write or append mode, not %r" % mode)
            if "a" in mode.lower():
                handle = open(filename, "ab")
            else:
                handle = open(filename, "wb")
        self._text = "b" not in mode.lower()
        self._handle = handle
        self._buffer = b""
        self.compresslevel = compresslevel

    def _write_block(self, block):
        # print("Saving %i bytes" % len(block))
        assert len(block) <= 65536
        # Giving a negative window bits means no gzip/zlib headers,
        # -15 used in samtools
        c = zlib.compressobj(self.compresslevel,
                             zlib.DEFLATED,
                             -15,
                             zlib.DEF_MEM_LEVEL,
                             0)
        compressed = c.compress(block) + c.flush()
        del c
        assert len(compressed) < 65536, \
            "TODO - Didn't compress enough, try less data in this block"
        crc = zlib.crc32(block)
        # Should cope with a mix of Python platforms...
        if crc < 0:
            crc = struct.pack("<i", crc)
        else:
            crc = struct.pack("<I", crc)
        bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
        crc = struct.pack("<I", zlib.crc32(block) & 0xffffffff)
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

    def write(self, data):
        # TODO - Check bytes vs unicode
        if isinstance(data, str):
            data = codecs.latin_1_encode(data)[0]
        # block_size = 2**16 = 65536
        data_len = len(data)
        if len(self._buffer) + data_len < 65536:
            # print("Cached %r" % data)
            self._buffer += data
            return
        else:
            # print("Got %r, writing out some data..." % data)
            self._buffer += data
            while len(self._buffer) >= 65536:
                self._write_block(self._buffer[:65536])
                self._buffer = self._buffer[65536:]

    def flush(self):
        while len(self._buffer) >= 65536:
            self._write_block(self._buffer[:65535])
            self._buffer = self._buffer[65535:]
        self._write_block(self._buffer)
        self._buffer = b""
        self._handle.flush()

    def close(self):
        """Flush data, write 28 bytes BGZF EOF marker, and close BGZF file.
        samtools will look for a magic EOF marker, just a 28 byte empty BGZF
        block, and if it is missing warns the BAM file may be truncated. In
        addition to samtools writing this block, so too does bgzip - so this
        implementation does too.
        """
        if self._buffer:
            self.flush()
        self._handle.write(_bgzf_eof)
        self._handle.flush()
        self._handle.close()

    def tell(self):
        """Returns a BGZF 64-bit virtual offset."""
        return make_virtual_offset(self._handle.tell(), len(self._buffer))

    def seekable(self):
        # Not seekable, but we do support tell...
        return False

    @classmethod
    def isatty(klass):
        return False

    def fileno(self):
        return self._handle.fileno()

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()
