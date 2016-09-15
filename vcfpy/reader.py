# -*- coding: utf-8 -*-
"""Parsing of VCF files from ``file``-like objects
"""

import gzip

from . import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# TODO: use context manager for making close-able
# TODO: allow some configuration to make warning into exceptions?


class VCFReader:
    """Class for parsing of files from ``file``-like objects

    Instead of using the constructor, use the class methods
    :py:meth:`~VCFReader.from_file` and
    :py:meth:`~VCFReader.from_path`.

    On construction, the header will be read from the file which can cause
    problems.  After construction, :py:class:`~VCFReader` can be used as
    an iterable of :py:class:`~pyvcf.record.VCFRecord`.

    :raises: :py:class:`~vcfpy.exceptions.InvalidHeaderException` in the case
        of problems reading the header
    """

    @classmethod
    def from_file(klass, stream, path=None):
        """Create new :py:class:`VCFReader` from file

        :param stream: ``file``-like object to read from
        :param path: optional string with path to store (for display only)
        """
        return VCFReader(stream=stream, path=path)

    @classmethod
    def from_path(klass, path):
        """Create new :py:class:`VCFReader` from path

        :param path: the path to load from (converted to ``str`` for
            compatibility with ``path.py``)
        """
        path = str(path)
        if path.endswith('.gz'):
            f = gzip.open(path, 'rt')
        else:
            f = open(path, 'rt')
        return klass.from_file(stream=f, path=path)

    def __init__(self, stream, path=None):
        #: stream (``file``-like object) to read from
        self.stream = stream
        #: optional ``str`` with the path to the stream
        self.path = path
        #: the parser to use
        self.parser = parser.VCFParser(stream)
        #: the VCFHeader
        self.header = self.parser.parse_header()
        #: the :py:class:`pyvcf.header.SamplesInfos` object with the sample
        #: name information
        self.samples = self.header.samples

    def jump_to(self, chrom, begin, end):
        """Jump to the start position of the given chromosomal position
        and limit iteration to the end position

        :param str chrom: name of the chromosome to jump to
        :param int begin: 0-based begin position (inclusive)
        :param int end: 0-based end position (exclusive)
        """
        # TODO: this remains as an exercise for later...
        raise NotImplementedError('Implement me!')

    def __iter__(self):
        return self

    def __next__(self):
        """Return next object from file

        :returns:
        :raises: ``vcfpy.exceptions.InvalidRecordException`` in the case of
            problems reading the record
        :raises: ``StopException`` if at end
        """
        result = self.parser.parse_next_record()
        if result is None:
            raise StopException
        else:
            return result
