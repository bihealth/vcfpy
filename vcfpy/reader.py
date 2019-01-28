# -*- coding: utf-8 -*-
"""Parsing of VCF files from ``file``-like objects
"""

import gzip
import os

import pysam

from . import parser

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class Reader:
    """Class for parsing of files from ``file``-like objects

    Instead of using the constructor, use the class methods
    :py:meth:`~Reader.from_stream` and
    :py:meth:`~Reader.from_path`.

    On construction, the header will be read from the file which can cause
    problems.  After construction, :py:class:`~Reader` can be used as
    an iterable of :py:class:`~vcfpy.record.Record`.

    :raises: :py:class:`~vcfpy.exceptions.InvalidHeaderException` in the case
        of problems reading the header

    .. note::
        It is important to note that the ``header`` member is used during
        the parsing of the file.  **If you need a modified version then
        create a copy, e.g., using :py:method:`~vcfpy.header.Header.copy`**.

    .. note::
        If you use the ``parsed_samples`` feature and you write out
        records then you must not change the ``FORMAT`` of the record.
    """

    @classmethod
    def from_stream(
        klass, stream, path=None, tabix_path=None, record_checks=None, parsed_samples=None
    ):
        """Create new :py:class:`Reader` from file

        .. note::
            If you use the ``parsed_samples`` feature and you write out
            records then you must not change the ``FORMAT`` of the record.

        :param stream: ``file``-like object to read from
        :param path: optional string with path to store (for display only)
        :param list record_checks: record checks to perform, can contain
            'INFO' and 'FORMAT'
        :param list parsed_samples: ``list`` of ``str`` values with names of
            samples to parse call information for (for speedup); leave to
            ``None`` for ignoring
        """
        record_checks = record_checks or []
        if tabix_path and not path:
            raise ValueError("Must give path if tabix_path is given")
        return Reader(
            stream=stream,
            path=path,
            tabix_path=tabix_path,
            record_checks=record_checks,
            parsed_samples=parsed_samples,
        )

    @classmethod
    def from_path(klass, path, tabix_path=None, record_checks=None, parsed_samples=None):
        """Create new :py:class:`Reader` from path

        .. note::
            If you use the ``parsed_samples`` feature and you write out
            records then you must not change the ``FORMAT`` of the record.

        :param path: the path to load from (converted to ``str`` for
            compatibility with ``path.py``)
        :param tabix_path: optional string with path to TBI index,
            automatic inferral from ``path`` will be tried on the fly
            if not given
        :param list record_checks: record checks to perform, can contain
            'INFO' and 'FORMAT'
        """
        record_checks = record_checks or []
        path = str(path)
        if path.endswith(".gz"):
            f = gzip.open(path, "rt")
            if not tabix_path:
                tabix_path = path + ".tbi"
                if not os.path.exists(tabix_path):
                    tabix_path = None  # guessing path failed
        else:
            f = open(path, "rt")
        return klass.from_stream(
            stream=f,
            path=path,
            tabix_path=tabix_path,
            record_checks=record_checks,
            parsed_samples=parsed_samples,
        )

    def __init__(self, stream, path=None, tabix_path=None, record_checks=None, parsed_samples=None):
        #: stream (``file``-like object) to read from
        self.stream = stream
        #: optional ``str`` with the path to the stream
        self.path = path
        #: optional ``str`` with path to tabix file
        self.tabix_path = tabix_path
        #: checks to perform on records, can contain 'FORMAT' and 'INFO'
        self.record_checks = tuple(record_checks or [])
        #: if set, list of samples to parse for
        self.parsed_samples = parsed_samples
        #: the ``pysam.TabixFile`` used for reading from index bgzip-ed VCF;
        #: constructed on the fly
        self.tabix_file = None
        # the iterator through the Tabix file to use
        self.tabix_iter = None
        #: the parser to use
        self.parser = parser.Parser(stream, self.path, self.record_checks)
        #: the Header
        self.header = self.parser.parse_header(parsed_samples)

    def fetch(self, chrom_or_region, begin=None, end=None):
        """Jump to the start position of the given chromosomal position
        and limit iteration to the end position

        :param str chrom_or_region: name of the chromosome to jump to if
            begin and end are given and a samtools region string otherwise
            (e.g. "chr1:123,456-123,900").
        :param int begin: 0-based begin position (inclusive)
        :param int end: 0-based end position (exclusive)
        """
        if begin is not None and end is None:
            raise ValueError("begin and end must both be None or neither")
        # close tabix file if any and is open
        if self.tabix_file and not self.tabix_file.closed:
            self.tabix_file.close()
        # open tabix file if not yet open
        if not self.tabix_file or self.tabix_file.closed:
            self.tabix_file = pysam.TabixFile(filename=self.path, index=self.tabix_path)
        # jump to the next position
        if begin is None:
            self.tabix_iter = self.tabix_file.fetch(region=chrom_or_region)
        else:
            self.tabix_iter = self.tabix_file.fetch(reference=chrom_or_region, start=begin, end=end)
        return self

    def close(self):
        """Close underlying stream"""
        if self.tabix_file and not self.tabix_file.closed:
            self.tabix_file.close()
        if self.stream:
            self.stream.close()

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        """Return next object from file

        :returns:
        :raises: ``vcfpy.exceptions.InvalidRecordException`` in the case of
            problems reading the record
        :raises: ``StopException`` if at end
        """
        if self.tabix_iter:
            return self.parser.parse_line(str(next(self.tabix_iter)))
        else:
            result = self.parser.parse_next_record()
            if result is None:
                raise StopIteration()
            else:
                return result
