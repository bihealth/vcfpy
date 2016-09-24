# -*- coding: utf-8 -*-
"""Parsing of VCF files from ``file``-like objects
"""

import gzip
import os

import pysam

from . import parser

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


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
    """

    @classmethod
    def from_stream(klass, stream, path=None, tabix_path=None,
                    record_checks=[]):
        """Create new :py:class:`Reader` from file

        :param stream: ``file``-like object to read from
        :param path: optional string with path to store (for display only)
        :param list record_checks: record checks to perform, can contain
            'INFO' and 'FORMAT'
        """
        if tabix_path and not path:
            raise ValueError('Must give path if tabix_path is given')
        return Reader(stream=stream, path=path, tabix_path=tabix_path,
                      record_checks=record_checks)

    @classmethod
    def from_path(klass, path, tabix_path=None, record_checks=[]):
        """Create new :py:class:`Reader` from path

        :param path: the path to load from (converted to ``str`` for
            compatibility with ``path.py``)
        :param tabix_path: optional string with path to TBI index,
            automatic inferral from ``path`` will be tried on the fly
            if not given
        :param list record_checks: record checks to perform, can contain
            'INFO' and 'FORMAT'
        """
        path = str(path)
        if path.endswith('.gz'):
            f = gzip.open(path, 'rt')
            if not tabix_path:
                tabix_path = path + '.tbi'
                if not os.path.exists(tabix_path):
                    tabix_path = None  # guessing path failed
        else:
            f = open(path, 'rt')
        return klass.from_stream(stream=f, path=path, tabix_path=tabix_path,
                                 record_checks=record_checks)

    def __init__(self, stream, path=None, tabix_path=None,
                 record_checks=[]):
        #: stream (``file``-like object) to read from
        self.stream = stream
        #: optional ``str`` with the path to the stream
        self.path = path
        #: optional ``str`` with path to tabix file
        self.tabix_path = tabix_path
        #: checks to perform on records, can contain 'FORMAT' and 'INFO'
        self.record_checks = tuple(record_checks)
        #: the ``pysam.TabixFile`` used for reading from index bgzip-ed VCF;
        #: constructed on the fly
        self.tabix_file = None
        # the iterator through the Tabix file to use
        self.tabix_iter = None
        #: the parser to use
        self.parser = parser.Parser(stream, self.path, self.record_checks)
        #: the Header
        self.header = self.parser.parse_header()
        #: the :py:class:`vcfpy.header.SamplesInfos` object with the sample
        #: name information
        self.samples = self.header.samples

    def fetch(self, chrom, begin, end):
        """Jump to the start position of the given chromosomal position
        and limit iteration to the end position

        :param str chrom: name of the chromosome to jump to
        :param int begin: 0-based begin position (inclusive)
        :param int end: 0-based end position (exclusive)
        """
        # close tabix file if any and is open
        if self.tabix_file and not self.tabix_file.closed:
            self.tabix_file.close()
        # open tabix file if not yet open
        if not self.tabix_file or self.tabix_file.closed:
            self.tabix_file = pysam.TabixFile(
                filename=self.path, index=self.tabix_path)
        # jump to the next position
        self.tabix_iter = self.tabix_file.fetch(chrom, begin, end)
        return self

    def close(self):
        """Close underlying stream"""
        if self.tabix_file and not self.tabix_file.closed:
            self.tabix_file.close()
        if self.stream:
            self.stream.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
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
