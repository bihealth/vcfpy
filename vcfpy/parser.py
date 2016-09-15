# -*- coding: utf-8 -*-
"""Parsing of VCF files from ``str``
"""

import ast
import collections
import itertools

from . import header
from . import exceptions

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# TODO: check whether the str.split() calls are a bottleneck, could write
#       faster parsers based on scanning (VCF follows a regular grammar
#       after all)


# expected "#CHROM" header prefix
REQUIRE_SAMPLE_HEADER = (
    '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
    'FORMAT')

#: Supported VCF versions, a warning will be issued otherwise
SUPPORTED_VCF_VERSIONS = (
    'VCFv4.0', 'VCFv4.1', 'VCFv4.2', 'VCFv4.3')


def split_quoted_string(s, delim=',', quote='"'):
    """Split string ``s`` at delimiter, correctly interpreting quotes"""
    # collect positions
    begins, ends = [0], []
    # run state automaton
    NORMAL, QUOTED, ESCAPED, DELIM = 0, 1, 2, 3
    state = NORMAL
    for pos, c in enumerate(s):
        if state == NORMAL:
            if c == delim:
                ends.append(pos)
                state = DELIM
            elif c == quote:
                state = QUOTED
            else:
                pass  # noop
        elif state == QUOTED:
            if c == '\\':
                state = ESCAPED
            elif c == quote:
                state = NORMAL
            else:
                pass  # noop
        elif state == DELIM:
            begins.append(pos)
            state = NORMAL
        else:  # state == ESCAPED
            state = QUOTED
    ends.append(len(s))
    assert len(begins) == len(ends)
    # Build resulting list
    return [s[start:end] for start, end in zip(begins, ends)]


def parse_mapping(value):
    """Parse the given VCF header line mapping

    Such a mapping consists of "key=value" pairs, separated by commas and
    wrapped into angular brackets ("<...>").  Strings are usually quoted,
    for certain known keys, exceptions are made, depending on the tag key.
    this, however, only gets important when serializing.

    :raises: :py:class:`pyvcf.exceptions.InvalidHeaderException` if
        there was a problem parsing the file
    """
    if not value.startswith('<') or not value.endswith('>'):
        raise exceptions.InvalidHeaderException(
            'Header mapping value was not wrapped in angular brackets')
    # split the comma-separated list into pairs, ignoring commas in quotes
    pairs = split_quoted_string(value[1:-1], delim=',', quote='"')
    # split these pairs into key/value pairs, converting flags to mappings
    # to True
    key_values = []
    for pair in pairs:
        if '=' in pair:
            key, value = pair.split('=', 1)
            if value.startswith('"') and value.endswith('"'):
                value = ast.literal_eval(value)
        else:
            key, value = pair, True
        key_values.append((key, value))
    # return completely parsed mapping as OrderedDict
    return collections.OrderedDict(key_values)


class VCFHeaderLineParserBase:
    """Parse into appropriate VCFHeaderLine"""

    def parse_key_value(self, key, value):
        """Parse the key/value pair

        :param str key: the key to use in parsing
        :param str value: the value to parse
        :returns: :py:class:`vcfpy.header.VCFHeaderLine` object
        """
        raise NotImplementedError('Must be overridden')


class StupidVCFHeaderLineParser(VCFHeaderLineParserBase):
    """Parse into VCFHeaderLine (no particular structure)"""

    def parse_key_value(self, key, value):
        return header.VCFHeaderLine(key, value)


class MappingVCFHeaderLineParser(VCFHeaderLineParserBase):
    """Parse into VCFHeaderLine (no particular structure)"""

    def __init__(self, line_class):
        """Initialize the parser"""
        super().__init__()
        #: the class to use for the VCF header line
        self.line_class = line_class

    def parse_key_value(self, key, value):
        return self.line_class(key, value, parse_mapping(value))


# Parsers to use for each VCF header type (given left of '=')
HEADER_PARSERS = {
    'FILTER': MappingVCFHeaderLineParser(header.VCFFilterHeaderLine),
    'FORMAT': MappingVCFHeaderLineParser(header.VCFFormatHeaderLine),
    'INFO': MappingVCFHeaderLineParser(header.VCFInfoHeaderLine),
    'contig': MappingVCFHeaderLineParser(header.VCFContigHeaderLine),
    '__default__': StupidVCFHeaderLineParser(),  # fallback
    # 'ALT': None,
    # 'assembly': None,
    # 'META': None,
    # 'SAMPLE': None,
    # 'PEDIGREE': None,
    # 'pedigreeDB': None,
}


class VCFHeaderParser:
    """Helper class for parsing a VCF header
    """

    def __init__(self, sub_parsers):
        self.sub_parsers = HEADER_PARSERS

    def parse_line(self, line):
        """Parse VCF header ``line`` (trailing '\r\n' or '\n' is ignored)

        :param str line: ``str`` with line to parse
        :param dict sub_parsers: ``dict`` mapping header line types to
            appropriate parser objects
        :returns: appropriate :py:class:`VCFHeaderLine` parsed from ``line``
        :raises: :py:class:`pyvcf.exceptions.InvalidHeaderException` if
            there was a problem parsing the file
        """
        if not line or not line.startswith('##'):
            raise exceptions.InvalidHeaderException(
                'Invalid VCF header line (must start with "##") {}'.format(
                    line))
        if '=' not in line:
            raise exceptions.InvalidHeaderException(
                'Invalid VCF header line (must contain "=") {}'.format(line))
        line = line[len('##'):].rstrip()  # trim '^##' and trailing whitespace
        # split key/value pair at "="
        key, value = line.split('=', 1)
        sub_parser = self.sub_parsers.get(key, self.sub_parsers['__default__'])
        return sub_parser.parse_key_value(key, value)


class VCFParser:
    """Class for line-wise parsing of VCF files

    In most cases, you want to use :py:class:`vcfpy.reader.VCFReader` instead.

    :param stream: ``file``-like object to read from
    :param str path: path the VCF is parsed from, for display purposes
        only, optional
    """

    def __init__(self, stream, path=None):
        self.stream = stream
        self.path = path
        #: header, once it has been read
        self.header = None
        # the currently read line
        self._line = stream.readline()  # trailing '\n'
        #: :py:class:`vcfpy.header.SamplesInfos` with sample information;
        #: set on parsing the header
        self.samples = None

    def _read_next_line(self):
        """Read and return next line, copy to self._line"""
        self._line = self.stream.readline()
        return self._line

    def parse_header(self):
        """Read and parse :py:class:`vcfpy.header.VCFHeader` from file, set
        into ``self.header`` and return it

        :returns: ``vcfpy.header.VCFHeader``
        :raises: ``vcfpy.exceptions.InvalidHeaderException`` in the case of
            problems reading the header
        """
        # parse header lines
        sub_parser = VCFHeaderParser(HEADER_PARSERS)
        header_lines = []
        while self._line and self._line.startswith('##'):
            header_lines.append(sub_parser.parse_line(self._line))
            self._read_next_line()
        # check first header line to be '##fileformat='
        self._check_header_lines(header_lines)  # raises InvalidHeaderException
        # parse sample info line
        if not self._line or not self._line.startswith('#CHROM'):
            raise exceptions.IncorrectVCFFormat(
                'Missing line starting with "#CHROM"')
        arr = self._line.rstrip().split('\t')
        if tuple(arr[:len(REQUIRE_SAMPLE_HEADER)]) != REQUIRE_SAMPLE_HEADER:
            raise exceptions.IncorrectVCFFormat(
                'Sample header line (starting with "#CHROM") does not '
                'start with required prefix {}'.format(
                    '\t'.join(REQUIRE_SAMPLE_HEADER)))
        self.samples = header.SamplesInfos(arr[len(REQUIRE_SAMPLE_HEADER):])
        # construct VCFHeader object
        self.header = header.VCFHeader(header_lines, self.samples)
        # read next line, must not be header
        self._read_next_line()
        if self._line and self._line.startswith('#'):
            raise exceptions.IncorrectVCFFormat(
                'Expecting non-header line or EOF after "#CHROM" line')
        return self.header

    def _check_header_lines(self, header_lines):
        """Check header lines, in particular for starting file "##fileformat"

        :raises: ``vcfpy.exceptions.InvalidHeaderException`` in the case of
            problems reading the header
        """
        if not header_lines:
            raise exceptions.InvalidHeaderException(
                'The VCF file did not contain any header lines!')
        first = header_lines[0]
        if first.key != 'fileformat':
            raise exceptions.InvalidHeaderException(
                'The VCF file did not start with ##fileformat')
        if first.value not in SUPPORTED_VCF_VERSIONS:
            print(('[vcfpy] WARNING: The VCF version {} is not known, '
                   'going on regardlessly').format(first.value),
                  file=sys.stderr)

    def parse_next_record(self):
        """Read, parse and return next :py:class:`vcfpy.record.VCFRecord`

        :returns: next VCF record or ``None`` if at end
        :raises: ``vcfpy.exceptions.InvalidRecordException`` in the case of
            problems reading the record
        """
        # XXX
        return None
