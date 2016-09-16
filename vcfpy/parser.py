# -*- coding: utf-8 -*-
"""Parsing of VCF files from ``str``
"""

import ast
import collections
import itertools
import sys

from . import header
from . import record
from . import exceptions

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


# expected "#CHROM" header prefix when there are samples
REQUIRE_SAMPLE_HEADER = (
    '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
# expected "#CHROM" header prefix when there are no samples
REQUIRE_NO_SAMPLE_HEADER = (
    '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

#: Supported VCF versions, a warning will be issued otherwise
SUPPORTED_VCF_VERSIONS = (
    'VCFv4.0', 'VCFv4.1', 'VCFv4.2', 'VCFv4.3')


def _warn(msg):
    """Print warning message"""
    print('[vcfpy] WARNING: {}'.format(msg), file=sys.stderr)


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

    :raises: :py:class:`vcfpy.exceptions.InvalidHeaderException` if
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
            key = key.strip()  # XXX lenient parsing
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


# Field value converters
_CONVERTERS = {
    'Integer': int,
    'Float': float,
    'Flag': lambda x: True,
    'Character': str,
    'String': str,
}


def convert_field_value(key, type_, value):
    """Convert atomic field value according to the type"""
    if value == '.':
        return None
    else:
        return _CONVERTERS[type_](value)


def parse_field_value(key, field_info, value):
    """Parse ``value`` according to ``field_info``
    """
    # XXX lenient parsing, ignoring if counts differ from field_info
    if field_info.type == 'Flag':
        assert value is True
        return True
    elif field_info.number == 1:
        return convert_field_value(key, field_info.type, value)
    else:
        if value == '.':
            return []
        else:
            return [convert_field_value(key, field_info.type, x)
                    for x in value.split(',')]


def process_alt(header, ref, alt_str):
    """Process alternative value using VCFHeader in ``header``"""
    # By its nature, this function contains a large number of case distinctions
    if ']' in alt_str or '[' in alt_str:
        return record.BreakEnd(record.BND, alt_str)
    elif alt_str.startswith('<') and alt_str.endswith('>'):
        inner = alt_str[1:-1]
        if any([inner.startswith(code) for code in record.SV_CODES]):
            return record.SV(record.SV, alt_str)
        else:
            return record.SymbolicAllele(record.SYMBOLIC, alt_str)
    else:  # substitution
        if len(ref) == len(alt_str):
            if len(ref) == 1:
                return record.Substitution(record.SNV, alt_str)
            else:
                return record.Substitution(record.MNV, alt_str)
        elif len(ref) > len(alt_str):
            if len(alt_str) == 0:
                raise exceptions.InvalidRecordException(
                    'Invalid VCF, empty ALT')
            elif len(alt_str) == 1:
                if ref[0] == alt_str[0]:
                    return record.Substitution(record.DEL, alt_str)
                else:
                    return record.Substitution(record.INDEL, alt_str)
            else:
                return record.Substitution(record.INDEL, alt_str)
        else:  # len(ref) < len(alt_str):
            if len(ref) == 0:
                raise exceptions.InvalidRecordException(
                    'Invalid VCF, empty REF')
            elif len(ref) == 1:
                if ref[0] == alt_str[0]:
                    return record.Substitution(record.INS, alt_str)
                else:
                    return record.Substitution(record.INDEL, alt_str)
            else:
                return record.Substitution(record.INDEL, alt_str)


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
        :raises: :py:class:`vcfpy.exceptions.InvalidHeaderException` if
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
        key = key.strip()  # XXX lenient parsing
        sub_parser = self.sub_parsers.get(key, self.sub_parsers['__default__'])
        return sub_parser.parse_key_value(key, value)


class VCFRecordParser:
    """Helper class for parsing VCF records"""

    def __init__(self, header, samples):
        #: VCFHeader with the meta information
        self.header = header
        #: SamplesInfos with sample information
        self.samples = samples
        # Expected number of fields
        if self.samples.names:
            self.expected_fields = 9 + len(self.samples.names)
        else:
            self.expected_fields = 8

    def parse_line(self, line_str):
        """Parse line from file (including trailing line break) and return
        resulting VCFRecord
        """
        line_str = line_str.rstrip()
        if not line_str:
            return None  # empty line, EOF
        arr = self._split_line(line_str)
        # CHROM
        chrom = arr[0]
        # POS
        pos = int(arr[1])
        # IDS
        if arr[2] == '.':
            ids = []
        else:
            ids = arr[2].split(';')
        # REF
        ref = arr[3]
        # ALT
        if arr[4] == '.':
            alts = []
        else:
            alts = list(map(lambda x: process_alt(self.header, ref, x),
                            arr[4].split(',')))
        # QUAL
        if arr[5] == '.':
            qual = None
        else:
            try:
                qual = int(arr[5])
            except ValueError:  # try as float
                qual = float(arr[5])
        # FILTER
        if arr[6] == '.':
            filt = []
        else:
            filt = arr[6].split(';')
        # INFO
        info = self._parse_info(arr[7])
        if not self.samples.names:
            format, calls = [], []
        else:
            # FORMAT
            format = arr[8].split(':')
            # per-sample calls
            calls = [record.Call(sample, data) for sample, data in
                    zip(self.samples.names,
                        self._parse_calls_data(format, arr[9:]))]
        return record.Record(
            chrom, pos, ids, ref, alts, qual, filt, info, format, calls)

    def _split_line(self, line_str):
        """Split line and check number of columns"""
        arr = line_str.rstrip().split('\t')
        if len(arr) != self.expected_fields:
            raise exceptions.InvalidRecordException(
                ('The line contains an invalid number of fields. Was '
                 '{} but expected {}\n{}'.format(
                     len(arr), 9 + len(self.samples.names),
                     line_str)))
        return arr

    def _parse_info(self, info_str):
        """Parse INFO column from string"""
        result = collections.OrderedDict()
        if info_str == '.':
            return result
        # The standard is very nice to parsers, we can simply split at
        # semicolon characters, although I (Manuel) don't know how strict
        # programs follow this
        for entry in info_str.split(';'):
            if '=' not in entry:  # flag
                result[entry] = parse_field_value(
                    entry, self.header.get_info_field_info(entry), True)
            else:
                key, value = entry.split('=', 1)
                key = key.strip()  # XXX lenient parsing
                result[key] = parse_field_value(
                    key, self.header.get_info_field_info(key), value)
        return result

    def _parse_calls_data(self, format, arr):
        """Parse genotype call information from arrays using format array

        :param list format: List of strings with format names
        :param list arr: List of strings with genotype information values
            for each sample
        """
        result = []
        for entry in arr:
            data = collections.OrderedDict()
            # The standard is very nice to parsers, we can simply split at
            # colon characters, although I (Manuel) don't know how strict
            # programs follow this
            for key, value in zip(format, entry.split(':')):
                data[key] = parse_field_value(
                    key, self.header.get_format_field_info(key), value)
            result.append(data)
        return result


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
        # helper for parsing the records
        self._record_parser = None

    def _read_next_line(self):
        """Read next line store in self._line and return old one"""
        prev_line = self._line
        self._line = self.stream.readline()
        return prev_line

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
        # check for space before INFO
        line = self._line.rstrip()
        pos = line.find('FORMAT') if ('FORMAT' in line) else line.find('INFO')
        if pos == -1:
            raise exceptions.IncorrectVCFFormat(
                'Ill-formatted line starting with "#CHROM"')
        if ' ' in line[:pos]:
            _warn('Found space in #CHROM line, splitting at whitespace '
                  'instead of tab; this VCF file is ill-formatted')
            arr = self._line.rstrip().split()
        else:
            arr = self._line.rstrip().split('\t')

        if len(arr) <= len(REQUIRE_NO_SAMPLE_HEADER):
            if tuple(arr) != REQUIRE_NO_SAMPLE_HEADER:
                raise exceptions.IncorrectVCFFormat(
                    'Sample header line indicates no sample but does not '
                    'equal required prefix {}'.format(
                        '\t'.join(REQUIRE_NO_SAMPLE_HEADER)))
        elif tuple(arr[:len(REQUIRE_SAMPLE_HEADER)]) != REQUIRE_SAMPLE_HEADER:
            raise exceptions.IncorrectVCFFormat(
                'Sample header line (starting with "#CHROM") does not '
                'start with required prefix {}'.format(
                    '\t'.join(REQUIRE_SAMPLE_HEADER)))
        self.samples = header.SamplesInfos(arr[len(REQUIRE_SAMPLE_HEADER):])
        # construct VCFHeader object
        self.header = header.VCFHeader(header_lines, self.samples)
        # construct record parser
        self._record_parser = VCFRecordParser(self.header, self.samples)
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
        return self._record_parser.parse_line(self._read_next_line())
