# -*- coding: utf-8 -*-
"""Parsing of VCF files from ``str``"""

import ast
import functools
import io
import math
import pathlib
import re
import warnings
from typing import Any, Callable, Iterable, Literal, cast

from vcfpy import exceptions, header, record

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


# expected "#CHROM" header prefix when there are samples
REQUIRE_SAMPLE_HEADER = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
# expected "#CHROM" header prefix when there are no samples
REQUIRE_NO_SAMPLE_HEADER = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

#: Supported VCF versions, a warning will be issued otherwise
SUPPORTED_VCF_VERSIONS = ("VCFv4.0", "VCFv4.1", "VCFv4.2", "VCFv4.3")


class QuotedStringSplitter:
    """Helper class for splitting quoted strings

    Has support for interpreting quoting strings but also brackets.  Meant
    for splitting the VCF header line dicts
    """

    #: state constant for normal
    NORMAL = 0
    #: state constant for quoted
    QUOTED = 1
    #: state constant for delimiter
    ESCAPED = 2
    #: state constant for array
    ARRAY = 3
    #: state constant for delimiter
    DELIM = 4

    def __init__(self, delim: str = ",", quote: str = '"', brackets: str = "[]"):
        #: string delimiter
        self.delim = delim
        #: quote character
        self.quote = quote
        #: two-character string with opening and closing brackets
        assert len(brackets) == 2
        self.brackets = brackets

    def run(self, s: str) -> list[str]:
        """Split string ``s`` at delimiter, correctly interpreting quotes

        Further, interprets arrays wrapped in one level of ``[]``.  No
        recursive brackets are interpreted (as this would make the grammar
        non-regular and currently this complexity is not needed).  Currently,
        quoting inside of braces is not supported either.  This is just to
        support the example from VCF v4.3.
        """
        begins: list[int] = [0]
        ends: list[int] = []
        # transition table
        DISPATCH: dict[Literal[0, 1, 2, 3, 4], Callable[[str, int, list[int], list[int]], Literal[0, 1, 2, 3, 4]]] = {
            self.NORMAL: self._handle_normal,
            self.QUOTED: self._handle_quoted,
            self.ARRAY: self._handle_array,
            self.DELIM: self._handle_delim,
            self.ESCAPED: self._handle_escaped,
        }
        # run state automaton
        state: Literal[0, 1, 2, 3, 4] = self.NORMAL
        for pos, c in enumerate(s):
            state = DISPATCH[state](c, pos, begins, ends)
        ends.append(len(s))
        assert len(begins) == len(ends)
        # Build resulting list
        return [s[start:end] for start, end in zip(begins, ends, strict=False)]

    def _handle_normal(
        self, c: str, pos: int, begins: list[int], ends: list[int]
    ) -> Literal[0, 1, 2, 3, 4]:  # pylint: disable=W0613
        if c == self.delim:
            ends.append(pos)
            return self.DELIM
        elif c == self.quote:
            return self.QUOTED
        elif c == self.brackets[0]:
            return self.ARRAY
        else:
            return self.NORMAL

    def _handle_quoted(
        self, c: str, pos: int, begins: list[int], ends: list[int]
    ) -> Literal[0, 1, 2, 3, 4]:  # pylint: disable=W0613
        if c == "\\":
            return self.ESCAPED
        elif c == self.quote:
            return self.NORMAL
        else:
            return self.QUOTED

    def _handle_array(
        self, c: str, pos: int, begins: list[int], ends: list[int]
    ) -> Literal[0, 1, 2, 3, 4]:  # pylint: disable=W0613
        if c == self.brackets[1]:
            return self.NORMAL
        else:
            return self.ARRAY

    def _handle_delim(
        self, c: str, pos: int, begins: list[int], ends: list[int]
    ) -> Literal[0, 1, 2, 3, 4]:  # pylint: disable=W0613
        begins.append(pos)
        return self.NORMAL

    def _handle_escaped(
        self, c: str, pos: int, begins: list[int], ends: list[int]
    ) -> Literal[0, 1, 2, 3, 4]:  # pylint: disable=W0613
        return self.QUOTED


def split_quoted_string(s: str, delim: str = ",", quote: str = '"', brackets: str = "[]") -> list[str]:
    return QuotedStringSplitter(delim, quote, brackets).run(s)


def split_mapping(pair_str: str) -> tuple[str, str]:
    """Split the ``str`` in ``pair_str`` at ``'='``

    Warn if key needs to be stripped
    """
    orig_key, value = pair_str.split("=", 1)
    key = orig_key.strip()
    if key != orig_key:
        warnings.warn(
            "Mapping key {} has leading or trailing space".format(repr(orig_key)),
            exceptions.LeadingTrailingSpaceInKey,
        )
    return key, value


def parse_mapping(value: str) -> dict[str, bool | str | list[str]]:
    """Parse the given VCF header line mapping

    Such a mapping consists of "key=value" pairs, separated by commas and
    wrapped into angular brackets ("<...>").  Strings are usually quoted,
    for certain known keys, exceptions are made, depending on the tag key.
    this, however, only gets important when serializing.

    :raises: :py:class:`vcfpy.exceptions.InvalidHeaderException` if
        there was a problem parsing the file
    """
    if not value.startswith("<") or not value.endswith(">"):
        raise exceptions.InvalidHeaderException("Header mapping value was not wrapped in angular brackets")
    # split the comma-separated list into pairs, ignoring commas in quotes
    pairs = split_quoted_string(value[1:-1], delim=",", quote='"')
    # split these pairs into key/value pairs, converting flags to mappings
    # to True
    key_values: list[tuple[str, bool | str | list[str]]] = []
    for pair in pairs:
        value_: bool | str | list[str]
        if "=" in pair:
            key, value = split_mapping(pair)
            if value.startswith('"') and value.endswith('"'):
                value_ = ast.literal_eval(value)
            elif value.startswith("[") and value.endswith("]"):
                value_ = [v.strip() for v in value[1:-1].split(",")]
            else:
                value_ = value
        else:
            key, value_ = pair, True
        key_values.append((key, value_))
    # return completely parsed mapping as OrderedDict
    return dict(key_values)


class HeaderLineParserBase:
    """Parse into appropriate HeaderLine"""

    def parse_key_value(self, key: str, value: str) -> header.HeaderLine:
        """Parse the key/value pair

        :param str key: the key to use in parsing
        :param str value: the value to parse
        :returns: :py:class:`vcfpy.header.HeaderLine` object
        """
        raise NotImplementedError("Must be overridden")


class StupidHeaderLineParser(HeaderLineParserBase):
    """Parse into HeaderLine (no particular structure)"""

    def parse_key_value(self, key: str, value: str) -> header.HeaderLine:
        return header.HeaderLine(key, value)


class MappingHeaderLineParser(HeaderLineParserBase):
    """Parse into HeaderLine (no particular structure)"""

    def __init__(self, line_class: Callable[[str, str, dict[str, bool | str | list[str]]], header.HeaderLine]):
        """Initialize the parser"""
        #: the class to use for the VCF header line
        self.line_class = line_class

    def parse_key_value(self, key: str, value: str) -> header.HeaderLine:
        return self.line_class(key, value, parse_mapping(value))


def build_header_parsers() -> dict[str, HeaderLineParserBase]:
    """Return mapping for parsers to use for each VCF header type

    Inject the WarningHelper into the parsers.
    """
    result: dict[str, HeaderLineParserBase] = {
        "ALT": MappingHeaderLineParser(header.AltAlleleHeaderLine),
        "contig": MappingHeaderLineParser(header.ContigHeaderLine),
        "FILTER": MappingHeaderLineParser(header.FilterHeaderLine),
        "FORMAT": MappingHeaderLineParser(header.FormatHeaderLine),
        "INFO": MappingHeaderLineParser(header.InfoHeaderLine),
        "META": MappingHeaderLineParser(header.MetaHeaderLine),
        "PEDIGREE": MappingHeaderLineParser(header.PedigreeHeaderLine),
        "SAMPLE": MappingHeaderLineParser(header.SampleHeaderLine),
        "__default__": StupidHeaderLineParser(),  # fallback
    }
    return result


# Field value converters
_CONVERTERS: dict[
    Literal["Integer", "Float", "Flag", "Character", "String"], Callable[[str], bool | int | float | str]
] = {
    "Integer": int,
    "Float": float,
    "Flag": lambda x: True,
    "Character": str,
    "String": str,
}


def convert_field_value(
    type_: Literal["Integer", "Float", "Flag", "Character", "String"], value: str
) -> bool | int | float | str | None:
    """Convert atomic field value according to the type"""
    if value == ".":
        return None
    elif type_ in ("Character", "String"):
        if "%" in value:
            for k, v in record.UNESCAPE_MAPPING:
                value = value.replace(k, v)
        return value
    else:
        try:
            return _CONVERTERS[type_](value)
        except ValueError:  # pragma: no cover
            warnings.warn(
                ("{} cannot be converted to {}, keeping as string.").format(value, type_),
                exceptions.CannotConvertValue,
            )
            return value


def parse_field_value(
    field_info: header.FieldInfo, value: str | bool
) -> bool | int | float | str | list[bool | int | float | str | None] | None:
    """Parse ``value`` according to ``field_info``"""
    if isinstance(value, bool) or field_info.type == "Flag":
        return True
    elif field_info.id in ("FORMAT/FT", "FT"):
        return [x for x in value.split(";") if x != "."]
    elif field_info.number == 1:
        return convert_field_value(field_info.type, value)
    else:
        if value == ".":
            return []
        else:
            return [convert_field_value(field_info.type, x) for x in value.split(",")]


# Regular expression for break-end
BREAKEND_PATTERN = re.compile(r"[\[\]]")


def parse_breakend(alt_str: str) -> tuple[str, int, str, Literal["+", "-"], str, bool]:
    """Parse breakend and return tuple with results, parameters for BreakEnd
    constructor
    """
    arr = BREAKEND_PATTERN.split(alt_str)
    assert isinstance(arr[1], str)
    mate_chrom, mate_pos = arr[1].split(":", 1)
    mate_pos = int(mate_pos)
    if mate_chrom[0] == "<":
        mate_chrom = mate_chrom[1:-1]
        within_main_assembly = False
    else:
        within_main_assembly = True
    FWD_REV: dict[bool, Literal["+", "-"]] = {True: record.FORWARD, False: record.REVERSE}
    orientation = FWD_REV[alt_str[0] == "[" or alt_str[0] == "]"]
    mate_orientation = FWD_REV["[" in alt_str]
    assert isinstance(arr[2], str) and isinstance(arr[0], str)
    if orientation == record.FORWARD:
        sequence = arr[2]
    else:
        sequence = arr[0]
    return (mate_chrom, mate_pos, orientation, mate_orientation, sequence, within_main_assembly)


def process_sub_grow(ref: str, alt_str: str) -> record.Substitution:
    """Process substution where the string grows"""
    if len(alt_str) == 0:
        raise exceptions.InvalidRecordException("Invalid VCF, empty ALT")
    elif len(alt_str) == 1:
        if ref[0] == alt_str[0]:
            return record.Substitution(record.DEL, alt_str)
        else:
            return record.Substitution(record.INDEL, alt_str)
    else:
        return record.Substitution(record.INDEL, alt_str)


def process_sub_shrink(ref: str, alt_str: str) -> record.Substitution:
    """Process substution where the string shrink"""
    if len(ref) == 0:  # pragma: no cover
        raise exceptions.InvalidRecordException("Invalid VCF, empty REF")
    elif len(ref) == 1:
        if ref[0] == alt_str[0]:
            return record.Substitution(record.INS, alt_str)
        else:
            return record.Substitution(record.INDEL, alt_str)
    else:
        return record.Substitution(record.INDEL, alt_str)


def process_sub(ref: str, alt_str: str) -> record.Substitution:
    """Process substitution"""
    if len(ref) == len(alt_str):
        if len(ref) == 1:
            return record.Substitution(record.SNV, alt_str)
        else:
            return record.Substitution(record.MNV, alt_str)
    elif len(ref) > len(alt_str):
        return process_sub_grow(ref, alt_str)
    else:  # len(ref) < len(alt_str):
        return process_sub_shrink(ref, alt_str)


def process_alt(header: header.Header, ref: str, alt_str: str) -> record.AltRecord:
    """Process alternative value using Header in ``header``"""
    # By its nature, this function contains a large number of case distinctions
    if "]" in alt_str or "[" in alt_str:
        return record.BreakEnd(*parse_breakend(alt_str))
    elif alt_str[0] == "." and len(alt_str) > 0:
        return record.SingleBreakEnd(record.FORWARD, alt_str[1:])
    elif alt_str[-1] == "." and len(alt_str) > 0:
        return record.SingleBreakEnd(record.REVERSE, alt_str[:-1])
    elif alt_str[0] == "<" and alt_str[-1] == ">":
        inner = alt_str[1:-1]
        return record.SymbolicAllele(inner)
    else:  # substitution
        return process_sub(ref, alt_str)


class HeaderParser:
    """Helper class for parsing a VCF header"""

    def __init__(self):
        #: Sub parsers to use for parsing the header lines
        self.sub_parsers = build_header_parsers()

    def parse_line(self, line: str) -> header.HeaderLine:
        """Parse VCF header ``line`` (trailing '\r\n' or '\n' is ignored)

        :param str line: ``str`` with line to parse
        :param dict sub_parsers: ``dict`` mapping header line types to
            appropriate parser objects
        :returns: appropriate :py:class:`HeaderLine` parsed from ``line``
        :raises: :py:class:`vcfpy.exceptions.InvalidHeaderException` if
            there was a problem parsing the file
        """
        if not line or not line.startswith("##"):  # pragma: no cover
            raise exceptions.InvalidHeaderException('Invalid VCF header line (must start with "##") {}'.format(line))
        if "=" not in line:  # pragma: no cover
            raise exceptions.InvalidHeaderException('Invalid VCF header line (must contain "=") {}'.format(line))
        line = line[len("##") :].rstrip()  # trim '^##' and trailing whitespace
        # split key/value pair at "="
        key, value = split_mapping(line)
        sub_parser = self.sub_parsers.get(key, self.sub_parsers["__default__"])
        return sub_parser.parse_key_value(key, value)


class RecordParser:
    """Helper class for parsing VCF records"""

    def __init__(
        self,
        header: header.Header,
        samples: header.SamplesInfos,
        record_checks: Iterable[Literal["FORMAT", "INFO"]] | None = None,
    ):
        #: Header with the meta information
        self.header = header
        #: SamplesInfos with sample information
        self.samples = samples
        #: The checks to perform, can contain 'INFO' and 'FORMAT'
        self.record_checks = tuple(record_checks or [])
        # Expected number of fields
        if self.samples.names:
            self.expected_fields = 9 + len(self.samples.names)
        else:
            self.expected_fields = 8
        # Cache of FieldInfo objects by FORMAT string
        self._format_cache: dict[str, list["header.FieldInfo"]] = {}
        # Cache of FILTER entries, also applied to FORMAT/FT
        self._filter_ids = set(self.header.filter_ids())
        # Helper for checking INFO fields
        if "INFO" in self.record_checks:
            self._info_checker = InfoChecker(self.header)
        else:
            self._info_checker = NoopInfoChecker()
        # Helper for checking FORMAT fields
        if "FORMAT" in self.record_checks:
            self._format_checker = FormatChecker(self.header)
        else:
            self._format_checker = NoopFormatChecker()

    def parse_line(self, line_str: str) -> "record.Record | None":
        """Parse line from file (including trailing line break) and return
        resulting Record
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
        if arr[2] == ".":
            ids = []
        else:
            ids = arr[2].split(";")
        # REF
        ref = arr[3]
        # ALT
        alts: list[record.AltRecord] = []
        if arr[4] != ".":
            for alt in arr[4].split(","):
                alts.append(process_alt(self.header, ref, alt))
        # QUAL
        if arr[5] == ".":
            qual = None
        else:
            try:
                qual = int(arr[5])
            except ValueError:  # try as float  # pragma: no cover
                qual = float(arr[5])
        # FILTER
        if arr[6] == ".":
            filt = []
        else:
            filt = arr[6].split(";")
        self._check_filters(filt, "FILTER")
        # INFO
        info = self._parse_info(arr[7], len(alts))
        if len(arr) == 9:
            raise exceptions.IncorrectVCFFormat("Expected 8 or 10+ columns, got 9!")  # pragma: no cover
        elif len(arr) == 8:
            format_ = None
            calls = None
        else:
            # FORMAT
            format_ = arr[8].split(":")
            # sample/call columns
            calls = self._handle_calls(alts, format_, arr[8], arr)
        return record.Record(chrom, pos, ids, ref, alts, qual, filt, info, format_, calls)

    def _handle_calls(
        self, alts: list[record.AltRecord], format_: list[str], format_str: str, arr: list[str]
    ) -> list["record.Call | record.UnparsedCall"]:
        """Handle FORMAT and calls columns, factored out of parse_line"""
        if format_str not in self._format_cache:
            self._format_cache[format_str] = list(map(self.header.get_format_field_info, format_))
        # per-sample calls
        calls: list["record.Call | record.UnparsedCall"] = []
        for sample, raw_data in zip(self.samples.names, arr[9:], strict=False):
            if self.samples.is_parsed(sample):
                data = self._parse_calls_data(format_, self._format_cache[format_str], raw_data)
                call = record.Call(sample, data)
                self._format_checker.run(call, len(alts))
                ft_value = call.data.get("FT") or []
                if not isinstance(ft_value, list):  # pragma: no cover
                    raise ValueError("FORMAT/FT field must be a list of strings but was {}".format(repr(ft_value)))
                ft_value_ = cast(list[str], ft_value)
                if not all(isinstance(x, str) for x in ft_value_):  # pragma: no cover
                    raise ValueError("FORMAT/FT field must be a list of strings but was {}".format(repr(ft_value_)))
                self._check_filters(ft_value_, "FORMAT/FT", call.sample)
                calls.append(call)
            else:
                calls.append(record.UnparsedCall(sample, raw_data))
        return calls

    def _check_filters(self, filt: list[str], source: str, sample: str | None = None):
        if not filt:
            return
        for f in filt:
            self._check_filter(f, source, sample)

    def _check_filter(self, f: str, source: str, sample: str | None):
        if f == "PASS":
            pass  # the PASS filter is implicitely defined
        elif f not in self._filter_ids:  # pragma: no cover
            if source == "FILTER":
                warnings.warn(
                    ("Filter not found in header: {}; problem in FILTER column").format(f),
                    exceptions.UnknownFilter,
                )
            else:
                assert source == "FORMAT/FT" and sample
                warnings.warn(
                    ("Filter not found in header: {}; problem in FORMAT/FT column of sample {}").format(f, sample),
                    exceptions.UnknownFilter,
                )

    def _split_line(self, line_str: str) -> list[str]:
        """Split line and check number of columns"""
        arr = line_str.rstrip().split("\t")
        if len(arr) != self.expected_fields:
            raise exceptions.InvalidRecordException(
                "The line contains an invalid number of fields. Was {} but expected {}\n{}".format(
                    len(arr), self.expected_fields, line_str
                )
            )
        return arr

    def _parse_info(self, info_str: str, num_alts: int) -> dict[str, Any]:
        """Parse INFO column from string"""
        result: dict[str, Any] = {}
        if info_str == ".":
            return result
        # The standard is very nice to parsers, we can simply split at
        # semicolon characters, although I (Manuel) don't know how strict
        # programs follow this
        for entry in info_str.split(";"):
            if "=" not in entry:  # flag
                key = entry
                result[key] = parse_field_value(self.header.get_info_field_info(key), True)
            else:
                key, value = split_mapping(entry)
                result[key] = parse_field_value(self.header.get_info_field_info(key), value)
            self._info_checker.run(key, result[key], num_alts)
        return result

    @classmethod
    def _parse_calls_data(
        cls, format_: list[str], infos: list["header.FieldInfo"], gt_str: str
    ) -> dict[str, bool | int | float | str | list[bool | int | float | str | None] | None]:
        """Parse genotype call information from arrays using format array

        :param list format: List of strings with format names
        :param gt_str arr: string with genotype information values
        """
        data: dict[str, bool | int | float | str | list[bool | int | float | str | None] | None] = {}
        # The standard is very nice to parsers, we can simply split at
        # colon characters, although I (Manuel) don't know how strict
        # programs follow this
        for key, info, value in zip(format_, infos, gt_str.split(":"), strict=False):
            data[key] = parse_field_value(info, value)
        return data


class HeaderChecker:
    """Helper class for checking a VCF header"""

    def run(self, header: header.Header) -> None:
        """Check the header

        Warnings will be printed using ``warnings`` while errors will raise
        an exception.

        :raises: ``vcfpy.exceptions.InvalidHeaderException`` in the case of
            severe errors reading the header
        """
        self._check_header_lines(header.lines)

    def _check_header_lines(self, header_lines: list[header.HeaderLine]) -> None:
        """Check header lines, in particular for starting file "##fileformat" """
        if not header_lines:
            raise exceptions.InvalidHeaderException(
                "The VCF file did not contain any header lines!"
            )  # pragma: no cover
        first = header_lines[0]
        if first.key != "fileformat":
            raise exceptions.InvalidHeaderException("The VCF file did not start with ##fileformat")
        if first.value not in SUPPORTED_VCF_VERSIONS:
            warnings.warn("Unknown VCF version {}".format(first.value), exceptions.UnknownVCFVersion)


@functools.lru_cache(maxsize=32)
def binomial(n: int, k: int):
    try:
        res = math.factorial(n) // math.factorial(k) // math.factorial(n - k)
    except ValueError:
        res = 0
    return res


class AbstractInfoChecker:
    """Abstract base class for INFO field checkers"""

    def run(self, key: str, value: str, num_alts: int) -> None:
        """Run the checker"""
        raise NotImplementedError  # pragma: no cover


class NoopInfoChecker(AbstractInfoChecker):
    """Helper class that performs no checks"""

    def run(self, key: str, value: str, num_alts: int) -> None:
        pass


class InfoChecker(AbstractInfoChecker):
    """Helper class for checking an INFO field"""

    def __init__(self, header: header.Header):
        #: VCFHeader to use for checking
        self.header = header

    def run(self, key: str, value: str, num_alts: int) -> None:
        """Check value in INFO[key] of record

        Currently, only checks for consistent counts are implemented

        :param str key: key of INFO entry to check
        :param value: value to check
        :param int alts: list of alternative alleles, for length
        """
        field_info = self.header.get_info_field_info(key)
        if not isinstance(value, list):
            return
        TABLE = {
            ".": len(value),
            "A": num_alts,
            "R": num_alts + 1,
            "G": binomial(num_alts + 1, 2),  # diploid only at the moment
        }
        expected = TABLE.get(str(field_info.number), field_info.number)
        if len(value) != expected:
            tpl = "Number of elements for INFO field {} is {} instead of {}"
            warnings.warn(tpl.format(key, len(value), field_info.number), exceptions.IncorrectListLength)


class AbstractNoopFormatChecker:
    """Abstract base class for FORMAT field checkers"""

    def run(self, call: "record.Call", num_alts: int) -> None:
        raise NotImplementedError  # pragma: no cover


class NoopFormatChecker(AbstractNoopFormatChecker):
    """Helper class that performs no checks"""

    def run(self, call: "record.Call", num_alts: int) -> None:
        pass


class FormatChecker(AbstractNoopFormatChecker):
    """Helper class for checking a FORMAT field"""

    def __init__(self, header: header.Header):
        #: VCFHeader to use for checking
        self.header = header

    def run(self, call: "record.Call", num_alts: int) -> None:
        """Check ``FORMAT`` of a record.Call

        Currently, only checks for consistent counts are implemented
        """
        for key, value in call.data.items():
            self._check_count(call, key, value, num_alts)

    def _check_count(self, call: "record.Call", key: str, value: str, num_alts: int) -> None:
        field_info = self.header.get_format_field_info(key)
        if field_info.id == "GT":
            return
        if isinstance(value, list):
            return
        num_alleles = len(call.gt_alleles or [])
        TABLE = {
            ".": len(value),
            "A": num_alts,
            "R": num_alts + 1,
            "G": binomial(num_alts + num_alleles, num_alleles),
        }
        expected = TABLE.get(str(field_info.number), field_info.number)
        if len(value) != expected:
            tpl = "Number of elements for FORMAT field {} is {} instead of {} (number specifier {})"
            warnings.warn(
                tpl.format(key, len(value), expected, field_info.number),
                exceptions.IncorrectListLength,
            )


class Parser:
    """Class for line-wise parsing of VCF files

    In most cases, you want to use :py:class:`vcfpy.reader.Reader` instead.

    :param stream: ``file``-like object to read from
    :param str path: path the VCF is parsed from, for display purposes
        only, optional
    """

    def __init__(
        self,
        stream: "io.TextIOWrapper",
        path: pathlib.Path | str | None = None,
        record_checks: Iterable[Literal["FORMAT", "INFO"]] | None = None,
    ):
        self.stream = stream
        self.path = None if path is None else str(path)
        #: checks to perform, can contain 'INFO' and 'FORMAT'
        self.record_checks = tuple(record_checks or []) or None
        #: header, once it has been read
        self.header = None
        # the currently read line
        self._line = stream.readline()  # trailing '\n'
        #: :py:class:`vcfpy.header.SamplesInfos` with sample information;
        #: set on parsing the header
        self.samples = None
        # helper for parsing the records
        self._record_parser = None
        # helper for checking the header
        self._header_checker = HeaderChecker()

    def _read_next_line(self):
        """Read next line store in self._line and return old one"""
        prev_line = self._line
        self._line = self.stream.readline()
        return prev_line

    def parse_header(self, parsed_samples: list[str] | None = None):
        """Read and parse :py:class:`vcfpy.header.Header` from file, set
        into ``self.header`` and return it

        :param list parsed_samples: ``list`` of ``str`` for subsetting the
            samples to parse
        :returns: ``vcfpy.header.Header``
        :raises: ``vcfpy.exceptions.InvalidHeaderException`` in the case of
            problems reading the header
        """
        # parse header lines
        sub_parser = HeaderParser()
        header_lines: list[header.HeaderLine] = []
        while self._line and self._line.startswith("##"):
            header_lines.append(sub_parser.parse_line(self._line))
            self._read_next_line()
        # parse sample info line
        self.samples = self._handle_sample_line(parsed_samples)
        # construct Header object
        self.header = header.Header(header_lines, self.samples)
        # check header for consistency
        self._header_checker.run(self.header)
        # construct record parser
        self._record_parser = RecordParser(self.header, self.samples, self.record_checks)
        # read next line, must not be header
        self._read_next_line()
        if self._line and self._line.startswith("#"):  # pragma: no cover
            raise exceptions.IncorrectVCFFormat('Expecting non-header line or EOF after "#CHROM" line')
        return self.header

    def _handle_sample_line(self, parsed_samples: list[str] | None = None):
        """ "Check and interpret the "##CHROM" line and return samples"""
        if not self._line or not self._line.startswith("#CHROM"):  # pragma: no cover
            raise exceptions.IncorrectVCFFormat('Missing line starting with "#CHROM"')
        # check for space before INFO
        line = self._line.rstrip()
        pos = line.find("FORMAT") if ("FORMAT" in line) else line.find("INFO")
        if pos == -1:  # pragma: no cover
            raise exceptions.IncorrectVCFFormat('Ill-formatted line starting with "#CHROM"')
        if " " in line[:pos]:
            warnings.warn(
                "Found space in #CHROM line, splitting at whitespace instead of tab; this VCF file is ill-formatted",
                exceptions.SpaceInChromLine,
            )
            arr = self._line.rstrip().split()
        else:
            arr = self._line.rstrip().split("\t")

        self._check_samples_line(arr)
        return header.SamplesInfos(arr[len(REQUIRE_SAMPLE_HEADER) :], parsed_samples)

    @classmethod
    def _check_samples_line(cls, arr: list[str]):
        """Peform additional check on samples line"""
        if len(arr) <= len(REQUIRE_NO_SAMPLE_HEADER):
            if tuple(arr) != REQUIRE_NO_SAMPLE_HEADER:
                raise exceptions.IncorrectVCFFormat(  # pragma: no cover
                    "Sample header line indicates no sample but does not equal required prefix {}".format(
                        "\t".join(REQUIRE_NO_SAMPLE_HEADER)
                    )
                )
        elif tuple(arr[: len(REQUIRE_SAMPLE_HEADER)]) != REQUIRE_SAMPLE_HEADER:
            raise exceptions.IncorrectVCFFormat(  # pragma: no cover
                'Sample header line (starting with "#CHROM") does not start with required prefix {}'.format(
                    "\t".join(REQUIRE_SAMPLE_HEADER)
                )
            )

    def parse_line(self, line: str):
        """Parse the given line without reading another one from the stream"""
        if self._record_parser is None:
            raise exceptions.InvalidRecordException("Cannot parse record before parsing header")  # pragma: no cover
        return self._record_parser.parse_line(line)

    def parse_next_record(self):
        """Read, parse and return next :py:class:`vcfpy.record.Record`

        :returns: next VCF record or ``None`` if at end
        :raises: ``vcfpy.exceptions.InvalidRecordException`` in the case of
            problems reading the record
        """
        return self.parse_line(self._read_next_line())

    def print_warn_summary(self):
        """If there were any warnings, print summary with warnings"""
        # TODO: remove?
