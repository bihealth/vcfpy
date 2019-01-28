# -*- coding: utf-8 -*-
"""Code for representing the VCF header part

The VCF header class structure is modeled after HTSJDK
"""

import json
import pprint
import warnings

from . import exceptions
from .compat import OrderedDict
from .exceptions import (
    DuplicateHeaderLineWarning,
    FieldInfoNotFound,
    FieldMissingNumber,
    FieldInvalidNumber,
    HeaderInvalidType,
    HeaderMissingDescription,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# Tuples of valid entries -----------------------------------------------------
#
#: valid INFO value types
INFO_TYPES = ("Integer", "Float", "Flag", "Character", "String")
#: valid FORMAT value types
FORMAT_TYPES = ("Integer", "Float", "Character", "String")
#: valid values for "Number" entries, except for integers
VALID_NUMBERS = ("A", "R", "G", ".")
#: header lines that contain an "ID" entry
LINES_WITH_ID = ("ALT", "contig", "FILTER", "FORMAT", "INFO", "META", "PEDIGREE", "SAMPLE")

# Constants for "Number" entries ----------------------------------------------
#
#: number of alleles excluding reference
HEADER_NUMBER_ALLELES = "A"
#: number of alleles including reference
HEADER_NUMBER_REF = "R"
#: number of genotypes
HEADER_NUMBER_GENOTYPES = "G"
#: unbounded number of values
HEADER_NUMBER_UNBOUNDED = "."


class FieldInfo:
    """Core information for describing field type and number"""

    # TODO: always put in id?
    def __init__(self, type_, number, description=None, id_=None):
        #: The type, one of INFO_TYPES or FORMAT_TYPES
        self.type = type_
        #: Number description, either an int or constant
        self.number = number
        #: Description for the header field, optional
        self.description = description
        #: The id of the field, optional.
        self.id = id_

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        return NotImplemented

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        return "FieldInfo({}, {}, {}, {})".format(
            *map(repr, [self.type, self.number, self.description, self.id])
        )

    def __repr__(self):
        return str(self)


# Reserved INFO keys ----------------------------------------------------------

#: Reserved fields for INFO from VCF v4.3
RESERVED_INFO = {
    # VCF v4.3, Section 1.6.1
    "AA": FieldInfo("String", 1, "Ancestral Allele"),
    "AC": FieldInfo(
        "Integer",
        "A",
        "Allele count in genotypes, for each ALT allele, in the " "same order as listed",
    ),
    "AD": FieldInfo("Integer", "R", "Total read depth for each allele"),
    "ADF": FieldInfo("Integer", "R", "Forward read depth for each allele"),
    "ADR": FieldInfo("Integer", "R", "Reverse read depth for each allele"),
    "AF": FieldInfo(
        "Float",
        "A",
        "Allele frequency for each ALT allele in the same order "
        "as listed: used for estimating from primary data not "
        "called genotypes",
    ),
    "AN": FieldInfo("Integer", 1, "Total number of alleles in called genotypes"),
    "BQ": FieldInfo("Float", 1, "RMS base quality at this position"),
    "CIGAR": FieldInfo(
        "String",
        "A",
        "CIGAR string describing how to align each ALT allele " "to the reference allele",
    ),
    "DB": FieldInfo("Flag", 0, "dbSNP membership"),
    "DP": FieldInfo(
        "Integer",
        1,
        "Combined depth across samples for small variants and "
        "Read Depth of segment containing breakend for SVs",
    ),
    "H2": FieldInfo("Flag", 0, "Membership in HapMap 2"),
    "H3": FieldInfo("Flag", 0, "Membership in HapMap 3"),
    "MQ": FieldInfo("Integer", 1, "RMS mapping quality"),
    "MQ0": FieldInfo("Integer", 1, "Number of MAPQ == 0 reads covering this record"),
    "NS": FieldInfo("Integer", 1, "Number of samples with data"),
    "SB": FieldInfo("Integer", 4, "Strand bias at this position"),
    "SOMATIC": FieldInfo(
        "Flag", 0, "Indicates that the record is a somatic mutation, " "for cancer genomics"
    ),
    "VALIDATED": FieldInfo("Flag", 0, "Validated by follow-up experiment"),
    "1000G": FieldInfo("Flag", 0, "Membership in 1000 Genomes"),
    # VCF v4.3, Section 3
    "IMPRECISE": FieldInfo("Flag", 0, "Imprecise structural variation"),
    "NOVEL": FieldInfo("Flag", 0, "Indicates a novel structural variation"),
    "END": FieldInfo(
        "Integer",
        1,
        "End position of the variant described in this record " "(for symbolic alleles)",
    ),
    "SVTYPE": FieldInfo("String", 1, "Type of structural variant"),
    "SVLEN": FieldInfo("Integer", 1, "Difference in length between REF and ALT alleles"),
    "CIPOS": FieldInfo("Integer", 2, "Confidence interval around POS for imprecise " "variants"),
    "CIEND": FieldInfo("Integer", 2, "Confidence interval around END for imprecise " "variants"),
    "HOMLEN": FieldInfo(
        "Integer", ".", "Length of base pair identical micro-homology at " "event breakpoints"
    ),
    "HOMSEQ": FieldInfo(
        "String", ".", "Sequence of base pair identical micro-homology at " "event breakpoints"
    ),
    "BKPTID": FieldInfo(
        "String", ".", "ID of the assembled alternate allele in the " "assembly file"
    ),
    "MEINFO": FieldInfo("String", 4, "Mobile element info of the form " "NAME,START,END,POLARITY"),
    "METRANS": FieldInfo(
        "String", 4, "Mobile element transduction info of the form " "CHR,START,END,POLARITY"
    ),
    "DGVID": FieldInfo("String", 1, "ID of this element in Database of Genomic Variation"),
    "DBVARID": FieldInfo("String", 1, "ID of this element in DBVAR"),
    "DBRIPID": FieldInfo("String", 1, "ID of this element in DBRIP"),
    "MATEID": FieldInfo("String", ".", "ID of mate breakends"),
    "PARID": FieldInfo("String", 1, "ID of partner breakend"),
    "EVENT": FieldInfo("String", 1, "ID of event associated to breakend"),
    "CILEN": FieldInfo(
        "Integer", 2, "Confidence interval around the inserted material " "between breakends"
    ),
    "DPADJ": FieldInfo("Integer", ".", "Read Depth of adjacency"),
    "CN": FieldInfo("Integer", 1, "Copy number of segment containing breakend"),
    "CNADJ": FieldInfo("Integer", ".", "Copy number of adjacency"),
    "CICN": FieldInfo("Integer", 2, "Confidence interval around copy number for the " "segment"),
    "CICNADJ": FieldInfo(
        "Integer", ".", "Confidence interval around copy number for the " "adjacency"
    ),
}

# Reserved FORMAT keys --------------------------------------------------------

RESERVED_FORMAT = {
    # VCF v 4.3, Section 1.6.2
    "AD": FieldInfo("Integer", "R", "Total, per-sample read depth"),
    "ADF": FieldInfo("Integer", "R", "Forward-strand, per-sample read depth"),
    "ADR": FieldInfo("Integer", "R", "Reverse-strand, per-sample read depth"),
    "DP": FieldInfo("Integer", 1, "Read depth at this position for this sample"),
    "EC": FieldInfo(
        "Integer", "A", "Expected alternate allele counts for each alternate " "allele"
    ),
    "FT": FieldInfo("String", "1", "Filters applied for this sample", "FORMAT/FT"),
    "GQ": FieldInfo("Integer", "G", "Phred-scale, conditional genotype quality"),
    "GP": FieldInfo("Float", "G", "Genotype posterior probabilities"),
    "GT": FieldInfo("String", 1, "Genotype call"),
    "GL": FieldInfo("Float", "G", "Log10-scaled likelihoods for genotypes"),
    "HQ": FieldInfo("Integer", 2, "Haplotype qualities"),
    "MQ": FieldInfo("Integer", 1, "RMS mapping quality"),
    "PL": FieldInfo("Integer", "G", "Phred-scaled genotype likelihoods, rounded to integers"),
    "PQ": FieldInfo("Integer", 1, "Phasing quality"),
    "PS": FieldInfo(
        "Integer",
        1,
        "Non-negative 32 bit integer giving phasing set " "for this sample and this chromosome",
    ),
    # VCF v4.3, Section 4
    "CN": FieldInfo("Integer", 1, "Copy number genotype for imprecise events"),
    "CNQ": FieldInfo("Float", 1, "Copy number genotype quality for imprecise events"),
    "CNL": FieldInfo("Float", "G", "Copy number genotype likelihood for imprecise events"),
    "CNP": FieldInfo("Float", "G", "Copy number posterior probabilities"),
    "NQ": FieldInfo("Integer", 1, "Phred style probability score that the variant is novel"),
    "HAP": FieldInfo("Integer", 1, "Unique haplotype identifier"),
    "AHAP": FieldInfo("Integer", 1, "Unique identifier of ancestral haplotype"),
}


# header files to enforce double-quoting for
QUOTE_FIELDS = ("Description", "Source", "Version")


def serialize_for_header(key, value):
    """Serialize value for the given mapping key for a VCF header line"""
    if key in QUOTE_FIELDS:
        return json.dumps(value)
    elif isinstance(value, str):
        if " " in value or "\t" in value:
            return json.dumps(value)
        else:
            return value
    elif isinstance(value, list):
        return "[{}]".format(", ".join(value))
    else:
        return str(value)


def header_without_lines(header, remove):
    """Return :py:class:`Header` without lines given in ``remove``

    ``remove`` is an iterable of pairs ``key``/``ID`` with the VCF header key
    and ``ID`` of entry to remove.  In the case that a line does not have
    a ``mapping`` entry, you can give the full value to remove.

    .. code-block:: python

        # header is a vcfpy.Header, e.g., as read earlier from file
        new_header = vcfpy.without_header_lines(
            header, [('assembly', None), ('FILTER', 'PASS')])
        # now, the header lines starting with "##assembly=" and the "PASS"
        # filter line will be missing from new_header
    """
    remove = set(remove)
    # Copy over lines that are not removed
    lines = []
    for line in header.lines:
        if hasattr(line, "mapping"):
            if (line.key, line.mapping.get("ID", None)) in remove:
                continue  # filter out
        else:
            if (line.key, line.value) in remove:
                continue  # filter out
        lines.append(line)
    return Header(lines, header.samples)


class Header:
    """Represent header of VCF file

    While this class allows mutating records, it should not be changed once it
    has been assigned to a writer.  Use :py:method:`~Header.copy` to create
    a copy that can be modified without problems.

    This class provides function for adding lines to a header and updating the
    supporting index data structures.  There is no explicit API for removing
    header lines, the best way is to reconstruct a new ``Header`` instance with
    a filtered list of header lines.
    """

    def __init__(self, lines=None, samples=None):
        #: ``list`` of :py:HeaderLine objects
        self.lines = lines or []
        #: :py:class:`SamplesInfo` object
        self.samples = samples
        # build indices for the different field types
        self._indices = self._build_indices()

    def _build_indices(self):
        """Build indices for the different field types"""
        result = {key: OrderedDict() for key in LINES_WITH_ID}
        for line in self.lines:
            if line.key in LINES_WITH_ID:
                result.setdefault(line.key, OrderedDict())
                if line.mapping["ID"] in result[line.key]:
                    warnings.warn(
                        ("Seen {} header more than once: {}, using first" "occurence").format(
                            line.key, line.mapping["ID"]
                        ),
                        DuplicateHeaderLineWarning,
                    )
                else:
                    result[line.key][line.mapping["ID"]] = line
            else:
                result.setdefault(line.key, [])
                result[line.key].append(line)
        return result

    def copy(self):
        """Return a copy of this header"""
        return Header([line.copy() for line in self.lines], self.samples.copy())

    def add_filter_line(self, mapping):
        """Add FILTER header line constructed from the given mapping

        :param mapping: ``OrderedDict`` with mapping to add.  It is
            recommended to use ``OrderedDict`` over ``dict`` as this makes
            the result reproducible
        :return: ``False`` on conflicting line and ``True`` otherwise
        """
        return self.add_line(FilterHeaderLine.from_mapping(mapping))

    def add_contig_line(self, mapping):
        """Add "contig" header line constructed from the given mapping

        :param mapping: ``OrderedDict`` with mapping to add.  It is
            recommended to use ``OrderedDict`` over ``dict`` as this makes
            the result reproducible
        :return: ``False`` on conflicting line and ``True`` otherwise
        """
        return self.add_line(ContigHeaderLine.from_mapping(mapping))

    def add_info_line(self, mapping):
        """Add INFO header line constructed from the given mapping

        :param mapping: ``OrderedDict`` with mapping to add.  It is
            recommended to use ``OrderedDict`` over ``dict`` as this makes
            the result reproducible
        :return: ``False`` on conflicting line and ``True`` otherwise
        """
        return self.add_line(InfoHeaderLine.from_mapping(mapping))

    def add_format_line(self, mapping):
        """Add FORMAT header line constructed from the given mapping

        :param mapping: ``OrderedDict`` with mapping to add.  It is
            recommended to use ``OrderedDict`` over ``dict`` as this makes
            the result reproducible
        :return: ``False`` on conflicting line and ``True`` otherwise
        """
        return self.add_line(FormatHeaderLine.from_mapping(mapping))

    def format_ids(self):
        """Return list of all format IDs"""
        return list(self._indices["FORMAT"].keys())

    def filter_ids(self):
        """Return list of all filter IDs"""
        return list(self._indices["FILTER"].keys())

    def info_ids(self):
        """Return list of all info IDs"""
        return list(self._indices["INFO"].keys())

    def get_lines(self, key):
        """Return header lines having the given ``key`` as their type"""
        if key in self._indices:
            return self._indices[key].values()
        else:
            return []

    def has_header_line(self, key, id_):
        """Return whether there is a header line with the given ID of the
        type given by ``key``

        :param key: The VCF header key/line type.
        :param id_: The ID value to compare fore

        :return: ``True`` if there is a header line starting with ``##${key}=``
            in the VCF file having the mapping entry ``ID`` set to ``id_``.
        """
        if key not in self._indices:
            return False
        else:
            return id_ in self._indices[key]

    def add_line(self, header_line):
        """Add header line, updating any necessary support indices

        :return: ``False`` on conflicting line and ``True`` otherwise
        """
        self.lines.append(header_line)
        self._indices.setdefault(header_line.key, OrderedDict())
        if not hasattr(header_line, "mapping"):
            return False  # no registration required
        if self.has_header_line(header_line.key, header_line.mapping["ID"]):
            warnings.warn(
                (
                    "Detected duplicate header line with type {} and ID {}. "
                    "Ignoring this and subsequent one"
                ).format(header_line.key, header_line.mapping["ID"]),
                DuplicateHeaderLineWarning,
            )
            return False
        else:
            self._indices[header_line.key][header_line.mapping["ID"]] = header_line
            return True

    def get_info_field_info(self, key):
        """Return :py:class:`FieldInfo` for the given INFO field"""
        return self._get_field_info("INFO", key)

    def get_format_field_info(self, key):
        """Return :py:class:`FieldInfo` for the given INFO field"""
        return self._get_field_info("FORMAT", key)

    def _get_field_info(self, type_, key):
        result = self._indices[type_].get(key)
        if result:
            return result
        if key in RESERVED_INFO:
            res = FieldInfo(RESERVED_INFO[key].type, RESERVED_INFO[key].number)
        else:
            res = FieldInfo("String", HEADER_NUMBER_UNBOUNDED)
        warnings.warn(
            "{} {} not found using {}/{} instead".format(type_, key, res.type, repr(res.number)),
            FieldInfoNotFound,
        )
        return res

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.lines, self.samples) == (other.lines, other.samples)
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return (self.lines, self.samples) != (other.lines, other.samples)
        return NotImplemented

    def __hash__(self):
        raise TypeError("Unhashable type: Header")

    def __str__(self):
        tpl = "Header(lines={}, samples={})"
        return tpl.format(*map(repr, (self.lines, self.samples)))

    def __repr__(self):
        return str(self)


class HeaderLine:
    """Base class for VCF header lines
    """

    def __init__(self, key, value):
        #: ``str`` with key of header line
        self.key = key
        # ``str`` with raw value of header line
        self._value = value

    def copy(self):
        """Return a copy"""
        return self.__class__(self.key, self.value)

    @property
    def value(self):
        return self._value

    def serialize(self):
        """Return VCF-serialized version of this header line"""
        return "".join(("##", self.key, "=", self.value))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.value) == (other.key, other.value)
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.value) != (other.key, other.value)
        return NotImplemented

    def __hash__(self):
        raise TypeError("Unhashable type: HeaderLine")

    def __str__(self):
        return "HeaderLine({}, {})".format(*map(repr, (self.key, self.value)))

    def __repr__(self):
        return str(self)


def mapping_to_str(mapping):
    """Convert mapping to string"""
    result = ["<"]
    for i, (key, value) in enumerate(mapping.items()):
        if i > 0:
            result.append(",")
        result += [key, "=", serialize_for_header(key, value)]
    result += [">"]
    return "".join(result)


class SimpleHeaderLine(HeaderLine):
    """Base class for simple header lines, currently contig and filter
    header lines

    Don't use this class directly but rather the sub classes.

    :raises: :py:class:`vcfpy.exceptions.InvalidHeaderException` in
        the case of missing key ``"ID"``
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value)
        # check existence of key "ID"
        if "ID" not in mapping:
            raise exceptions.InvalidHeaderException(
                'Missing key "ID" in header line "{}={}"'.format(key, value)
            )
        #: ``collections.OrderedDict`` with key/value mapping of the attributes
        self.mapping = OrderedDict(mapping.items())

    def copy(self):
        """Return a copy"""
        mapping = OrderedDict(self.mapping.items())
        return self.__class__(self.key, self.value, mapping)

    @property
    def value(self):
        return mapping_to_str(self.mapping)

    def serialize(self):
        return "".join(map(str, ["##", self.key, "=", self.value]))

    def __str__(self):
        return "SimpleHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.value, self.mapping) == (other.key, other.value, other.mapping)
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.value, self.mapping) != (other.key, other.value, other.mapping)
        return NotImplemented


class AltAlleleHeaderLine(SimpleHeaderLine):
    """Alternative allele header line

    Mostly used for defining symbolic alleles for structural variants and
    IUPAC ambiguity codes
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return AltAlleleHeaderLine("ALT", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: name of the alternative allele
        self.id = self.mapping["ID"]

    def __hash__(self):
        raise TypeError("Unhashable type: AltAlleleHeaderLine")

    def __str__(self):
        return "AltAlleleHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class ContigHeaderLine(SimpleHeaderLine):
    """Contig header line

    Most importantly, parses the ``'length'`` key into an integer
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return ContigHeaderLine("contig", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        # convert 'length' entry to integer if possible
        if "length" in self.mapping:
            mapping["length"] = int(mapping["length"])
        else:
            warnings.warn(
                'Field "length" not found in header line {}={}'.format(key, value),
                FieldInfoNotFound,
            )
        #: name of the contig
        self.id = self.mapping["ID"]
        #: length of the contig, ``None`` if missing
        self.length = self.mapping.get("length")

    def __hash__(self):
        raise TypeError("Unhashable type: ContigHeaderLine")

    def __str__(self):
        return "ContigHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class FilterHeaderLine(SimpleHeaderLine):
    """FILTER header line
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return FilterHeaderLine("FILTER", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        # check for "Description" key
        if "Description" not in self.mapping:
            warnings.warn(
                'Field "Description" not found in header line {}={}'.format(key, value),
                FieldInfoNotFound,
            )
        #: token for the filter
        self.id = self.mapping["ID"]
        #: description for the filter, ``None`` if missing
        self.description = self.mapping.get("Description")

    def __hash__(self):
        raise TypeError("Unhashable type: FilterHeaderLine")

    def __str__(self):
        return "FilterHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class MetaHeaderLine(SimpleHeaderLine):
    """Alternative allele header line

    Used for defining set of valid values for samples keys
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return MetaHeaderLine("META", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: name of the alternative allele
        self.id = self.mapping["ID"]

    def __hash__(self):
        raise TypeError("Unhashable type: MetaHeaderLine")

    def __str__(self):
        return "MetaHeaderLine({}, {}, {})".format(*map(repr, (self.key, self.value, self.mapping)))


class PedigreeHeaderLine(SimpleHeaderLine):
    """Header line for defining a pedigree entry
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return PedigreeHeaderLine("PEDIGREE", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: name of the alternative allele
        self.id = self.mapping["ID"]

    def __hash__(self):
        raise TypeError("Unhashable type: PedigreeHeaderLine")

    def __str__(self):
        return "PedigreeHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class SampleHeaderLine(SimpleHeaderLine):
    """Header line for defining a SAMPLE entry
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return SampleHeaderLine("SAMPLE", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: name of the alternative allele
        self.id = self.mapping["ID"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.value, self.mapping) == (other.key, other.value, other.mapping)
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.value, self.mapping) != (other.key, other.value, other.mapping)
        return NotImplemented

    def __hash__(self):
        raise TypeError("Unhashable type: SampleHeaderLine")

    def __str__(self):
        return "SampleHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class CompoundHeaderLine(HeaderLine):
    """Base class for compound header lines, currently format and header lines

    Compound header lines describe fields that can have more than one entry.

    Don't use this class directly but rather the sub classes.
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value)
        #: OrderedDict with key/value mapping
        self.mapping = OrderedDict(mapping.items())
        # check that 'Number' is given and use "." otherwise
        if "Number" not in self.mapping:
            warnings.warn(
                '[vcfpy] WARNING: missing number, using unbounded/"." instead', FieldMissingNumber
            )
            self.mapping["Number"] = "."
        try:
            self.mapping["Number"] = self._parse_number(self.mapping["Number"])
        except ValueError:
            warnings.warn(
                ("[vcfpy] WARNING: invalid number {}, using " 'unbounded/"." instead').format(
                    self.mapping["Number"]
                ),
                FieldInvalidNumber,
            )
            self.mapping["Number"] = "."

    def copy(self):
        """Return a copy"""
        mapping = OrderedDict(self.mapping.items())
        return self.__class__(self.key, self.value, mapping)

    @classmethod
    def _parse_number(klass, number):
        """Parse ``number`` into an ``int`` or return ``number`` if a valid
        expression for a INFO/FORMAT "Number".

        :param str number: ``str`` to parse and check
        """
        try:
            return int(number)
        except ValueError as e:
            if number in VALID_NUMBERS:
                return number
            else:
                raise e

    @property
    def value(self):
        return mapping_to_str(self.mapping)

    def serialize(self):
        return "".join(map(str, ["##", self.key, "=", self.value]))

    def __str__(self):
        return "CompoundHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class InfoHeaderLine(CompoundHeaderLine):
    """Header line for INFO fields

    Note that the ``Number`` field will be parsed into an ``int`` if
    possible.  Otherwise, the constants ``HEADER_NUMBER_*`` will be used.
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return InfoHeaderLine("INFO", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: key in the INFO field
        self.id = self.mapping["ID"]
        # check for "Number" field
        self.number = self.mapping["Number"]
        # check for "Type" field
        type_ = self.mapping.get("Type")
        if "Type" not in self.mapping:
            warnings.warn(
                ('Field "Type" not found in header line, using String ' "instead {}={}").format(
                    key, value
                ),
                HeaderInvalidType,
            )
            type_ = "String"
        if "Type" in self.mapping and type_ not in INFO_TYPES:
            warnings.warn(
                (
                    "Invalid INFO value type {} in header line, using String " "instead, {}={}"
                ).format(self.mapping["Type"], key, value),
                HeaderInvalidType,
            )
            type_ = "String"
        #: value type
        self.type = type_
        # check for "Description" key
        if "Description" not in self.mapping:
            warnings.warn(
                'Field "Description" not found in header line {}={}'.format(key, value),
                HeaderMissingDescription,
            )
        #: description, should be given, ``None`` if not given
        self.description = self.mapping.get("Description")
        #: source of INFO field, ``None`` if not given
        self.source = self.mapping.get("Source")
        #: version of INFO field, ``None`` if not given
        self.version = self.mapping.get("Version")

    def __hash__(self):
        raise TypeError("Unhashable type: InfoHeaderLine")

    def __str__(self):
        return "InfoHeaderLine({}, {}, {})".format(*map(repr, (self.key, self.value, self.mapping)))


class FormatHeaderLine(CompoundHeaderLine):
    """Header line for FORMAT fields
    """

    @classmethod
    def from_mapping(klass, mapping):
        """Construct from mapping, not requiring the string value"""
        return FormatHeaderLine("FORMAT", mapping_to_str(mapping), mapping)

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: key in the INFO field
        self.id = self.mapping["ID"]
        # check for "Number" field
        self.number = self.mapping["Number"]
        # check for "Type" field
        type_ = self.mapping.get("Type")
        if "Type" not in self.mapping:
            warnings.warn(
                ('Field "Type" not found in header line, using String ' "instead {}={}").format(
                    key, value
                ),
                HeaderInvalidType,
            )
            type_ = "String"
        if "Type" in self.mapping and type_ not in FORMAT_TYPES:
            warnings.warn(
                (
                    "Invalid FORMAT value type {} in header line, using String " "instead, {}={}"
                ).format(self.mapping["Type"], key, value),
                HeaderInvalidType,
            )
            type_ = "String"
        #: value type
        self.type = type_
        # check for "Description" key
        if "Description" not in self.mapping:
            warnings.warn(
                'Field "Description" not found in header line {}={}'.format(key, value),
                HeaderMissingDescription,
            )
        #: description, should be given, ``None`` if not given
        self.description = self.mapping.get("Description")
        #: source of INFO field, ``None`` if not given
        self.source = self.mapping.get("Source")
        #: version of INFO field, ``None`` if not given
        self.version = self.mapping.get("Version")

    def __hash__(self):
        raise TypeError("Unhashable type: FormatHeaderLine")

    def __str__(self):
        return "FormatHeaderLine({}, {}, {})".format(
            *map(repr, (self.key, self.value, self.mapping))
        )


class SamplesInfos:
    """Helper class for handling the samples in VCF files

    The purpose of this class is to decouple the sample name list somewhat
    from :py:class:`Header`.  This encapsulates subsetting samples for which
    the genotype should be parsed and reordering samples into output files.

    Note that when subsetting is used and the records are to be written out
    again then the ``FORMAT`` field must not be touched.
    """

    def __init__(self, sample_names, parsed_samples=None):
        #: list of sample that are read from/written to the VCF file at
        #: hand in the given order
        self.names = list(sample_names)
        #: ``set`` with the samples for which the genotype call fields should
        #: be read; can be used for partial parsing (speedup) and defaults
        #: to the full list of samples, None if all are parsed
        self.parsed_samples = parsed_samples
        if self.parsed_samples:
            self.parsed_samples = set(self.parsed_samples)
            assert self.parsed_samples <= set(self.names), "Must be subset!"
        #: mapping from sample name to index
        self.name_to_idx = dict([(name, idx) for idx, name in enumerate(self.names)])

    def copy(self):
        """Return a copy of the object"""
        return SamplesInfos(self.names)

    def is_parsed(self, name):
        """Return whether the sample name is parsed"""
        return (not self.parsed_samples) or name in self.parsed_samples

    def __hash__(self):
        raise TypeError("Unhashable type: SamplesInfos")

    def __str__(self):
        tpl = "SamplesInfos(names={}, name_to_idx={})"
        return tpl.format(self.names, pprint.pformat(self.name_to_idx, width=10 ** 10))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.names == other.names
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return self.names != other.names
        return NotImplemented
