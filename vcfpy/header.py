# -*- coding: utf-8 -*-
"""Code for representing the VCF header part

The VCF header class structure is modeled after HTSJDK
"""

import collections
import sys

from . import exceptions

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# TODO: ##ALT=
# TODO: ##assembly=
# TODO: ##PEDIGREE=
# TODO: ##META=
# TODO: ##SAMPLE=
# TODO: use @property instead of assignment to self.<property>
# TODO: id_ => id?
# TODO: explicit tests for header module

# Tuples of valid entries -----------------------------------------------------
#
#: valid INFO value types
INFO_TYPES = ('Integer', 'Float', 'Flag', 'Character', 'String')
#: valid FORMAT value types
FORMAT_TYPES = ('Integer', 'Float', 'Character', 'String')
#: valid values for "Number" entries, except for integers
VALID_NUMBERS = ('A', 'R', 'G', '.')

# Constants for "Number" entries ----------------------------------------------
#
#: number of alleles excluding reference
HEADER_NUMBER_ALLELES = 'A'
#: number of alleles including reference
HEADER_NUMBER_REF = 'R'
#: number of genotypes
HEADER_NUMBER_GENOTYPES = 'G'
#: unbounded number of values
HEADER_NUMBER_UNBOUNDED = '.'


def _warn_missing(msg):
    """Print warning message in case of missing attributes"""
    print('[vcfpy] WARNING: {}'.format(msg), file=sys.stderr)


class VCFHeader:
    """Represent header of VCF file

    While this class allows mutating records, it should not be changed once it
    has been assigned to
    """

    def __init__(self, lines=[], samples=None):
        #: ``list`` of :py:VCFHeaderLine objects
        self.lines = list(lines)
        #: :py:class:`SamplesInfo` object
        self.samples = samples

    def __str__(self):
        tpl = 'VCFHeader(lines={}, samples={})'
        return tpl.format(*map(repr, (self.lines, self.samples)))

    def __repr__(self):
        return str(self)


class VCFHeaderLine:
    """Base class for VCF header lines
    """

    def __init__(self, key, value):
        #: ``str`` with key of header line
        self.key = key
        #: ``str`` with raw value of header line
        self.value = value

    def __str__(self):
        return 'VCFHeaderLine({}, {})'.format(
            *map(repr, (self.key, self.value)))

    def __repr__(self):
        return str(self)


class VCFSimpleHeaderLine(VCFHeaderLine):
    """Base class for simple header lines, currently contig and filter
    header lines

    :raises: :py:class:`vcfpy.exceptions.InvalidHeaderException` in
        the case of missing key ``"ID"``
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value)
        # check existence of key "ID"
        if 'ID' not in mapping:
            raise exceptions.InvalidHeaderException(
                'Missing key "ID" in header line "{}={}"'.format(
                    key, value))
        #: ``collections.OrderedDict`` with key/value mapping of the attributes
        self.mapping = collections.OrderedDict(mapping.items())

    def __str__(self):
        return 'VCFSimpleHeaderLine({}, {}, {})'.format(
            *map(repr, (self.key, self.value, self.mapping)))


class VCFContigHeaderLine(VCFSimpleHeaderLine):
    """Contig header line

    Most importantly, parses the ``'length'`` key into an integer
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        # convert 'length' entry to integer if possible
        if 'length' in self.mapping:
            mapping['length'] = int(mapping['length'])
        else:
            _warn_missing(
                'Field "length" not found in header line {}={}'.format(
                    key, value))
        #: name of the contig
        self.id_ = self.mapping['ID']
        #: length of the contig, ``None`` if missing
        self.length = self.mapping.get('length')

    def __str__(self):
        return 'VCFContigHeaderLine({}, {}, {})'.format(
            *map(repr, (self.key, self.value, self.mapping)))


class VCFFilterHeaderLine(VCFSimpleHeaderLine):
    """FILTER header line
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        # check for "Description" key
        if 'Description' not in self.mapping:
            _warn_missing(
                'Field "Description" not found in header line {}={}'.format(
                    key, value))
        #: token for the filter
        self.id_ = self.mapping['ID']
        #: description for the filter, ``None`` if missing
        self.description = self.mapping.get('Description')

    def __str__(self):
        return 'VCFFilterHeaderLine({}, {}, {})'.format(
            *map(repr, (self.key, self.value, self.mapping)))


class VCFCompoundHeaderLine(VCFHeaderLine):
    """Base class for compound header lines, currently format and header lines

    Compound header lines describe fields that can have more than one entry.
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value)
        #: OrderedDict with key/value mapping
        self.mapping = collections.OrderedDict(mapping.items())
        # check that 'Number' is given and use "." otherwise
        if 'Number' not in self.mapping:
            print(('[pyvcf] WARNING: missing number, using '
                   'unbounded/"." instead'), file=sys.stderr)
            self.mapping['Number'] = '.'
        try:
            self.mapping['Number'] = self._parse_number(
                self.mapping['Number'])
        except ValueError:
            print(('[pyvcf] WARNING: invalid number {}, using '
                   'unbounded/"." instead').format(self.mapping['Number']),
                  file=sys.stderr)
            self.mapping['Number'] = '.'

    def _parse_number(self, number):
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

    def __str__(self):
        return 'VCFCompoundHeaderLine({}, {}, {})'.format(
            *map(repr, (self.key, self.value, self.mapping)))


class VCFInfoHeaderLine(VCFCompoundHeaderLine):
    """Header line for INFO fields

    Note that the ``Number`` field will be parsed into an ``int`` if
    possible.  Otherwise, the constants ``HEADER_NUMBER_*`` will be used.
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: key in the INFO field
        self.id_ = self.mapping['ID']
        # check for "Number" field
        self.number = self.mapping['Number']
        # check for "Type" field
        type_ = self.mapping.get('Type')
        if 'Type' not in self.mapping:
            _warn_missing(
                ('Field "Type" not found in header line, using String '
                 'instead {}={}').format(key, value))
            type_ = 'String'
        if 'Type' in self.mapping and type_ not in INFO_TYPES:
            _warn_missing(
                ('Invalid INFO value type {} in header line, using String '
                 'instead, {}={}').format(self.mapping['Type'], key, value))
            type_ = 'String'
        #: value type
        self.type = type_
        # check for "Description" key
        if 'Description' not in self.mapping:
            _warn_missing(
                'Field "Description" not found in header line {}={}'.format(
                    key, value))
        #: description, should be given, ``None`` if not given
        self.description = self.mapping.get('Description')
        #: source of INFO field, ``None`` if not given
        self.source = self.mapping.get('Source')
        #: version of INFO field, ``None`` if not given
        self.version = self.mapping.get('Version')

    def __str__(self):
        return 'VCFInfoHeaderLine({}, {}, {})'.format(
            *map(repr, (self.key, self.value, self.mapping)))


class VCFFormatHeaderLine(VCFCompoundHeaderLine):
    """Header line for FORMAT fields
    """

    def __init__(self, key, value, mapping):
        super().__init__(key, value, mapping)
        #: key in the INFO field
        self.id_ = self.mapping['ID']
        # check for "Number" field
        self.number = self.mapping['Number']
        # check for "Type" field
        type_ = self.mapping.get('Type')
        if 'Type' not in self.mapping:
            _warn_missing(
                ('Field "Type" not found in header line, using String '
                 'instead {}={}').format(key, value))
            type_ = 'String'
        if 'Type' in self.mapping and type_ not in FORMAT_TYPES:
            _warn_missing(
                ('Invalid INFO value type {} in header line, using String '
                 'instead, {}={}').format(self.mapping['Type'], key, value))
            type_ = 'String'
        #: value type
        self.type = type_
        # check for "Description" key
        if 'Description' not in self.mapping:
            _warn_missing(
                'Field "Description" not found in header line {}={}'.format(
                    key, value))
        #: description, should be given, ``None`` if not given
        self.description = self.mapping.get('Description')
        #: source of INFO field, ``None`` if not given
        self.source = self.mapping.get('Source')
        #: version of INFO field, ``None`` if not given
        self.version = self.mapping.get('Version')

    def __str__(self):
        return 'VCFFormatHeaderLine({}, {}, {})'.format(
            *map(repr, (self.key, self.value, self.mapping)))


class SamplesInfos:
    """Helper class for handling and mapping of sample names to numeric indices
    """

    def __init__(self, sample_names):
        #: list of sample names
        self.names = list(sample_names)
        #: mapping from sample name to index
        self.name_to_idx = dict([
            (name, idx) for idx, name in enumerate(self.names)])

    def __str__(self):
        tpl = 'SampleInfo(names={}, name_to_idx={})'
        return tpl.format(self.names, self.name_to_idx)

    def __repr__(self):
        return str(self)