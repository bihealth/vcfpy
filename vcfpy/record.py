# -*- coding: utf-8 -*-
"""Code for representing a VCF record

The VCF record structure is modeled after the one of PyVCF
"""

import re


#: Code for single nucleotide variant allele
SNV = 'SNV'
#: Code for a multi nucleotide variant allele
MNV = 'MNV'
#: Code for "clean" deletion allele
DEL = 'DEL'
#: Code for "clean" insertion allele
INS = 'INS'
#: Code for indel allele, includes substitutions of unequal length
INDEL = 'INDEL'
#: Code for structural variant allele
SV = 'SV'
#: Code for break-end allele
BND = 'BND'
#: Code for symbolic allele that is neither SV nor BND
SYMBOLIC = 'SYMBOLIC'

#: Code for mixed variant type
MIXED = 'MIXED'

#: Code for homozygous reference
HOM_REF = 0
#: Code for heterozygous
HET = 1
#: Code for homozygous alternative
HOM_ALT = 2

#: Codes for structural variants
SV_CODES = ('DEL', 'INS', 'DUP', 'INV', 'CNV')

#: Characters reserved in VCF, have to be escaped
RESERVED_CHARS = ':;=%,\r\n\t'
#: Mapping for escaping reserved characters
ESCAPE_MAPPING = [
    ('%', '%25'), (':', '%3A'), (';', '%3B'), ('=', '%3D'), (',', '%2C'),
    ('\r', '%0D'), ('\n', '%0A'), ('\t', '%09')
]
#: Mapping from escaped characters to reserved one
UNESCAPE_MAPPING = [(v, k) for k, v in ESCAPE_MAPPING]


class Record:
    """Represent one record from the VCF file

    Record objects are iterators of their calls
    """

    def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
                 calls):
        #: A ``str`` with the chromosome name
        self.CHROM = CHROM
        #: An ``int`` with a 1-based begin position
        self.POS = POS
        #: An ``int`` with a 0-based begin position
        self.begin = POS - 1
        #: An ``int`` with a 0-based end position
        self.end = None  # XXX
        #: A list of the semicolon-separated values of the ID column
        self.ID = list(ID)
        #: A ``str`` with the REF value
        self.REF = REF
        #: A list of alternative allele records of type :py:class:`AltRecord`
        self.ALT = list(ALT)
        #: The quality value, can be ``None``
        self.QUAL = QUAL
        #: A list of strings for the FILTER column
        self.FILTER = FILTER
        #: An OrderedDict giving the values of the INFO column, flags are
        #: mapped to ``True``
        self.INFO = INFO
        #: A list of strings for the FORMAT column
        self.FORMAT = FORMAT
        #: A list of genotype :py:class:`Call` objects
        self.calls = list(calls)
        for call in self.calls:
            call.site = self
        #: A mapping from sample name to entry in self.calls
        self.call_for_sample = {call.sample: call for call in self.calls}

    def is_snv(self):
        """Return ``True`` if it is a SNV"""
        return (len(self.REF) == 1 and all(a.type == 'SNV' for a in self.ALT))

    @property
    def affected_start(self):
        """Return affected start position in 0-based coordinates

        For SNVs, MNVs, and deletions, the behaviour is the start position.
        In the case of insertions, the position behind the insert position is
        returned, yielding a 0-length interval together with
        :py:method:`affected_end`
        """
        types = {alt.type for alt in self.ALT}  # set!
        BAD_MIX = {INS, SV, BND, SYMBOLIC}  # don't mix well with others
        if (BAD_MIX & types) and len(types) == 1 and list(types)[0] == INS:
            # Only insertions, return 0-based position right of first base
            return self.POS  # right of first base
        else:  # Return 0-based start position of first REF base
            return (self.POS - 1)  # left of first base

    @property
    def affected_end(self):
        """Return affected start position in 0-based coordinates

        For SNVs, MNVs, and deletions, the behaviour is based on the start
        position and the length of the REF.  In the case of insertions, the
        position behind the insert position is returned, yielding a 0-length
        interval together with :py:method:`affected_start`
        """
        types = {alt.type for alt in self.ALT}  # set!
        BAD_MIX = {INS, SV, BND, SYMBOLIC}  # don't mix well with others
        if (BAD_MIX & types) and len(types) == 1 and list(types)[0] == INS:
            # Only insertions, return 0-based position right of first base
            return self.POS  # right of first base
        else:  # Return 0-based end position, behind last REF base
            return (self.POS - 1) + len(self.REF)

    def add_filter(self, label):
        """Add label to FILTER if not set yet, removing ``PASS`` entry if
        present
        """
        if label not in self.FILTER:
            if 'PASS' in self.FILTER:
                self.FILTER = [f for f in self.FILTER if f != 'PASS']
            self.FILTER.append(label)

    def add_format(self, key, value=None):
        """Add an entry to format

        The record's calls ``data[key]`` will be set to ``value`` if not yet
        set and value is not ``None``.  If key is already in FORMAT then
        nothing is done.
        """
        if key in self.FORMAT:
            return
        self.FORMAT.append(key)
        if value is not None:
            for call in self:
                call.data.setdefault(key, value)

    def __iter__(self):
        """Return generator yielding from ``self.calls``"""
        yield from self.calls

    def __str__(self):
        tpl = 'Record({})'
        lst = [self.CHROM, self.POS, self.ID, self.REF, self.ALT, self.QUAL,
               self.FILTER, self.INFO, self.FORMAT, self.calls]
        return tpl.format(', '.join(map(repr, lst)))

    def __repr__(self):
        return str(self)


ALLELE_DELIM = re.compile(r'[|/]')


class Call:
    """The information for a genotype callable

    By VCF, this should always include the genotype information and
    can contain an arbitrary number of further annotation, e.g., the
    coverage at the variant position.
    """

    def __init__(self, sample, data, site=None):
        #: the name of the sample for which the call was made
        self.sample = sample
        #: an OrderedDict with the key/value pair information from the
        #: call's data
        self.data = data
        #: the :py:class:`Record` of this :py:class:`Call`
        self.site = site
        #: the allele numbers (0, 1, ...) in this calls or None for no-call
        self.gt_alleles = None
        #: whether or not the variant is fully called
        self.called = None
        #: the number of alleles in this sample's call
        self.plodity = None
        if self.data.get('GT', None) is not None:
            self.gt_alleles = []
            for allele in ALLELE_DELIM.split(self.data['GT']):
                if allele == '.':
                    self.gt_alleles.append(allele)
                else:
                    self.gt_alleles.append(int(allele))
            self.called = all([al is not None for al in self.gt_alleles])
            self.ploidty = len(self.gt_alleles)

    @property
    def is_phased(self):
        """Return boolean indicating whether this call is phased"""
        return '|' in self.data.get('GT', '')

    @property
    def gt_phase_char(self):
        """Return character to use for phasing"""
        return '/' if not self.is_phased else '|'

    @property
    def gt_bases(self):
        """Return the actual genotype bases, e.g. if VCF genotype is 0/1,
        could return ('A', 'T')
        """
        result = []
        for a in self.gt_alleles:
            if a is None:
                result.append(None)
            elif a == 0:
                result.append(self.site.REF)
            else:
                result.append(self.site.ALT[a - 1].value)
        return tuple(result)

    @property
    def gt_type(self):
        """The type of genotype, returns one of ``HOM_REF``, ``HOM_ALT``, and
        ``HET``.
        """
        if not self.called:
            return None  # not called
        elif all(a == 0 for a in self.gt_alleles):
            return HOM_REF
        elif len(set(self.gt_alleles)) == 1:
            return HOM_ALT
        else:
            return HET

    def is_filtered(self, require=None, ignore=['PASS']):
        """Return ``True`` for filtered calls

        :param iterable ignore: if set, the filters to ignore, make sure to
            include 'PASS', when setting
        :param iterable require: if set, the filters to require for returning
            ``True``
        """
        if 'FT' not in self.data or not self.data['FT']:
            return False
        for ft in self.data['FT']:
            if ft in ignore:
                continue  # skip
            if not require:
                return True
            elif ft in require:
                return True
        return False

    @property
    def is_het(self):
        """Return ``True`` for heterozygous calls"""
        return self.gt_type == HET

    @property
    def is_variant(self):
        """Return ``True`` for non-hom-ref calls"""
        return self.gt_type != HOM_REF

    def __str__(self):
        tpl = 'Call({})'
        lst = [self.sample, self.data]
        return tpl.format(', '.join(map(repr, lst)))

    def __repr__(self):
        return str(self)


class AltRecord:
    """An alternative allele Record

    Currently, can be a substitution, an SV placeholder, or breakend
    """

    def __init__(self, type_=None):
        #: String describing the type of the variant, could be one of
        #: SNV, MNV, could be any of teh types described in the ALT
        #: header lines, such as DUP, DEL, INS, ...
        self.type = type_


class Substitution(AltRecord):
    """A basic alternative allele record describing a REF->AltRecord
    substitution

    Note that this subsumes MNVs, insertions, and deletions.
    """

    def __init__(self, type_, value):
        super().__init__(type_)
        #: The alternative base sequence to use in the substitution
        self.value = value

    def __str__(self):
        tpl = 'Substitution(type_={}, value={})'
        return tpl.format(*map(repr, [self.type, self.value]))

    def __repr__(self):
        return str(self)


class SV(AltRecord):
    """A placeholder for an SV

    For SVs, simply the ``type`` is written out as ``"<type>"`` as alternative
    """

    def __init__(self, type_, value):
        super().__init__(type_)
        #: The alternative base sequence to use in the substitution
        self.value = value

    def __str__(self):
        tpl = 'Substitution(type_={}, value={})'
        return tpl.format(*map(repr, [self.type, self.value]))

    def __repr__(self):
        return str(self)


class BreakEnd(AltRecord):
    """A placeholder for a breakend"""

    def __init__(self, type_, value):
        super().__init__(type_)
        #: The alternative base sequence to use in the substitution
        self.value = value

    def __str__(self):
        tpl = 'Substitution(type_={}, value={})'
        return tpl.format(*map(repr, [self.type, self.value]))

    def __repr__(self):
        return str(self)


class SingleBreakEnd(AltRecord):
    """A placeholder for a single breakend"""

    def __init__(self, type_, value):
        super().__init__(type_)
        #: The alternative base sequence to use in the substitution
        self.value = value

    def __str__(self):
        tpl = 'Substitution(type_={}, value={})'
        return tpl.format(*map(repr, [self.type, self.value]))

    def __repr__(self):
        return str(self)


class SymbolicAllele(AltRecord):
    """A placeholder for a symbolic allele"""

    def __init__(self, type_, value):
        super().__init__(type_)
        #: The alternative base sequence to use in the substitution
        self.value = value

    def __str__(self):
        tpl = 'Substitution(type_={}, value={})'
        return tpl.format(*map(repr, [self.type, self.value]))

    def __repr__(self):
        return str(self)
