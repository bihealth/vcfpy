# -*- coding: utf-8 -*-
"""Code for representing a VCF record

The VCF record structure is modeled after the one of PyVCF
"""

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

#: Codes for structural variants
SV_CODES = ('DEL', 'INS', 'DUP', 'INV', 'CNV')


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

    @property
    def affected_start(self):
        """Return affected start position in 0-based coordinates
        
        For SNVs, MNVs, and deletions, the behaviour is the start position.
        In the case of insertions, the position behind the insert position is
        returned, yielding a 0-length interval together with
        :py:method:`affected_end`
        """
        types = {alt.type for alt in self.ALT}  # set!
        BAD_MIX = [INS, SV, BND, SYMBOLIC]  # don't mix well with others
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
        BAD_MIX = [INS, SV, BND, SYMBOLIC]  # don't mix well with others
        if (BAD_MIX & types) and len(types) == 1 and list(types)[0] == INS:
            # Only insertions, return 0-based position right of first base
            return self.POS  # right of first base
        else:  # Return 0-based end position, behind last REF base
            return (self.POS - 1) + len(self.REF)

    def add_filter(self, label):
        """Add label to FILTER if not set yet"""
        if label not in self.FILTER:
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

    @property
    def gt_bases(self):
        """Return the actual genotype alleles, e.g. if VCF genotype is 0/1,
        could return A/T"""
        raise NotImplementedError('Implement me!')

    @property
    def gt_type(self):
        """The type of genotype, mapping is

        - hom_ref = 0
        - het = 1
        - hom_alt = 2 (which alt is untracked)
        - uncalled = ``None``
        """
        raise NotImplementedError('Implement me!')

    @property
    def is_filtered(self):
        """Return ``True`` for filtered calls"""
        raise NotImplementedError('Implement me!')

    @property
    def is_het(self):
        """Return ``True`` for filtered calls"""
        raise NotImplementedError('Implement me!')

    @property
    def is_variant(self):
        """Return ``True`` for filtered calls"""
        raise NotImplementedError('Implement me!')

    @property
    def is_phased(self):
        """Return ``True`` for phased calls"""
        raise NotImplementedError('Implement me!')

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
