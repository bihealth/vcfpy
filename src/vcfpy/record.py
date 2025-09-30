# -*- coding: utf-8 -*-
"""Code for representing a VCF record

The VCF record structure is modeled after the one of PyVCF
"""

import re
import warnings
from typing import Any, Iterator, Literal, Sequence

from vcfpy.exceptions import CannotModifyUnparsedCallWarning

#: Code for single nucleotide variant allele
SNV = "SNV"
#: Code for a multi nucleotide variant allele
MNV = "MNV"
#: Code for "clean" deletion allele
DEL = "DEL"
#: Code for "clean" insertion allele
INS = "INS"
#: Code for indel allele, includes substitutions of unequal length
INDEL = "INDEL"
#: Code for structural variant allele
SV = "SV"
#: Code for break-end allele
BND = "BND"
#: Code for symbolic allele that is neither SV nor BND
SYMBOLIC = "SYMBOLIC"

#: Code for mixed variant type
MIXED = "MIXED"

#: Code for homozygous reference
HOM_REF = 0
#: Code for heterozygous
HET = 1
#: Code for homozygous alternative
HOM_ALT = 2

#: Characters reserved in VCF, have to be escaped in INFO fields
RESERVED_CHARS: dict[Literal["INFO", "FORMAT"], str] = {"INFO": ";=%,\r\n\t", "FORMAT": ":=%,\r\n\t"}
#: Mapping for escaping reserved characters
ESCAPE_MAPPING = [
    ("%", "%25"),
    (":", "%3A"),
    (";", "%3B"),
    ("=", "%3D"),
    (",", "%2C"),
    ("\r", "%0D"),
    ("\n", "%0A"),
    ("\t", "%09"),
]
#: Mapping from escaped characters to reserved one
UNESCAPE_MAPPING = [(v, k) for k, v in ESCAPE_MAPPING]


class Record:
    """Represent one record from the VCF file

    Record objects are iterators of their calls
    """

    def __init__(
        self,
        CHROM: str,
        POS: int,
        ID: list[str],
        REF: str,
        ALT: list["AltRecord"],
        QUAL: float | None,
        FILTER: list[str],
        INFO: dict[str, Any],
        FORMAT: list[str] | None = None,
        calls: Sequence["Call | UnparsedCall"] | None = None,
    ):
        if bool(FORMAT) != bool(calls):
            raise ValueError("Either provide both FORMAT and calls or none.")
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
        self.ALT: list["AltRecord"] = list(ALT)
        #: The quality value, can be ``None``
        self.QUAL = QUAL
        #: A list of strings for the FILTER column
        self.FILTER = FILTER
        #: An OrderedDict giving the values of the INFO column, flags are
        #: mapped to ``True``
        self.INFO = INFO
        #: A list of strings for the FORMAT column.  Optional, must be given if
        #: and only if ``calls`` is also given.
        self.FORMAT = FORMAT or []
        #: A list of genotype :py:class:`Call` objects.  Optional, must be given if
        #: and only if ``FORMAT`` is also given.
        self.calls = calls or []
        #: A mapping from sample name to entry in self.calls.
        self.call_for_sample = {}
        self.update_calls(self.calls)

    def update_calls(self, calls: Sequence["Call | UnparsedCall"]):
        """Update ``self.calls`` and other fields as necessary."""
        for call in calls:
            call.site = self
        self.call_for_sample = {call.sample: call for call in calls}

    def is_snv(self):
        """Return ``True`` if it is a SNV"""
        return len(self.REF) == 1 and all(a.type == "SNV" for a in self.ALT)

    @property
    def affected_start(self):
        """Return affected start position in 0-based coordinates

        For SNVs, MNVs, and deletions, the behaviour is the start position.
        In the case of insertions, the position behind the insert position is
        returned, yielding a 0-length interval together with
        :py:meth:`~Record.affected_end`
        """
        types = {alt.type for alt in self.ALT}  # set!
        BAD_MIX = {INS, SV, BND, SYMBOLIC}  # don't mix well with others
        if (BAD_MIX & types) and len(types) == 1 and list(types)[0] == INS:
            # Only insertions, return 0-based position right of first base
            return self.POS  # right of first base
        else:  # Return 0-based start position of first REF base
            return self.POS - 1  # left of first base

    @property
    def affected_end(self):
        """Return affected start position in 0-based coordinates

        For SNVs, MNVs, and deletions, the behaviour is based on the start
        position and the length of the REF.  In the case of insertions, the
        position behind the insert position is returned, yielding a 0-length
        interval together with :py:meth:`~Record.affected_start`
        """
        types = {alt.type for alt in self.ALT}  # set!
        BAD_MIX = {INS, SV, BND, SYMBOLIC}  # don't mix well with others
        if (BAD_MIX & types) and len(types) == 1 and list(types)[0] == INS:
            # Only insertions, return 0-based position right of first base
            return self.POS  # right of first base
        else:  # Return 0-based end position, behind last REF base
            return (self.POS - 1) + len(self.REF)

    def add_filter(self, label: str):
        """Add label to FILTER if not set yet, removing ``PASS`` entry if
        present
        """
        if label not in self.FILTER:
            if "PASS" in self.FILTER:
                self.FILTER = [f for f in self.FILTER if f != "PASS"]
            self.FILTER.append(label)

    def add_format(self, key: str, value: Any | None = None):
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
                if isinstance(call, UnparsedCall):
                    warnings.warn("UnparsedCall encountered, skipping", CannotModifyUnparsedCallWarning)
                else:
                    call.data.setdefault(key, value)

    def __iter__(self) -> Iterator["Call | UnparsedCall"]:
        """Return generator yielding from ``self.calls``"""
        yield from self.calls

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        tpl = "Record({})"
        lst: list[Any] = [
            self.CHROM,
            self.POS,
            self.ID,
            self.REF,
            self.ALT,
            self.QUAL,
            self.FILTER,
            self.INFO,
            self.FORMAT,
            self.calls,
        ]
        return tpl.format(", ".join(map(repr, lst)))

    def __repr__(self):
        return str(self)


class UnparsedCall:
    """Placeholder for :py:class:`Call` when parsing only a subset of fields"""

    def __init__(self, sample: str, unparsed_data: Any, site: Record | None = None):
        #: the name of the sample for which the call was made
        self.sample = sample
        #: ``str`` with the unparsed data
        self.unparsed_data = unparsed_data
        #: the :py:class:`Record` of this :py:class:`Call`
        self.site = site


#: Regular expression for splitting alleles
ALLELE_DELIM = re.compile(r"[|/]")


class Call:
    """The information for a genotype callable

    By VCF, this should always include the genotype information and
    can contain an arbitrary number of further annotation, e.g., the
    coverage at the variant position.
    """

    def __init__(self, sample: str, data: dict[str, Any], site: Record | None = None):
        #: the name of the sample for which the call was made
        self.sample = sample
        #: an OrderedDict with the key/value pair information from the
        #: call's data
        self.data: dict[str, Any] = data
        #: the :py:class:`Record` of this :py:class:`Call`
        self.site = site
        #: the allele numbers (0, 1, ...) in this calls or None for no-call
        self.gt_alleles: list[int | None] | None = None
        #: whether or not the variant is fully called
        self.called = None
        #: the number of alleles in this sample's call
        self.ploidy = None
        self._genotype_updated()

    def set_genotype(self, genotype: str | None):
        """Set ``self.data["GT"]`` to ``genotype`` and properly update related
        properties.
        """
        self.data["GT"] = genotype
        self._genotype_updated()

    def _genotype_updated(self):
        """Update fields related to ``self.data["GT"]``."""
        if self.data.get("GT", None) is None:
            self.gt_alleles = None
            self.called = None
            self.ploidy = None
        else:
            self.gt_alleles = []
            for allele in ALLELE_DELIM.split(str(self.data["GT"])):
                if allele == ".":
                    self.gt_alleles.append(None)
                else:
                    self.gt_alleles.append(int(allele))
            self.called = all(al is not None for al in self.gt_alleles)
            self.ploidy = len(self.gt_alleles)

    @property
    def is_phased(self):
        """Return boolean indicating whether this call is phased"""
        return "|" in self.data.get("GT", "")

    @property
    def gt_phase_char(self):
        """Return character to use for phasing"""
        return "/" if not self.is_phased else "|"

    @property
    def gt_bases(self) -> tuple[str | None, ...]:
        """Return the actual genotype bases, e.g. if VCF genotype is 0/1,
        could return ('A', 'T')
        """
        result: list[str | None] = []
        for a in self.gt_alleles or []:
            if a is None:
                result.append(None)
            elif a == 0:
                if self.site is None:
                    raise ValueError("Cannot determine bases without site being set")
                result.append(self.site.REF)
            else:
                if self.site is None:  # pragma: no cover
                    raise ValueError("Cannot determine bases without site being set")
                result.append(getattr(self.site.ALT[a - 1], "value", None))
        return tuple(result)

    @property
    def gt_type(self) -> Literal[0, 1, 2] | None:
        """The type of genotype, returns one of ``HOM_REF``, ``HOM_ALT``, and
        ``HET``.
        """
        if not self.called or not self.gt_alleles:
            return None  # not called
        elif all(a == 0 for a in self.gt_alleles):
            return HOM_REF
        elif len(set(self.gt_alleles)) == 1:
            return HOM_ALT
        else:
            return HET

    def is_filtered(self, require: list[str] | None = None, ignore: list[str] | None = None):
        """Return ``True`` for filtered calls

        :param iterable ignore: if set, the filters to ignore, make sure to
            include 'PASS', when setting, default is ``['PASS']``
        :param iterable require: if set, the filters to require for returning
            ``True``
        """
        ignore = ignore or ["PASS"]
        if "FT" not in self.data or not self.data["FT"]:
            return False
        for ft in self.data["FT"]:
            if ft in ignore:
                continue  # skip
            if not require:
                return True
            elif ft in require:
                return True
        return False

    @property
    def is_het(self) -> bool:
        """Return ``True`` for heterozygous calls"""
        return self.gt_type == HET

    @property
    def is_variant(self) -> bool:
        """Return ``True`` for non-hom-ref calls"""
        return bool(self.gt_type)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self) -> int:
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        tpl = "Call({})"
        lst: list[str | dict[str, Any]] = [self.sample, self.data]
        return tpl.format(", ".join(map(repr, lst)))

    def __repr__(self):
        return str(self)


class AltRecord:
    """An alternative allele Record

    Currently, can be a substitution, an SV placeholder, or breakend
    """

    def __init__(
        self, type_: Literal["SNV", "MNV", "DEL", "INS", "INDEL", "SV", "BND", "SYMBOLIC", "MIXED"] | None = None
    ):
        #: String describing the type of the variant, could be one of
        #: SNV, MNV, could be any of teh types described in the ALT
        #: header lines, such as DUP, DEL, INS, ...
        self.type = type_

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self) -> int:
        return hash(tuple(sorted(self.__dict__.items())))

    def serialize(self) -> str:
        """Return ``str`` with representation for VCF file"""
        raise NotImplementedError("Abstract class, implemented in sub class")


class Substitution(AltRecord):
    """A basic alternative allele record describing a REF->AltRecord
    substitution

    Note that this subsumes MNVs, insertions, and deletions.
    """

    def __init__(
        self, type_: Literal["SNV", "MNV", "DEL", "INS", "INDEL", "SV", "BND", "SYMBOLIC", "MIXED"], value: str
    ):
        super().__init__(type_)
        #: The alternative base sequence to use in the substitution
        self.value = value

    def serialize(self) -> str:
        return self.value

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self) -> int:
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self) -> str:
        tpl = "Substitution(type_={}, value={})"
        return tpl.format(*map(repr, [self.type, self.value]))

    def __repr__(self):
        return str(self)


#: code for five prime orientation :py:class:`BreakEnd`
FIVE_PRIME = "5"
#: code for three prime orientation :py:class:`BreakEnd`
THREE_PRIME = "3"

#: code for forward orientation
FORWARD = "+"
#: code for reverse orientation
REVERSE = "-"


class BreakEnd(AltRecord):
    """A placeholder for a breakend"""

    def __init__(
        self,
        mate_chrom: str | None,
        mate_pos: int | None,
        orientation: str | None,
        mate_orientation: Literal["+", "-"] | None,
        sequence: str,
        within_main_assembly: bool | None,
    ):
        super().__init__("BND")
        #: chromosome of the mate breakend
        self.mate_chrom = mate_chrom
        #: position of the mate breakend
        self.mate_pos = mate_pos
        #: orientation of this breakend
        self.orientation = orientation
        #: orientation breakend's mate
        self.mate_orientation = mate_orientation
        #: breakpoint's connecting sequence
        self.sequence = sequence
        #: ``bool`` specifying if the breakend mate is within the assembly
        #: (``True``) or in an ancillary assembly (``False``)
        self.within_main_assembly = within_main_assembly

    def serialize(self):
        """Return string representation for VCF"""
        if self.mate_chrom is None:
            remote_tag = "."
        else:
            if self.within_main_assembly:
                mate_chrom = self.mate_chrom
            else:
                mate_chrom = "<{}>".format(self.mate_chrom)
            if self.mate_orientation is None:  # pragma: no cover
                raise ValueError("mate_orientation must be set if mate_chrom is set")
            tpl = {FORWARD: "[{}:{}[", REVERSE: "]{}:{}]"}[self.mate_orientation]
            remote_tag = tpl.format(mate_chrom, self.mate_pos)
        if self.orientation == FORWARD:
            return remote_tag + self.sequence
        else:
            return self.sequence + remote_tag

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        tpl = "BreakEnd({})"
        vals: list[Any] = [
            self.mate_chrom,
            self.mate_pos,
            self.orientation,
            self.mate_orientation,
            self.sequence,
            self.within_main_assembly,
        ]
        return tpl.format(", ".join(map(repr, vals)))

    def __repr__(self):
        return str(self)


class SingleBreakEnd(BreakEnd):
    """A placeholder for a single breakend"""

    def __init__(self, orientation: str, sequence: str):
        super().__init__(None, None, orientation, None, sequence, None)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        tpl = "SingleBreakEnd({})"
        vals: list[Any] = [self.orientation, self.sequence]
        return tpl.format(", ".join(map(repr, vals)))


class SymbolicAllele(AltRecord):
    """A placeholder for a symbolic allele

    The allele symbol must be defined in the header using an ``ALT`` header
    before being parsed.  Usually, this is used for succinct descriptions of
    structural variants or IUPAC parameters.
    """

    def __init__(self, value: str):
        super().__init__(SYMBOLIC)
        #: The symbolic value, e.g. 'DUP'
        self.value = value

    def serialize(self):
        return "<{}>".format(self.value)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        raise NotImplementedError  # pragma: no cover

    def __ne__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return not self.__eq__(other)
        raise NotImplementedError  # pragma: no cover

    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    def __str__(self):
        return "SymbolicAllele({})".format(repr(self.value))

    def __repr__(self):
        return str(self)
