# -*- coding: utf-8 -*-

from .compat import OrderedDict
from .exceptions import (
    HeaderNotFound,
    IncorrectVCFFormat,
    InvalidHeaderException,
    InvalidRecordException,
    VCFPyException,
)
from .header import (
    AltAlleleHeaderLine,
    CompoundHeaderLine,
    ContigHeaderLine,
    FieldInfo,
    FilterHeaderLine,
    FormatHeaderLine,
    Header,
    HeaderLine,
    InfoHeaderLine,
    MetaHeaderLine,
    PedigreeHeaderLine,
    SampleHeaderLine,
    SamplesInfos,
    SimpleHeaderLine,
    header_without_lines,
)
from .reader import Reader
from .record import (
    BND,
    DEL,
    FIVE_PRIME,
    FORWARD,
    HET,
    HOM_ALT,
    HOM_REF,
    INDEL,
    INS,
    MIXED,
    MNV,
    REVERSE,
    SNV,
    SV,
    SYMBOLIC,
    THREE_PRIME,
    AltRecord,
    BreakEnd,
    Call,
    Record,
    SingleBreakEnd,
    Substitution,
    SymbolicAllele,
    UnparsedCall,
)
from .version import __version__
from .writer import Writer
