# -*- coding: utf-8 -*-

try:
    from cyordereddict import OrderedDict
except ImportError:
    from collections import OrderedDict

from .exceptions import VCFPyException, InvalidHeaderException, \
    InvalidRecordException, IncorrectVCFFormat, HeaderNotFound

from .header import Header, HeaderLine, SimpleHeaderFile, \
    SimpleHeaderFile, ContigHeaderLine, FilterHeaderLine, \
    CompoundHeaderLine, InfoHeaderLine, FormatHeaderLine, \
    SamplesInfos

from .record import Record, Call, AltRecord, Substitution, SV, BreakEnd, \
    SymbolicAllele

from .reader import Reader

from .writer import Writer

__version__ = '0.2.0'
