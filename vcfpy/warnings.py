# -*- coding: utf-8 -*-
"""Warnings for the vcfpy module
"""

# TODO: improve names
# TODO: add "once" filter

class VCFPyWarning(Warning):
    """Base class for module's warnings"""


class DuplicateHeaderLineWarning(VCFPyWarning):
    """A header line occurs twice in a header"""


class FieldInfoNotFound(VCFPyWarning):
    """A header field is not found, default is used"""


class FieldMissingNumber(VCFPyWarning):
    """Raised when compound heade misses number"""


class FieldInvalidNumber(VCFPyWarning):
    """Raised when compound header has invalid number"""


class HeaderInvalidType(VCFPyWarning):
    """Raised when compound header has invalid type"""


class HeaderMissingDescription(VCFPyWarning):
    """Raised when compound header has missing description"""


class LeadingTrailingSpaceInKey(VCFPyWarning):
    """Leading or trailing space in key"""


class UnknownFilter(VCFPyWarning):
    """Missing FILTER"""


class UnknownVCFVersion(VCFPyWarning):
    """Unknown VCF version"""


class IncorrectListLength(VCFPyWarning):
    """Wrong length of multi-element field"""


class SpaceInChromLine(VCFPyWarning):
    """Space instead of TAB in ##CHROM line"""
