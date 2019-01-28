# -*- coding: utf-8 -*-
"""Code for defining types depending on fast implementations

Namely, defines the type ``OrderedDict``:

- on CPython >=3.6, define as ``dict``
- otherwise, if ``cyordereddict`` is available use
  ``cyordereddict.OrderedDict``
- otherwise, use ``collections.OrderedDict``
"""

import sys


if sys.version_info[:2] >= (3, 6):
    OrderedDict = dict
else:
    try:
        from cyordereddict import OrderedDict
    except ImportError:
        from collections import OrderedDict  # noqa: ignore=F401
