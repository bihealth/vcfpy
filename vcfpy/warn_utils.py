# -*- coding: utf-8 -*-
"""Code for printing warnings"""

import sys

try:
    from cyordereddict import OrderedDict
except ImportError:
    from collections import OrderedDict

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'


class WarningHelper:
    """Helper class for checkers

    This class implements a "warn_once" function that allows to print warnings
    only once and a "print_summary" function that, in the end, allows to print
    a summary table with number of warnings.
    """

    def __init__(self, prefix='[vcfpy] ', stream=sys.stderr):
        #: string to prepend before all warnings
        self.prefix = prefix
        #: the stream to write warnings to
        self.stream = stream
        #: mapping from warning string to counter
        self.warning_counter = OrderedDict()

    def warn_once(self, message):
        """Warn once with message"""
        if message in self.warning_counter:
            self.warning_counter[message] += 1
        else:
            print('{}{}'.format(self.prefix, message), file=self.stream)
            print('(Subsequent identical messages will not be printed)',
                  file=self.stream)
            self.warning_counter[message] = 1

    def print_summary(self, title='WARNINGS', format_='{: 6}\t{}'):
        """Print warning messages and count to ``self.stream``"""
        print('{}\n'.format(title), file=self.stream)
        for msg, count in self.warning_counter.items():
            print(format_.format(count, msg), file=self.stream)
