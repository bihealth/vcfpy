# -*- coding: utf-8 -*-
"""Writing of VCF files to ``file``-like objects

Currently, only writing to plain-text files is supported
"""

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# TODO: use context manager for making close-able
# TODO: allow some configuration to make warning into exceptions?


class VCFWriter:
    """Class for writing VCF files to ``file``-like objects

    Instead of using the constructor, use the class methods
    :py:meth:`~VCFWriter.from_file` and
    :py:meth:`~VCFWriter.from_path`.

    The writer has to be constructed with a :py:class:`~vcfpy.header.VCFHeader`
    and a :py:class:`~vcfpy.header.SamplesInfos` object and the full VCF
    header will be written immediately on construction.  This, of course,
    implies that modifying the header after construction is illegal.
    """

    @classmethod
    def from_file(klass, stream, path=None):
        """Create new :py:class:`VCFWriter` from file

        :param stream: ``file``-like object to write to
        :param path: optional string with path to store (for display only)
        """
        return VCFWriter(stream=stream, path=path)

    @classmethod
    def from_path(klass, path):
        """Create new :py:class:`VCFWriter` from path

        :param path: the path to load from (converted to ``str`` for
            compatibility with ``path.py``)
        """
        path = str(path)
        if path.endswith('.gz'):
            raise NotImplementedError('Writing to bgzf not supported')
        else:
            f = open(path, 'wt')
        return klass.from_file(stream=f, path=path)

    def __init__(self, header, samples, stream, path=None):
        #: the :py:class:~vcfpy.header.VCFHeader` written out
        self.header = header
        #: the :py:class:~vcfpy.header.SamplesInfos` written out
        self.samples = samples
        #: stream (``file``-like object) to read from
        self.stream = stream
        #: optional ``str`` with the path to the stream
        self.path = path

    def write_record(self, record):
        """Write out the given :py:class:`vcfpy.record.VCFRecord` to this
        VCFWriter"""
        self._serialize_record(record)

    def _serialize_record(self, record):
        """Serialize whole VCFRecord"""
        f = self._empty_to_dot
        row = [record.CHROM, record.POS]
        row += f(';'.join(record.ID))
        row.append(f(record.REF_))
        row += [f(a.value) for a in record.ALT]
        row.append(f(record.QUAL))
        row.append(f(';'.join(record.FILTER)))
        row.append(f(self._serialize_info(record)))
        row.append(':'.join(record.FORMAT))
        row += [self._serialize_call(record.format, c) for c in record.calls]
        print(*row, sep='\t', file=self.stream)

    def _serialize_info(self, record):
        """Return serialized version of record.INFO"""
        raise NotImplementedError('Implement me!')

    def _serialize_call(self, format, call):
        """Return serialized version of the Call using the record's FORMAT'"""
        raise NotImplementedError('Implement me!')

    def _empty_to_dot(self, val):
        """Return val or '.' if empty value"""
        if val == '' or val is None or val == []:
            return '.'
        else:
            return val
