# -*- coding: utf-8 -*-
"""Writing of VCF files to ``file``-like objects

Currently, only writing to plain-text files is supported
"""

from . import parser
from . import record
from . import bgzf

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


def format_atomic(value, section):
    """Format atomic value

    This function also takes care of escaping the value in case one of the
    reserved characters occurs in the value.
    """
    # Perform escaping
    if isinstance(value, str):
        if any(r in value for r in record.RESERVED_CHARS[section]):
            for k, v in record.ESCAPE_MAPPING:
                value = value.replace(k, v)
    # String-format the given value
    if value is None:
        return "."
    else:
        return str(value)


def format_value(field_info, value, section):
    """Format possibly compound value given the FieldInfo"""
    if section == "FORMAT" and field_info.id == "FT":
        if not value:
            return "."
        elif isinstance(value, list):
            return ";".join(map(lambda x: format_atomic(x, section), value))
    elif field_info.number == 1:
        if value is None:
            return "."
        else:
            return format_atomic(value, section)
    else:
        if not value:
            return "."
        else:
            return ",".join(map(lambda x: format_atomic(x, section), value))


class Writer:
    """Class for writing VCF files to ``file``-like objects

    Instead of using the constructor, use the class methods
    :py:meth:`~Writer.from_stream` and
    :py:meth:`~Writer.from_path`.

    The writer has to be constructed with a :py:class:`~vcfpy.header.Header`
    object and the full VCF header will be written immediately on construction.
    This, of course, implies that modifying the header after construction is
    illegal.
    """

    @classmethod
    def from_stream(klass, stream, header, path=None, use_bgzf=None):
        """Create new :py:class:`Writer` from file

        Note that for getting bgzf support, you have to pass in a stream
        opened in binary mode.  Further, you either have to provide a ``path``
        ending in ``".gz"`` or set ``use_bgzf=True``.  Otherwise, you will
        get the notorious "TypeError: 'str' does not support the buffer
        interface".

        :param stream: ``file``-like object to write to
        :param header: VCF header to use, lines and samples are deep-copied
        :param path: optional string with path to store (for display only)
        :param use_bgzf: indicator whether to write bgzf to ``stream``
            if ``True``, prevent if ``False``, interpret ``path`` if ``None``
        """
        if use_bgzf or (use_bgzf is None and path and path.endswith(".gz")):
            stream = bgzf.BgzfWriter(fileobj=stream)
        return Writer(stream, header, path)

    @classmethod
    def from_path(klass, path, header):
        """Create new :py:class:`Writer` from path

        :param path: the path to load from (converted to ``str`` for
            compatibility with ``path.py``)
        :param header: VCF header to use, lines and samples are deep-copied
        """
        path = str(path)
        use_bgzf = False  # we already interpret path
        if path.endswith(".gz"):
            f = bgzf.BgzfWriter(filename=path)
        else:
            f = open(path, "wt")
        return klass.from_stream(f, header, path, use_bgzf=use_bgzf)

    def __init__(self, stream, header, path=None):
        #: stream (``file``-like object) to read from
        self.stream = stream
        #: the :py:class:~vcfpy.header.Header` to write out, will be
        #: deep-copied into the ``Writer`` on initialization
        self.header = header.copy()
        #: optional ``str`` with the path to the stream
        self.path = path
        # write out headers
        self._write_header()

    def _write_header(self):
        """Write out the header"""
        for line in self.header.lines:
            print(line.serialize(), file=self.stream)
        if self.header.samples.names:
            print(
                "\t".join(list(parser.REQUIRE_SAMPLE_HEADER) + self.header.samples.names),
                file=self.stream,
            )
        else:
            print("\t".join(parser.REQUIRE_NO_SAMPLE_HEADER), file=self.stream)

    def close(self):
        """Close underlying stream"""
        self.stream.close()

    def write_record(self, record):
        """Write out the given :py:class:`vcfpy.record.Record` to this
        Writer"""
        self._serialize_record(record)

    def _serialize_record(self, record):
        """Serialize whole Record"""
        f = self._empty_to_dot
        row = [record.CHROM, record.POS]
        row.append(f(";".join(record.ID)))
        row.append(f(record.REF))
        if not record.ALT:
            row.append(".")
        else:
            row.append(",".join([f(a.serialize()) for a in record.ALT]))
        row.append(f(record.QUAL))
        row.append(f(";".join(record.FILTER)))
        row.append(f(self._serialize_info(record)))
        if record.FORMAT:
            row.append(":".join(record.FORMAT))
        row += [
            self._serialize_call(record.FORMAT, record.call_for_sample[s])
            for s in self.header.samples.names
        ]
        print(*row, sep="\t", file=self.stream)

    def _serialize_info(self, record):
        """Return serialized version of record.INFO"""
        result = []
        for key, value in record.INFO.items():
            info = self.header.get_info_field_info(key)
            if info.type == "Flag":
                result.append(key)
            else:
                result.append("{}={}".format(key, format_value(info, value, "INFO")))
        return ";".join(result)

    def _serialize_call(self, format_, call):
        """Return serialized version of the Call using the record's FORMAT'"""
        if isinstance(call, record.UnparsedCall):
            return call.unparsed_data
        else:
            result = [
                format_value(self.header.get_format_field_info(key), call.data.get(key), "FORMAT")
                for key in format_
            ]
            return ":".join(result)

    @classmethod
    def _empty_to_dot(klass, val):
        """Return val or '.' if empty value"""
        if val == "" or val is None or val == []:
            return "."
        else:
            return val

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.close()
