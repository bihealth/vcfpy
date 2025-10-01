import gzip
import logging
import struct
from io import RawIOBase
from typing import Dict, Generator, Iterable, List, Tuple, Union

from puretabix.fsm import FSMachine
from typing_extensions import Self

from .bgzip import BlockGZipReader
from .vcf import LINE_START, VCFAccumulator, VCFLine, get_vcf_fsm

logger = logging.getLogger(__name__)


class TabixIndex:
    def __init__(
        self,
        file_format: int,
        column_sequence: int,
        column_begin: int,
        column_end: int,
        meta: str,
        headerlines_count: int,
        indexes: Dict[str, Tuple[Dict[int, Tuple[Tuple[int, int], ...]], Tuple[int, ...]]],
    ):
        """
        In-memory representation of a Tabix index. See https://samtools.github.io/hts-specs/tabix.pdf
        for more information.

        Generally these are pretty small files that need to be read entirely, thus downloading them
        locally before processing is recommented e.g. io.BytesIO
        """

        # 0 = generic tab-delemited
        # 1 = SAM
        # 2 = VCF
        self.file_format = file_format
        # column for sequence ids, 1-based
        self.column_sequence = column_sequence
        # column for region start, 1-based
        self.column_begin = column_begin
        # column for region end, 1-based
        self.column_end = column_end
        self.meta = meta
        self.headerlines_count = headerlines_count

        # a dictionary of names to (bin_index, interval_index)
        self.indexes = indexes

    @classmethod
    def from_file(cls, fileobj: RawIOBase) -> Self:
        """
        Generally these are pretty small files that need to be read entirely, thus downloading them
        locally before processing is recommented e.g. io.BytesIO
        """
        # the index file is block-gzipped but small enough we can
        # load it into memory and process like a regular gzip file
        with gzip.GzipFile(fileobj=fileobj) as f:
            header_pattern = "<4siiiii4sii"
            header = struct.unpack(header_pattern, f.read(struct.calcsize(header_pattern)))

            magic = header[0]
            if magic != b"TBI\01":  # check magic
                raise RuntimeError(f"invalid tabix index magic {magic}.")

            # number of named sequences (e.g. chromosomes)
            n_sequences = header[1]

            file_format = header[2]
            # 0 = generic tab-delemited
            # 1 = SAM
            # 2 = VCF
            if file_format not in (0, 1, 2):
                raise RuntimeError(f"invalid tabix index format {file_format}.")

            # these are 1 based
            # value of 0 states not included in file
            # e.g. VCF has no explicit end column
            column_sequence = header[3]  # Column for the sequence name
            column_begin = header[4]  # Column for the start of a region
            column_end = header[5]  # Column for the end of a region

            # this is the comment marker, usually #
            meta = header[6].decode("ascii")[0]
            assert meta == "#", (header[6], meta)

            # number of lines of header at the start of the file
            # this does not include lines marked as comments
            headerlines_count = header[7]

            # sequence names are a series of bytes followed by a null byte
            names = tuple(map(bytes.decode, f.read(header[8]).split(b"\x00")[:-1]))  # throw the last empty one away
            if len(names) != n_sequences:
                raise RuntimeError(f"unexpected number of sequences {n_sequences} vs {len(names)}")

            indexes: Dict[str, Tuple[Dict[int, Tuple[Tuple[int, int], ...]], Tuple[int, ...]]] = {}
            # for each sequence
            for name in names:
                # each sequence has a bin index and an interval index

                # parse the bin index
                n_bins = struct.unpack("<i", f.read(4))[0]
                bins: Dict[int, Tuple[Tuple[int, int], ...]] = {}
                for _ in range(n_bins):
                    # each bin has a key, and a series of chunks
                    bin_key, n_chunks = struct.unpack("<Ii", f.read(8))
                    chunks: Tuple[Tuple[int, int], ...] = tuple(
                        (i, j) for i, j in struct.iter_unpack("<QQ", f.read(16 * n_chunks))
                    )

                    assert bin_key not in bins
                    bins[bin_key] = chunks

                # parse the interval index
                n_intervals = struct.unpack("<i", f.read(4))[0]
                intervals: Tuple[int, ...] = tuple(i[0] for i in struct.iter_unpack("<Q", f.read(8 * n_intervals)))

                if name in indexes:
                    raise RuntimeError(f"duplicate sequence name {name}")
                indexes[name] = (bins, intervals)

        return cls(
            file_format,
            column_sequence,
            column_begin,
            column_end,
            meta,
            headerlines_count,
            indexes,
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.file_format}, {self.column_sequence}, {self.column_begin}, {self.column_end}, {self.meta}, {self.headerlines_count}, {self.indexes})"

    def _lookup_linear(self, sequence_name: str, start: int) -> Union[None, int]:
        """
        For each tiling 16 kb window keep the virtual file offset of the leftmost record (i.e.
        having the smallest start coordinate) that overlaps the window.

        Given a location, get the smallest start location of all records that overlap the 16kb window
        containing the location
        """
        # linear index is in 16kb intervals
        # 16kb = 16 * (2 ** 10) = 2 ** 14 = 1 << 14
        # throw away the first 14 bits to get index position
        i = start >> 14
        # if this sequence_name isn't valid, say that
        if sequence_name not in self.indexes:
            return None
        linear_index = self.indexes[sequence_name][1]
        # if it would be beyond the index, say that
        if i >= len(linear_index):
            return None
        # its a valid sequnce name and a valid interval window
        return linear_index[i]

    def _lookup_bin_chunks(self, sequence_name: str, start: int, end: int) -> Generator[Tuple[int, int], None, None]:
        """
        Records are assigned to a bin if the entirely fit in the bin.
        So we want all the records in all the bins that overlap with the region of interest.
        These records *might* overlap with the region of interest.
        """
        bin_index = self.indexes[sequence_name][0]
        for chunks_bin_index in reversed(tuple(self.region_to_bins(start, end))):
            if chunks_bin_index in bin_index:
                for chunk in bin_index[chunks_bin_index]:
                    yield chunk

    def lookup_virtual(self, sequence_name: str, start: int, end: int) -> Union[Tuple[None, None], Tuple[int, int]]:
        virtual_start = None
        virtual_end = None

        linear_start = self._lookup_linear(sequence_name, start)
        # if this is not in the linear index, cant return anything
        if not linear_start:
            return None, None

        for chunk_start, chunk_end in self._lookup_bin_chunks(sequence_name, start, end):
            if chunk_end <= linear_start:
                # if the chunk finished before this section of the linear starts, skip the chunk
                # rare, but does happen sometimes
                continue

            # move the chunk start to where the linear start begins
            chunk_start = min(chunk_start, linear_start)

            if virtual_start is None or chunk_start < virtual_start:
                virtual_start = chunk_start

            if virtual_end is None or chunk_end > virtual_end:
                virtual_end = chunk_end

        # either both or neither must be set
        assert (virtual_start is None) == (virtual_end is None)
        assert (virtual_start is not None) == (virtual_end is not None)

        return virtual_start, virtual_end  # type: ignore

    def write_to(self, outfile: RawIOBase) -> None:
        # header
        outfile.write(
            struct.pack(
                "<4siiiii4si",
                b"TBI\01",  # magic number
                len(self.indexes),  # n sequences
                self.file_format,  # file format 0 generic, 1 sam, 2 vcf
                self.column_sequence,  # column for sequence ids, 1-based
                self.column_begin,  # column for region start, 1-based
                self.column_end,  # column for region end, 1-based
                self.meta.encode("ascii") + b"\x00\x00\x00",  # this is a character, but represented as a int
                self.headerlines_count,
            )
        )

        # length of concatenated zero terminated names
        names = tuple(self.indexes.keys())  # ensure consistent order
        names_concat = b"".join((i.encode("ascii") + b"\0" for i in names))
        outfile.write(struct.pack(f"<i{len(names_concat)}s", len(names_concat), names_concat))

        for name in names:
            # n_bin
            #   bin
            #   n_chunk
            #     chunk_begin
            #     chunk_end
            # n_intv
            #   ioff
            bin_index = self.indexes[name][0]
            outfile.write(struct.pack("<i", len(bin_index)))
            for bin_i in bin_index.keys():
                # unsigned bin number, n_chunk
                outfile.write(struct.pack("<Ii", bin_i, len(bin_index[bin_i])))
                for chunk_begin, chunk_end in bin_index[bin_i]:
                    outfile.write(struct.pack("<QQ", chunk_begin, chunk_end))

            intv_index = self.indexes[name][1]
            outfile.write(struct.pack("<i", len(intv_index)))
            for ioff in intv_index:
                outfile.write(struct.pack("<Q", ioff))

    @classmethod
    def build_from(cls, rawfile: RawIOBase) -> Self:
        """
        read a vcf file in blockgzip format and create an index object for it
        """

        file_format = 2  # 0 generic, 1 sam, 2 vcf
        column_sequence = 1  # column for sequence ids, 1-based
        column_begin = 2  # column for region start, 1-based
        column_end = 0  # column for region end, 1-based
        meta = "#"

        # dictionary of names to (bin_index, interval_index)
        # bin_index is a dictionary of bin numbers to [(chunk_start, chunk_end),...]
        # interval_index is a list of 16kbase interval start locations
        indexes: Dict[str, Tuple[Dict[int, List[Tuple[int, int]]], List[int]]] = {}

        # go through each line in turn
        # parse it into a structured object
        vcf_fsm = get_vcf_fsm()
        accumulator = VCFAccumulator()

        # these are the internal two index types
        bin_index = {}
        interval_index: List[int] = []

        bgzipped = BlockGZipReader(rawfile)
        bgzipped.seek(0)

        for (
            start_block,
            start_offset,
            end_block,
            end_offset,
            line,
        ) in bgzipped.generate_lines_offset():
            line_str = line.decode() + "\n"
            vcf_fsm.run(line_str, LINE_START, accumulator)
            vcf_line = accumulator.to_vcfline()
            accumulator.reset()
            # TODO add better debugging for unexpected lines

            # skip any comment lines with # or ##
            if vcf_line.comment_raw or vcf_line.comment_key:
                continue

            # have we started a new chromosome?
            if vcf_line.chrom not in indexes:
                chrom_bin_index: Dict[int, List[Tuple[int, int]]] = {}
                chrom_interval_index: List[int] = []
                indexes[vcf_line.chrom] = (chrom_bin_index, chrom_interval_index)
                bin_index = indexes[vcf_line.chrom][0]
                interval_index = indexes[vcf_line.chrom][1]

            # get the combined number for the block & offset
            start_virtual = start_block << 16 | start_offset
            end_virtual = end_block << 16 | end_offset

            # subtract 1 because its 0 offset
            record_start = vcf_line.pos - 1
            # subtract another 1 because half-open end
            record_end = vcf_line.pos - 1 + len(vcf_line.ref) - 1

            # bin index
            # smallest bin that completely contains the record

            bin_i = cls.region_to_bin(record_start, record_end)

            if bin_i not in bin_index.keys():
                bin_index[bin_i] = [(start_virtual, end_virtual + 1)]
            else:
                # extend chunk if directly continuous
                if bin_index[bin_i][-1][1] == start_virtual:
                    bin_index[bin_i][-1] = (
                        bin_index[bin_i][-1][0],
                        end_virtual + 1,
                    )
                else:
                    bin_index[bin_i].append((start_virtual, end_virtual + 1))

            # interval index
            # is the lowest virtual offset of all records that overlap interval
            # half-closed half-open records

            # throw away the first 14 bits to get interval index
            # e.g. 0 => 0, 16384 => 1, etc
            interval_i = record_start >> 14

            # line fully in block
            if start_block == end_block:
                # pad if necessary
                while len(interval_index) <= interval_i:
                    interval_index.append(start_virtual)

        # now they have been built, freeze into immutability
        # turn Dict[str, Tuple[Dict[int, List[ Tuple[int, int]]]     , List[int]]]
        # into Dict[str, Tuple[Dict[int, Tuple[Tuple[int, int], ...]], Tuple[int, ...]]]
        indexes_frozen: Dict[str, Tuple[Dict[int, Tuple[Tuple[int, int], ...]], Tuple[int, ...]]] = {}
        for chrom in indexes:
            bin_index_frozen: Dict[int, Tuple[Tuple[int, int], ...]] = {}
            for i in indexes[chrom][0]:
                bin_index_frozen[i] = tuple(indexes[chrom][0][i])
            interval_index_frozen = tuple(indexes[chrom][1])
            indexes_frozen[chrom] = (bin_index_frozen, interval_index_frozen)

        return cls(
            file_format,
            column_sequence,
            column_begin,
            column_end,
            meta,
            0,
            indexes_frozen,
        )

    @staticmethod
    def region_to_bins(begin: int, end: int, n_levels: int = 5, min_shift: int = 14) -> Generator[int, None, None]:
        """
        generator of keys to bins of records which *may* overlap the given region

        n_levels: int, optional
            cluster level, 5 for tabix
        min_shift: int, optional
            minimum shift, 14 for tabix
        """
        t = 0
        s = min_shift + (n_levels << 1) + n_levels
        for level in range(n_levels + 1):
            b = t + (begin >> s)
            e = t + (end >> s)
            n = e - b + 1
            for k in range(b, e + 1):
                yield k
                n += 1
            t += 1 << ((level << 1) + level)
            s -= 3

    @staticmethod
    def region_to_bin(begin: int, end: int) -> int:
        """
        returns the index of the smallest bin that contains the region
        as a half-closed half-open interval
        """
        if begin >> 14 == end >> 14:
            return ((1 << 15) - 1) // 7 + (begin >> 14)
        elif begin >> 17 == end >> 17:
            return ((1 << 12) - 1) // 7 + (begin >> 17)
        elif begin >> 20 == end >> 20:
            return ((1 << 9) - 1) // 7 + (begin >> 20)
        elif begin >> 23 == end >> 23:
            return ((1 << 6) - 1) // 7 + (begin >> 23)
        elif begin >> 26 == end >> 26:
            return ((1 << 3) - 1) // 7 + (begin >> 26)
        return 0

    @staticmethod
    def bin_start(k: int) -> int:
        # bit_length-1 is int log2
        lvl = (((7 * k) + 1).bit_length() - 1) // 3
        ol: int = ((2 ** (3 * lvl)) - 1) // 7
        sl: int = 2 ** (29 - (3 * lvl))
        start = (k - ol) * sl
        return start

    @staticmethod
    def bin_size(k: int) -> int:
        # bit_length-1 is int log2
        lvl = (((7 * k) + 1).bit_length() - 1) // 3
        sl: int = 2 ** (29 - (3 * lvl))
        return sl


class TabixIndexedFile:
    def __init__(self, fileobj: RawIOBase, index: TabixIndex):
        self.index = index
        self.bgzipped = BlockGZipReader(fileobj)

    @classmethod
    def from_files(cls, fileobj: RawIOBase, index_fileobj: RawIOBase) -> Self:
        return cls(fileobj, TabixIndex.from_file(index_fileobj))

    def fetch_bytes_block_offset(self, block_start: int, offset_start: int, block_end: int, offset_end: int) -> bytes:
        value = b""
        self.bgzipped.seek(block_start)
        block = block_start
        while block <= block_end:
            _, _, decompressed, _ = self.bgzipped.get_block()
            # print(decompressed)
            # empty block at end of file
            if not decompressed:
                break

            if block == block_end:
                # end block, drop beyond offset end
                decompressed = decompressed[:offset_end]
            if block == block_start:
                # start block, drop before offset start
                decompressed = decompressed[offset_start:]
            value = value + decompressed
            # update ready for next loop
            block = self.bgzipped.tell()
        return value

    def fetch_bytes_virtual(self, virtual_start: int, virtual_end: int) -> bytes:
        # the lower 16 bits store the offset of the byte inside the gzip block
        # the rest store the offset of gzip block
        block_start = virtual_start >> 16
        offset_start = virtual_start & 0xFFFF
        block_end = virtual_end >> 16
        offset_end = virtual_end & 0xFFFF
        return self.fetch_bytes_block_offset(block_start, offset_start, block_end, offset_end)

    def fetch_bytes(self, name: str, start: int, end: Union[None, int]) -> bytes:
        """
        Returns bytes in the region of interest
        """
        # quick check
        if name not in self.index.indexes.keys():
            return b""

        # default if only start specified
        if not end:
            end = start

        # use the index to get "virtual" file offsets that include all of the region of interest
        virtual_start, virtual_end = self.index.lookup_virtual(name, start, end)

        # location not indexed, return empty string
        if not virtual_start and not virtual_end:
            return b""

        # check for weirdness if only one is missing
        if virtual_start is None or virtual_end is None:
            raise RuntimeError()

        return self.fetch_bytes_virtual(virtual_start, virtual_end - 1)

    def fetch_lines(self, name: str, start: int, end: Union[None, int]) -> Generator[str, None, None]:
        """
        Returns lines in the region of interest
        """
        # default if only start specified
        if not end:
            end = start

        region = self.fetch_bytes(name, start, end)
        # no region, no lines
        if not region:
            return ("" for _ in [])

        # turn the bytes into a list of strings
        lines: Iterable[str] = region.decode("utf-8").splitlines()

        # filter out comments
        lines = (line for line in lines if not line.startswith(self.index.meta))

        # split lines
        lines_split = ((line, line.split("\t")) for line in lines)

        # filter lines of wrong lengths i.e. cut off around chunk boundries
        expected_len = max(
            self.index.column_sequence,
            self.index.column_begin,
            self.index.column_end,
        )

        # filter lines before start and after end
        column_begin = self.index.column_begin
        # default to using begin column again
        column_end = self.index.column_end if self.index.column_end else self.index.column_begin

        return (
            line
            for line, line_split in lines_split
            if (
                len(line_split) >= expected_len  # right length
                and int(line_split[column_begin - 1]) >= start  # after start
                and int(line_split[column_end - 1]) <= end  # before end
            )
        )

    def fetch(self, name: str, start: int, end: Union[None, int] = None) -> str:
        """
        Returns region of interest
        """
        return "\n".join(self.fetch_lines(name, start, end))


class TabixIndexedVCFFile(TabixIndexedFile):
    accumulator: VCFAccumulator
    vcf_fsm: FSMachine

    def __init__(self, fileobj: RawIOBase, index: TabixIndex):
        super().__init__(fileobj, index)
        self.vcf_fsm = get_vcf_fsm()
        self.accumulator = VCFAccumulator()

    def fetch_vcf_lines(self, name: str, start: int, end: Union[None, int] = None) -> Generator[VCFLine, None, None]:
        # default if only start specified
        if not end:
            end = start

        for line in self.fetch_bytes(name, start, end).decode("utf-8").splitlines():
            try:
                self.vcf_fsm.run(line, LINE_START, self.accumulator)
            except ValueError as e:
                logger.error(f"Error parsing {line}")
                raise e
            vcfline = self.accumulator.to_vcfline()
            self.accumulator.reset()
            if vcfline.pos >= start and vcfline.pos <= end:  # after start and before end
                yield vcfline


from typing import Dict, Generator, Mapping, Optional, Tuple

from typing_extensions import Self

from .fsm import (
    FSMachine,
    RegexTransition,
    SetInTransition,
    SetNotInTransition,
)


class VCFLine:
    """
    Representation of a single VCF line.

    Comment lines (starting #) will have:
      comment_raw

    Meta-information (starting ##) will have
      comment_key
      comment_value_str
      comment_value_dict

    Data lines will have:
        chrom       VCF CHROM
        pos         VCF POS
        _id         VCF IDs (if any)
        ref         VCF REF
        alt         VCF ALTs (if any)
        qual        floating point representation of VCF QUAL (if possible)
        qual_str    raw string representation of VCF QUAL
        _filter     VCF FILTERs (if any)
        info        VCF INFO as a dictionary

    Data lines may have:
        sample      Iterable of dictionaries for the samples

    For output, a VCFLine can be converted to a string.

    For convenience, several additional features exist

        is_comment  checks if either comment_raw or comment_key is defined
        get_gt()    returns the ref and/or alt alleles for each sample as a tuple of tuples of strings
    """

    comment_raw: str  # non-empty if line starts with #
    comment_key: str  # non-empty if line starts with ##
    comment_value_str: str  # non-empty if line starts with ##
    comment_value_dict: Dict[str, Optional[str]]
    chrom: str
    pos: int
    _id: Tuple[str, ...]
    ref: str
    alt: Tuple[str, ...]
    qual: Optional[float]
    qual_str: str
    _filter: Tuple[str, ...]
    info: Dict[str, Optional[Iterable[str]]]
    sample: Tuple[Dict[str, str], ...]

    __slots__ = (
        "comment_raw",
        "comment_key",
        "comment_value_str",
        "comment_value_dict",
        "chrom",
        "pos",
        "_id",
        "ref",
        "alt",
        "qual",
        "qual_str",
        "_filter",
        "info",
        "sample",
    )

    def __init__(
        self,
        comment_raw: str,
        comment_key: str,
        comment_value_str: str,
        comment_value_dict: Dict[str, Optional[str]],
        chrom: str,
        pos: int,
        _id: Iterable[str],
        ref: str,
        alt: Iterable[str],
        qual_str: str,
        _filter: Iterable[str],
        info: Mapping[str, Optional[Iterable[str]]],
        sample: Iterable[Mapping[str, str]],
    ):
        # TODO validation of permitted characters
        self.comment_raw = comment_raw
        self.comment_key = comment_key
        self.comment_value_str = comment_value_str
        self.comment_value_dict = comment_value_dict
        self.chrom = chrom
        self.pos = pos
        self._id = tuple(_id)
        self.ref = ref
        self.alt = tuple(alt)
        self.qual_str = qual_str
        self._filter = tuple(_filter)
        self.info = dict(info)
        # this may be zero to many
        self.sample = tuple(dict(s) for s in sample)

        try:
            self.qual = float(qual_str)
        except ValueError:
            # if we can't do the conversion, don't worry
            # value will be None but original in qual_str
            self.qual = None

    def __str__(self) -> str:
        # CHROM POS ID REF ALT QUAL FILTER INFO (FORMAT) (SAMPLE) (SAMPLE) ...
        if self.comment_raw:
            return "#" + self.comment_raw
        elif self.comment_key and self.comment_value_str:
            # meta information line
            # ##fileformat=VCFv4.2
            return "".join(("##", self.comment_key, "=", self.comment_value_str))
        elif self.comment_key and self.comment_value_dict:
            # meta information line
            # ##FILTER=<ID=ID,Description="description">
            return "".join(
                (
                    "##",
                    self.comment_key,
                    "=<",
                    ",".join((f"{key}={value}" if value else key for key, value in self.comment_value_dict.items())),
                    ">",
                )
            )
        else:
            # data line
            required = "\t".join(
                (
                    self.chrom,
                    str(self.pos),
                    ";".join(self._id),
                    self.ref,
                    ",".join(self.alt),
                    self.qual_str,  # use non-lossy stored version
                    ";".join(self._filter),
                    ";".join(
                        (f"{key}={','.join(value)}" if value else key for key, value in self.info.items())
                        if self.info
                        else "."
                    ),
                )
            )

            if not self.sample:
                return required
            else:
                # get the superset of all keys of all samples
                keyset = set()
                keylist = []
                for sample in self.sample:
                    for key in sample.keys():
                        if key not in keyset:
                            keylist.append(key)
                            keyset.add(key)
                sample_keys = ":".join(keylist)
                line = "\t".join((required, sample_keys))

                # get the values for each key in superset or .
                for sample in self.sample:
                    values = []
                    for key in keylist:
                        values.append(sample.get(key, "."))
                    line = line + "\t" + ":".join(values)
                return line

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({repr(self.comment_raw)},{repr(self.comment_key)},{repr(self.comment_value_str)},"
            + f"{repr(self.comment_value_dict)},{repr(self.chrom)},{repr(self.pos)},{repr(self._id)},{repr(self.ref)},{repr(self.alt)},"
            + f"{repr(self.qual_str)},{repr(self._filter)},{repr(self.info)},{repr(self.sample)})"
        )

    @property
    def is_comment(self) -> bool:
        return bool(self.comment_raw or self.comment_key)

    def get_genotype(self) -> Tuple[Tuple[str, ...], ...]:
        """
        Convenience function to get the GT of each sample, and use that to match the allels in ref & alt.

        Returns a tuple of strings.
        """
        if self.is_comment:
            raise ValueError("Cannot get_gt() of a comment")

        result: List[Tuple[str, ...]] = [() for _ in range(len(self.sample))]
        refalt = (self.ref,) + self.alt
        for i, sample in enumerate(self.sample):
            if "GT" in sample:
                # split into separate alleles
                allele_pos = sample["GT"].replace("|", "/").split("/")
                # for each allele either keep a dot or do ref+alt lookup
                alleles = [i if i == "." else refalt[int(i)] for i in allele_pos]
                result[i] = tuple(alleles)
        return tuple(result)

    @classmethod
    def as_comment_raw(cls, comment_raw: str) -> "VCFLine":
        return cls(
            comment_raw,
            comment_key="",
            comment_value_str="",
            comment_value_dict={},
            chrom="",
            pos=0,
            _id=[],
            ref="",
            alt=[],
            qual_str="",
            _filter="",
            info={},
            sample=[],
        )

    @classmethod
    def as_comment_key_string(cls, key: str, value_str: str) -> "VCFLine":
        # for example ##fileformat=VCFv4.3
        return cls(
            comment_raw="",
            comment_key=key,
            comment_value_str=value_str,
            comment_value_dict={},
            chrom="",
            pos=0,
            _id=[],
            ref="",
            alt=[],
            qual_str="",
            _filter="",
            info={},
            sample=[],
        )

    @classmethod
    def as_comment_key_dict(cls, key: str, value_dict: Mapping[str, str]) -> "VCFLine":
        # for example
        # ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
        return cls(
            comment_raw="",
            comment_key=key,
            comment_value_str="",
            comment_value_dict=dict(value_dict),
            chrom="",
            pos=0,
            _id=[],
            ref="",
            alt=[],
            qual_str="",
            _filter="",
            info={},
            sample=[],
        )

    @classmethod
    def as_data(
        cls,
        chrom: str,
        pos: int,
        _id: Iterable[str],
        ref: str,
        alt: Iterable[str],
        qual_str: str,
        _filter: Iterable[str],
        info: Mapping[str, Optional[Iterable[str]]],
        sample: Iterable[Mapping[str, str]],
    ) -> Self:
        return cls(
            comment_raw="",
            comment_key="",
            comment_value_str="",
            comment_value_dict={},
            chrom=chrom,
            pos=pos,
            _id=_id,
            ref=ref,
            alt=alt,
            qual_str=qual_str,
            _filter=_filter,
            info=info,
            sample=sample,
        )


class VCFAccumulator:
    __characters: List[str]
    comment_struct: Dict[str, Optional[str]]
    _id: List[str]
    alt: List[str]
    _filter: List[str]
    info: Dict[str, Optional[List[str]]]
    format: List[str]
    samples: List[Dict[str, str]]

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        self.__characters = []

        self.comment_raw = ""
        self.comment_key = ""
        self.comment_value = ""
        self.comment_struct = {}
        self._comment_struct_key = ""
        self._comment_struct_value = ""
        self.chrom = ""
        self.pos = 0
        self._id = []
        self.ref = ""
        self.alt = []
        self._alt_option = ""
        self.qual = ""
        self._filter = []
        self.info = {}
        self.__infokey = ""
        self.format = []
        self.samples = []

    def to_vcfline(self) -> VCFLine:
        return VCFLine(
            self.comment_raw,
            self.comment_key,
            self.comment_value,
            self.comment_struct,
            self.chrom,
            self.pos,
            self._id,
            self.ref,
            self.alt,
            self.qual,
            self._filter,
            self.info,
            self.samples,
        )

    def append_character(self, _char: str) -> None:
        self.__characters.append(_char)

    def end_comment(self, _char: str) -> None:
        self.comment_raw = "".join(self.__characters)
        self.__characters = []

    def comment_value_to_end(self, _char: str) -> None:
        self.comment_value = "".join(self.__characters)
        self.__characters = []

    def comment_key_to_comment_value(self, _char: str) -> None:
        self.comment_key = "".join(self.__characters)
        self.__characters = []

    def comment_value_to_comment_struct_key(self, _char: str) -> None:
        self.comment_value = "".join(self.__characters)
        self.__characters = []

    def comment_struct_key_to_comment_struct_value(self, _char: str) -> None:
        comment_struct_key = "".join(self.__characters)
        assert self._comment_struct_key not in self.comment_struct, (
            self._comment_struct_key,
            self.comment_struct,
        )
        self.comment_struct[comment_struct_key] = None
        self.__characters = []

    def comment_struct_value_to_comment_struct_key(self, _char: str) -> None:
        self._comment_struct_value = "".join(self.__characters)
        self.__characters = []
        # this needs dict keys to be insertion ordered
        self.comment_struct[tuple(self.comment_struct.keys())[-1]] = self._comment_struct_value

    def chrom_to_pos(self, _char: str) -> None:
        self.chrom = "".join(self.__characters)
        self.__characters = []

    def pos_to_id(self, _char: str) -> None:
        self.pos = int("".join(self.__characters))
        self.__characters = []

    def id_to_ref(self, _char: str) -> None:
        self._id.append("".join(self.__characters))
        self.__characters = []

    def ref_to_alt(self, _char: str) -> None:
        self.ref = "".join(self.__characters)
        self.__characters = []

    def alt_to_alt(self, _char: str) -> None:
        self.alt.append("".join(self.__characters))
        self.__characters = []

    def alt_to_qual(self, _char: str) -> None:
        self.alt.append("".join(self.__characters))
        self.__characters = []

    def qual_to_filter(self, _char: str) -> None:
        self.qual = "".join(self.__characters)
        self.__characters = []

    def filter_to_info_key(self, _char: str) -> None:
        self._filter.append("".join(self.__characters))
        self.__characters = []

    def info_key_to_info_key(self, _char: str) -> None:
        self.__infokey = "".join(self.__characters)
        self.info[self.__infokey] = None
        self.__characters = []

    def info_key_to_info_value(self, _char: str) -> None:
        self.__infokey = "".join(self.__characters)
        self.__characters = []

    def info_value_to_format(self, _char: str) -> None:
        __infovalue = "".join(self.__characters)
        if self.__infokey in self.info and self.info[self.__infokey] is not None:
            # cannot get here if none, so this won't error
            self.info[self.__infokey].append(__infovalue)  # type: ignore[union-attr]
        else:
            self.info[self.__infokey] = [__infovalue]
        self.__characters = []

    def format_to_sample(self, _char: str) -> None:
        self.format.append("".join(self.__characters))
        self.__characters = []

    def sample_to_sample(self, _char: str) -> None:
        sample_str = "".join(self.__characters)
        sample_parts = sample_str.split(":")
        sample = {}
        for key, value in zip(self.format, sample_parts, strict=False):
            sample[key] = value
        self.samples.append(sample)
        self.__characters = []


# states for the finite state machine to parse comments
LINE_START = "LINE_START"
COMMENT = "COMMENT"
COMMENT_KEY = "COMMENT_KEY"
COMMENT_VALUE = "COMMENT_VALUE"
COMMENT_STRUCT_KEY = "COMMENT_STRUCT_KEY"
COMMENT_STRUCT_VALUE = "COMMENT_STRUCT_VALUE"
COMMENT_STRUCT_VALUE_QUOTED = "COMMENT_STRUCT_VALUE_QUOTED"

# ##CHROM POS ID REF ALT QUAL FILTER INFO (FORMAT) (SAMPLE)*
CHROM = "CHROM"
POS = "POS"
ID = "ID"
REF = "REF"
ALT = "ALT"
QUAL = "QUAL"
FILTER = "FILTER"
INFO_KEY = "INFO_KEY"
INFO_VALUE = "INFO_VALUE"
FORMAT = "FORMAT"
SAMPLE = "SAMPLE"


def get_vcf_fsm() -> FSMachine:
    fsm_vcf = FSMachine()

    fsm_vcf.add_transition(LINE_START, CHROM, SetNotInTransition, "#", VCFAccumulator.append_character)
    fsm_vcf.add_transition(CHROM, CHROM, SetNotInTransition, "\t", VCFAccumulator.append_character)
    # TODO verify CHROM in contig or assembly from comments?
    fsm_vcf.add_transition(CHROM, POS, SetInTransition, "\t", VCFAccumulator.chrom_to_pos)
    fsm_vcf.add_transition(POS, POS, SetInTransition, "0123456789", VCFAccumulator.append_character)
    fsm_vcf.add_transition(POS, ID, SetInTransition, "\t", VCFAccumulator.pos_to_id)
    fsm_vcf.add_transition(ID, ID, RegexTransition, r"\S", VCFAccumulator.append_character)
    fsm_vcf.add_transition(ID, ID, SetInTransition, ";", VCFAccumulator.id_to_ref)
    fsm_vcf.add_transition(ID, REF, SetInTransition, "\t", VCFAccumulator.id_to_ref)
    fsm_vcf.add_transition(REF, REF, SetInTransition, "ACGTN", VCFAccumulator.append_character)
    fsm_vcf.add_transition(REF, ALT, SetInTransition, "\t", VCFAccumulator.ref_to_alt)
    fsm_vcf.add_transition(ALT, ALT, SetNotInTransition, ",\t", VCFAccumulator.append_character)
    fsm_vcf.add_transition(ALT, ALT, SetInTransition, ",", VCFAccumulator.alt_to_alt)
    fsm_vcf.add_transition(ALT, QUAL, SetInTransition, "\t", VCFAccumulator.alt_to_qual)
    fsm_vcf.add_transition(QUAL, QUAL, SetInTransition, "0123456789.-", VCFAccumulator.append_character)
    fsm_vcf.add_transition(QUAL, FILTER, SetInTransition, "\t", VCFAccumulator.qual_to_filter)
    fsm_vcf.add_transition(FILTER, FILTER, SetNotInTransition, "\t;", VCFAccumulator.append_character)
    fsm_vcf.add_transition(FILTER, FILTER, SetInTransition, ";", VCFAccumulator.filter_to_info_key)
    fsm_vcf.add_transition(FILTER, INFO_KEY, SetInTransition, "\t", VCFAccumulator.filter_to_info_key)
    fsm_vcf.add_transition(
        INFO_KEY,
        INFO_KEY,
        SetNotInTransition,
        "=\t;\n",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        INFO_KEY,
        INFO_VALUE,
        SetInTransition,
        "=",
        VCFAccumulator.info_key_to_info_value,
    )
    fsm_vcf.add_transition(INFO_KEY, INFO_KEY, SetInTransition, ";", VCFAccumulator.info_key_to_info_key)
    fsm_vcf.add_transition(INFO_KEY, FORMAT, SetInTransition, "\t", VCFAccumulator.info_key_to_info_key)
    fsm_vcf.add_transition(
        INFO_KEY,
        None,
        SetInTransition,
        ("\n", None),
        VCFAccumulator.info_key_to_info_key,
    )
    fsm_vcf.add_transition(
        INFO_VALUE,
        INFO_VALUE,
        SetNotInTransition,
        "\t;,\n",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        INFO_VALUE,
        INFO_VALUE,
        SetInTransition,
        ",",
        VCFAccumulator.info_value_to_format,
    )
    fsm_vcf.add_transition(INFO_VALUE, INFO_KEY, SetInTransition, ";", VCFAccumulator.info_value_to_format)
    fsm_vcf.add_transition(INFO_VALUE, FORMAT, SetInTransition, "\t", VCFAccumulator.info_value_to_format)
    fsm_vcf.add_transition(
        INFO_VALUE,
        None,
        SetInTransition,
        ("\n", None),
        VCFAccumulator.info_value_to_format,
    )
    fsm_vcf.add_transition(FORMAT, FORMAT, SetNotInTransition, "\t:", VCFAccumulator.append_character)
    fsm_vcf.add_transition(FORMAT, FORMAT, SetInTransition, ":", VCFAccumulator.format_to_sample)
    fsm_vcf.add_transition(FORMAT, SAMPLE, SetInTransition, "\t", VCFAccumulator.format_to_sample)
    fsm_vcf.add_transition(
        SAMPLE,
        SAMPLE,
        SetNotInTransition,
        ("\t", "\n", None),
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(SAMPLE, SAMPLE, SetInTransition, "\t", VCFAccumulator.sample_to_sample)
    fsm_vcf.add_transition(SAMPLE, None, SetInTransition, ("\n", None), VCFAccumulator.sample_to_sample)
    fsm_vcf.add_transition(LINE_START, COMMENT, SetInTransition, "#", None)
    fsm_vcf.add_transition(
        COMMENT,
        COMMENT,
        SetNotInTransition,
        ("#", "\n", None),
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(COMMENT, None, SetInTransition, ("\n", None), VCFAccumulator.end_comment)
    fsm_vcf.add_transition(COMMENT, COMMENT_KEY, SetInTransition, "#", None)
    fsm_vcf.add_transition(
        COMMENT_KEY,
        COMMENT_KEY,
        SetNotInTransition,
        "=",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_KEY,
        COMMENT_VALUE,
        SetInTransition,
        "=",
        VCFAccumulator.comment_key_to_comment_value,
    )
    fsm_vcf.add_transition(
        COMMENT_VALUE,
        COMMENT_VALUE,
        SetNotInTransition,
        "<\n",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_VALUE,
        COMMENT_STRUCT_KEY,
        SetInTransition,
        "<",
        VCFAccumulator.comment_value_to_comment_struct_key,
    )
    fsm_vcf.add_transition(
        COMMENT_VALUE,
        None,
        SetInTransition,
        ("\n", None),
        VCFAccumulator.comment_value_to_end,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_KEY,
        COMMENT_STRUCT_KEY,
        SetNotInTransition,
        "=",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_KEY,
        COMMENT_STRUCT_VALUE,
        SetInTransition,
        "=",
        VCFAccumulator.comment_struct_key_to_comment_struct_value,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT_STRUCT_VALUE,
        SetNotInTransition,
        '",>',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT_STRUCT_VALUE_QUOTED,
        SetInTransition,
        '"',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE_QUOTED,
        COMMENT_STRUCT_VALUE_QUOTED,
        SetNotInTransition,
        '"',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE_QUOTED,
        COMMENT_STRUCT_VALUE,
        SetInTransition,
        '"',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT_STRUCT_KEY,
        SetInTransition,
        ",",
        VCFAccumulator.comment_struct_value_to_comment_struct_key,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT,
        SetInTransition,
        ">",
        VCFAccumulator.comment_struct_value_to_comment_struct_key,
    )

    return fsm_vcf


def read_vcf_lines(input_: Iterable[str], header_only: bool = False) -> Generator[VCFLine, None, None]:
    """
    Convenience function for parsing a source of VCF lines
    """
    vcf_fsm = get_vcf_fsm()
    accumulator = VCFAccumulator()
    for line in input_:
        if line:
            vcf_fsm.run(line, LINE_START, accumulator)
            vcfline = accumulator.to_vcfline()
            accumulator.reset()

            # if we only want to get the initial header
            # and this is not part of the header
            # then stop
            # Note: this does include the vcf column headings (sample names)
            if not vcfline.is_comment and header_only:
                break

            yield vcfline


"""
Finite State Machine
--------------------

Simple implementation of a finite state machine, used for parsing in a
character-by-character manner


based on https://gist.github.com/brianray/8d3e697dbbf150f725291d74ac0cee8b
"""

import logging
import re
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    Mapping,
    Optional,
    Pattern,
    Type,
    Union,
)

logger = logging.getLogger(__name__)

# note: Optional[Callable[[Any, str], None]] isn't the optimal type hint for the callbacks
# In theory, they can take an arbitrary number of parameters as long as the last position is a str
# and the others are matched via callback_args/callback_kwargs
# This might be able to be captured in newer Python versions with Protocol and/or ParamSpec


class RegexTransition:
    dst: Any
    condition: Pattern[str]
    callback: Optional[Callable[[Any, str], None]]
    __slots__ = ["dst", "condition", "callback"]

    def __init__(
        self,
        destination_state: Any,
        condition: str,
        callback: Callable[[Any, str], None],
    ):
        self.dst = destination_state
        self.condition = re.compile(condition)
        self.callback = callback

    def match(self, _input: str) -> bool:
        return bool(self.condition.match(_input))


class SetInTransition:
    dst: Any
    condition: FrozenSet[str]
    callback: Optional[Callable[[Any, str], None]]
    __slots__ = ["dst", "condition", "callback"]

    def __init__(
        self,
        destination_state: Any,
        condition: Iterable[str],
        callback: Optional[Callable[[Any, str], None]],
    ):
        self.dst = destination_state
        self.condition = frozenset(condition)
        self.callback = callback

    def match(self, _input: Union[str, bytes]) -> bool:
        return bool(_input in self.condition)


class SetNotInTransition:
    dst: Any
    condition: FrozenSet[str]
    callback: Optional[Callable[[Any, str], None]]

    __slots__ = ["dst", "condition", "callback"]

    def __init__(
        self,
        destination_state: Any,
        condition: Iterable[str],
        callback: Optional[Callable[[Any, str], None]],
    ):
        self.dst = destination_state
        self.condition = frozenset(condition)
        self.callback = callback

    def match(self, _input: Union[str, bytes]) -> bool:
        return bool(_input not in self.condition)


class FSMachine:
    transitions: Dict[Any, Any]

    def __init__(self) -> None:
        self.transitions = {}

    def add_transition(
        self,
        start_state: Any,
        end_state: Any,
        transition_class: Type[Any],
        condition: Any,
        callback: Optional[Callable[[Any, str], None]],
    ) -> None:
        """

        callback should be a function that accepts the current character of the state
        machine as the last positional arguments i.e. it is called like:
        callback(*args, current_char, **kwargs)

        because callback uses the last positional argument as the character, the first
        positional argment can be used by a class instance as self
        """
        if start_state not in self.transitions:
            self.transitions[start_state] = []
        self.transitions[start_state].append(transition_class(end_state, condition, callback))

    def run(
        self,
        inputs: Iterable[str],
        initial_state: Any,
        *args: Any,
        **kwargs: Dict[Any, Any],
    ) -> None:
        self.current_state = initial_state
        for c in inputs:
            self.process_next(c, args, kwargs)
            # if state is None, early exit
            if not self.current_state:
                break

        # process that we reached the end of the input
        if self.current_state:
            self.process_next(None, args, kwargs)

        # check at a valid end state
        assert not self.current_state, f"Unexpected ending at {self.current_state}"

    def process_next(
        self,
        _input: Optional[str],
        callback_args: Any,
        callback_kwargs: Mapping[Any, Any],
    ) -> bool:
        frozen_state = self.current_state
        for transition in self.transitions[frozen_state]:
            if transition.match(_input):
                # found a transition that matches
                # update the state
                self.current_state = transition.dst
                # call the callback, if it exists
                # because callback uses the last positional argument as the input, the first
                # positional argment can be used by a class instance as self
                if transition.callback:
                    transition.callback(*callback_args, _input, **callback_kwargs)
                # say that we matched this input
                return True
        raise ValueError(f"Unrecognized input {_input} in state {frozen_state}")


import io
import logging
import zlib
from typing import Generator, Optional, Tuple

from typing_extensions import Buffer

logger = logging.getLogger(__name__)

headerpattern = "<BBBBIBBHBBHH"
headersize = struct.calcsize(headerpattern)
tailpattern = "<II"
tailsize = struct.calcsize(tailpattern)


class BlockGZipReader:
    raw: io.IOBase

    def __init__(self, raw: io.IOBase):
        assert raw.seekable()
        self.raw = raw
        assert self.check_is_block_gzip()

    def seek(self, offset: int) -> int:
        return self.raw.seek(offset)

    def tell(self) -> int:
        return self.raw.tell()

    def check_is_gzip(self) -> bool:
        """
        returns a boolean for if it has a gzip header
        will seek to start of file
        """
        self.raw.seek(0)
        bytes_data = self.raw.read(3)
        header = struct.unpack("<BBB", bytes_data)
        return bool(header[0] == 31 and header[1] == 139 and header[2] == 8)

    def check_is_block_gzip(self) -> bool:
        """
        returns a boolean for if it has a block gzip header
        also checks if it has a gzip header
        will seek to start of file
        """
        if not self.check_is_gzip():
            return False
        # NOTE assumes there is only one extra header
        # not sure if this is required by block gzip spec or not
        self.raw.seek(12)
        bytes_data = self.raw.read(4)
        header = struct.unpack("<ccH", bytes_data)
        return bool(header[0] == b"B" and header[1] == b"C" and header[2] == 2)

    @staticmethod
    def check_is_header(header: Tuple[int, ...]) -> bool:
        """
        tests if a series of integers matches a blockgzip header
        """
        if header[0] != 31:
            return False
        if header[1] != 139:
            return False
        if header[2] != 8:
            return False
        if header[3] != 4:
            return False
        if header[8] != 66:
            return False
        if header[9] != 67:
            return False
        if header[10] != 2:
            return False
        return True

    def get_header(self) -> Tuple[int, ...]:
        """
        reads the next header from the file from the current point, assuming file is currently
        at start of a block
        """
        bytesread = self.raw.read(headersize)
        if len(bytesread) < headersize:
            raise EOFError(f"Expected to read {headersize} read {len(bytesread)}")
        header = struct.unpack(headerpattern, bytesread)
        assert self.check_is_header(header)
        return header

    def get_cdata(self, header: Tuple[int, ...]) -> bytes:
        """
        reads the compressed data of a block from the current point, given the bytes from the
        header of that block to determine size
        """
        blocksize = header[11] - header[7] - 19
        cdata = self.raw.read(blocksize)
        assert len(cdata) == blocksize, f"Unable to read up to {blocksize} of cdata"
        return bytes(cdata)

    def get_tail(self, decompressed: Optional[bytes] = None) -> Tuple[int, int]:
        """
        reads the tail of the block from the current point as a tuple of crc and isize
        if the decompressed bytes are provided, will use the crc checksum in the tail
        to validate the bytes match expectation
        """
        # read isize and crc check
        tailbytes = self.raw.read(tailsize)
        if len(tailbytes) != tailsize:
            raise ValueError(f"Unable to read {tailsize} bytes for tail")
        tail_crc, tail_isize = struct.unpack(tailpattern, tailbytes)
        # if we were given the decompressed data, check it matches expectation
        if decompressed:
            # check decompressed size is expected
            assert len(decompressed) == tail_isize
            # check crc check is expected
            assert zlib.crc32(decompressed) == tail_crc
        return tail_crc, tail_isize

    def get_cdata_decompressed(self, header: Tuple[int, ...]) -> Tuple[bytes, bytes]:
        """
        reads the compressed data of a block from the current point, given the bytes from the
        header of that block to determine size

        decompresses the data and returns both compressed and decompressed forms
        """
        cdata = self.get_cdata(header)
        # now do the actual decompression
        decompressor = zlib.decompressobj(wbits=-15)  # we've alread read the header, so ignore it
        decompressed = decompressor.decompress(cdata)
        assert not decompressor.unconsumed_tail, f"unconsumed tail of {len(decompressor.unconsumed_tail)}"
        assert not decompressor.unused_data, f"unused data present of {len(decompressor.unused_data)}"
        return cdata, decompressed

    def get_block(
        self, header: Optional[Tuple[int, ...]] = None
    ) -> Tuple[Tuple[int, ...], bytes, bytes, Tuple[int, int]]:
        """
        reads the block at the current point in the file
        includes decompression and validation

        optionally accepts a header that has already been read from the file
        e.g. when scanning for the next block
        """
        if not header:
            header = self.get_header()
        cdata, decompressed = self.get_cdata_decompressed(header)
        tail = self.get_tail(decompressed)
        return header, cdata, decompressed, tail

    def get_block_lines(
        self, header: Optional[Tuple[int, ...]] = None
    ) -> Tuple[Tuple[int, ...], bytes, bytes, bytes, Tuple[bytes, ...], bytes, Tuple[int, int]]:
        """
        reads the block at the current point in the file
        includes decompression and validation

        will split by line and return the partial first and last lines, as well as the complete lines in between
        note these are returned as bytes not true string objects

        optionally accepts a header that has already been read from the file
        e.g. when scanning for the next block
        """
        start = self.raw.tell()
        header, cdata, decompressed, tail = self.get_block(header)
        # empty block
        if len(decompressed) == 0:
            return header, cdata, decompressed, b"", (), b"", tail

        # line endings can abut block ending so keep them
        decompressedlines = decompressed.splitlines(keepends=True)
        if start == 0:
            # first block has no partial start line
            firstline = b""
            lines = decompressedlines[0:-1]
            lastline = decompressedlines[-1]
        else:
            firstline = decompressedlines[0]
            lines = decompressedlines[1:-1]
            lastline = decompressedlines[-1]
        return header, cdata, decompressed, firstline, tuple(lines), lastline, tail

    def get_block_lines_offset(self, header: Optional[Tuple[int, ...]] = None) -> Tuple[
        Tuple[int, ...],
        bytes,
        bytes,
        bytes,
        Tuple[bytes, ...],
        bytes,
        Tuple[int, ...],
        Tuple[int, ...],
        Tuple[int, int],
    ]:
        """
        reads the block at the current point in the file
        includes decompression and validation

        will split by line and return the partial first and last lines, as well as the complete lines in between
        note these are returned as bytes not true string objects

        also returns the offsets within the block of start and end of each line (inclusive, including separator)

        optionally accepts a header that has already been read from the file
        e.g. when scanning for the next block
        """
        (
            header,
            cdata,
            decompressed,
            firstline,
            lines,
            lastline,
            tail,
        ) = self.get_block_lines(header)
        offsetstarts = []
        offsetends = []
        offset = len(firstline)
        for line in lines:
            offsetstarts.append(offset)
            offset += len(line)
            offsetends.append(offset)
        return (
            header,
            cdata,
            decompressed,
            firstline,
            lines,
            lastline,
            tuple(offsetstarts),
            tuple(offsetends),
            tail,
        )

    def scan_block_lines_offset(self, end: int = -1) -> Tuple[
        int,
        int,
        Tuple[int, ...],
        bytes,
        bytes,
        bytes,
        Tuple[bytes, ...],
        bytes,
        Tuple[int, ...],
        Tuple[int, ...],
        Tuple[int, int],
    ]:
        """
        starting from the current position, scan forward through the file for the next block start

        will read the block and leave the file pointing at the start of the next block
        """
        buffer = b""
        blockstart = self.raw.tell()

        # as long as there is more file to read
        while end < 0 or blockstart < end:
            # populate buffer
            if len(buffer) < headersize:
                bytesread = self.raw.read(headersize - len(buffer))
                buffer = buffer + bytesread
            # check not at end
            if len(buffer) < headersize:
                logger.warning(f"Unable to read up to {headersize}")
                raise EOFError()

            header = struct.unpack(headerpattern, buffer)
            # this is a valid location for a block
            if not self.check_is_header(header):
                # move ahead a byte
                buffer = buffer[1:]
                blockstart += 1
            else:
                # this is a valid location for a block
                (
                    header,
                    cdata,
                    decompressed,
                    firstline,
                    lines,
                    lastline,
                    offsetstarts,
                    offsetends,
                    tail,
                ) = self.get_block_lines_offset(header)
                blockend = self.raw.tell()
                return (
                    blockstart,
                    blockend,
                    header,
                    cdata,
                    decompressed,
                    firstline,
                    lines,
                    lastline,
                    offsetstarts,
                    offsetends,
                    tail,
                )
        # reach the end of the file without finding a block
        raise EOFError()

    def generate_lines_offset(self, end: int = -1) -> Generator[Tuple[int, int, int, int, bytes], None, None]:
        """
        starting from the current position, scan forward through the file
        generator that yields for each line in the file
        will stop at the end of the block that includes the end point, or
        at an empty block that indicates the end of the file
        """
        partialline = b""
        blockstart_previous = 0
        offsetstart_previous = 0
        more = True
        while more:
            (
                blockstart,
                _,
                _,
                _,
                decompressed,
                firstline,
                lines,
                lastline,
                offsetstarts,
                offsetends,
                _,
            ) = self.scan_block_lines_offset(end=end)

            if not decompressed:
                # empty block is end of file
                more = False
                # process the last partial line first
                blockstarts: Tuple[int, ...] = (blockstart_previous,)
                blockends: Tuple[int, ...] = (blockstart_previous,)
                offsetstarts = (offsetstart_previous,)
                offsetends = (offsetstart_previous + len(partialline),)
                lines = (partialline,)
            else:
                blockstarts = (blockstart_previous,) + ((blockstart,) * len(lines))
                blockends = (blockstart,) + ((blockstart,) * len(lines))
                offsetstarts = (offsetstart_previous,) + offsetstarts
                offsetends = (len(partialline),) + offsetends
                # append holdover partial line to initial line to make a new line
                lines = (partialline + firstline,) + lines
                # keep the last partial line for next block
                # last block ended on a line ending no rollover needed
                if lastline.endswith(b"\n"):
                    blockstarts = blockstarts + (blockstart,)
                    blockends = blockends + (blockstart,)
                    offsetstarts = offsetstarts + (offsetends[-1],)
                    offsetends = offsetends + (offsetends[-1] + len(lastline),)
                    lines = lines + (lastline,)
                    partialline = b""
                else:
                    partialline = lastline

            for line, start_block, start_offset, end_block, end_offset in zip(
                lines, blockstarts, offsetstarts, blockends, offsetends, strict=False
            ):
                if line:
                    yield start_block, start_offset, end_block, end_offset, line

            # prepare to do next block
            blockstart_previous = blockstart
            offsetstart_previous = offsetends[-1] + 1

    def generate_lines(self, end: int = -1) -> Generator[bytes, None, None]:
        for _, _, _, _, line in self.generate_lines_offset(end):
            yield line


class BlockGZipWriter(io.BufferedIOBase):
    # buffer size is 64kb which is 1 block
    # 65280 is what bgzip uses, for some reason?
    def __init__(self, raw: io.RawIOBase, block_size: int = 65536 - 256):
        self.raw = raw
        assert self.raw.writable()
        self.block_size = block_size
        self.block_buffer = b""

    def write(self, data: Buffer) -> int:
        self.block_buffer = self.block_buffer + data
        write_size = 0
        while len(self.block_buffer) > self.block_size:
            content = self.block_buffer[: self.block_size]
            self.block_buffer = self.block_buffer[len(content) :]
            block = self.make_block(content)
            self.raw.write(block)
            write_size += len(block)
        return write_size

    def flush(self) -> None:
        while len(self.block_buffer):
            content = self.block_buffer[: self.block_size]
            self.block_buffer = self.block_buffer[len(content) :]
            block = self.make_block(content)
            self.raw.write(block)
        self.raw.flush()

    def close(self) -> None:
        self.flush()
        # add an empty block at the end
        self.raw.write(self.make_block(b""))
        self.raw.close()

    @staticmethod
    def compress_content(content: bytes) -> bytes:
        # make a new compressor each time
        compressor = zlib.compressobj(wbits=-15)
        compressed = compressor.compress(content)
        compressed = compressed + compressor.flush()

        return compressed

    @staticmethod
    def generate_header(compressed: bytes) -> bytes:
        header = [0] * 12
        header[0] = 31  # ID1
        header[1] = 139  # ID2
        header[2] = 8  # compression method
        header[3] = 4  # flags bit2 FEXTRA
        header[4] = 0  # MTIME
        header[5] = 0  # eXtra FLags
        header[6] = 255  # OS 255 is default unspecified
        header[7] = 6  # XLEN
        header[8] = 66
        header[9] = 67
        header[10] = 2
        header[11] = len(compressed) + 6 + 19
        headerbytes = struct.pack(headerpattern, *header)
        return headerbytes

    @staticmethod
    def generate_tail(content: bytes) -> bytes:
        tail_crc = zlib.crc32(content)
        tail_isize = len(content)
        tail = [tail_crc, tail_isize]
        tailbytes = struct.pack(tailpattern, *tail)
        return tailbytes

    @classmethod
    def make_block(cls, content: bytes) -> bytes:
        compressed = cls.compress_content(content)
        headerbytes = cls.generate_header(compressed)
        tailbytes = cls.generate_tail(content)
        return headerbytes + compressed + tailbytes
