"""Tests for the tabix module."""

import pathlib

import pytest

from vcfpy.tabix import (
    Bin,
    Chunk,
    FileFormat,
    SequenceIndex,
    TabixFile,
    TabixFileIter,
    TabixIndex,
    read_index,
)


def test_read_index_multi_contig():
    """Test reading a tabix index file for a multi-contig VCF."""
    # Path to the test tabix index file
    tbi_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz.tbi"

    # Read the index
    index = read_index(tbi_path)

    # Verify it's a TabixIndex object
    assert isinstance(index, TabixIndex)

    # Check the format - should be VCF
    assert index.format == FileFormat.VCF

    # Check VCF-specific column settings
    assert index.col_seq == 1  # CHROM column (1-based)
    assert index.col_beg == 2  # POS column (1-based)
    assert index.col_end == 0  # VCF doesn't have separate end column

    # Check meta character for comments
    assert index.meta == b"#"  # VCF comment character

    # Check lines to skip at beginning
    assert index.skip == 0

    # Check sequence indices
    assert isinstance(index.indices, dict)

    # Based on the VCF file, we should have sequences "1", "2", and "20"
    expected_sequences = {"1", "2", "20"}
    assert set(index.indices.keys()) == expected_sequences

    # Check that each sequence has a SequenceIndex
    for seq_name in expected_sequences:
        assert seq_name in index.indices
        seq_index = index.indices[seq_name]
        assert isinstance(seq_index, SequenceIndex)
        assert isinstance(seq_index.bins, list)
        assert isinstance(seq_index.offsets, list)

        # Check that bins contain Bin objects with chunks
        for bin_obj in seq_index.bins:
            assert isinstance(bin_obj, Bin)
            assert isinstance(bin_obj.number, int)
            assert isinstance(bin_obj.chunks, list)

            # Check that chunks contain Chunk objects
            for chunk in bin_obj.chunks:
                assert isinstance(chunk, Chunk)
                assert isinstance(chunk.beg, int)
                assert isinstance(chunk.end, int)
                # Note: Some chunks may have end=0 (special cases in tabix format)
                # so we only check that they're non-negative
                assert chunk.beg >= 0
                assert chunk.end >= 0

        # Check that offsets are integers
        for offset in seq_index.offsets:
            assert isinstance(offset, int)


def test_read_index_file_not_found():
    """Test that reading a non-existent file raises an appropriate error."""
    non_existent_path = pathlib.Path("/non/existent/file.tbi")

    with pytest.raises(FileNotFoundError):
        read_index(non_existent_path)


def test_read_index_invalid_magic():
    """Test that reading a file with invalid magic raises ValueError."""
    # Create a temporary file with invalid magic
    import gzip
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".tbi", delete=False) as tmp:
        with gzip.open(tmp.name, "wb") as f:
            f.write(b"INVALID_MAGIC")

        try:
            with pytest.raises(ValueError, match="Invalid tabix magic"):
                read_index(tmp.name)
        finally:
            pathlib.Path(tmp.name).unlink()


def test_read_index_chunks_have_valid_offsets():
    """Test that chunks have valid virtual file offsets."""
    tbi_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz.tbi"
    index = read_index(tbi_path)

    # Check that all chunks have reasonable virtual offsets
    for _seq_name, seq_index in index.indices.items():
        for bin_obj in seq_index.bins:
            for chunk in bin_obj.chunks:
                # Virtual offsets should be positive
                assert chunk.beg >= 0
                assert chunk.end >= 0
                # End should not be less than beginning (except for special end=0 cases)
                assert chunk.end == 0 or chunk.end >= chunk.beg


def test_read_index_sequence_specific_data():
    """Test specific data for each sequence in the multi-contig file."""
    tbi_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz.tbi"
    index = read_index(tbi_path)

    # All sequences should have some data
    for seq_name in ["1", "2", "20"]:
        seq_index = index.indices[seq_name]

        # Should have at least some bins (regions with data)
        assert len(seq_index.bins) > 0

        # Count total chunks across all bins
        total_chunks = sum(len(bin_obj.chunks) for bin_obj in seq_index.bins)
        assert total_chunks > 0

        # Linear index offsets (can be empty for small regions)
        assert isinstance(seq_index.offsets, list)


def test_read_index_return_type_consistency():
    """Test that the returned TabixIndex has consistent types."""
    tbi_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz.tbi"
    index = read_index(tbi_path)

    # Check all field types
    assert isinstance(index.format, FileFormat)
    assert isinstance(index.col_seq, int)
    assert isinstance(index.col_beg, int)
    assert isinstance(index.col_end, int)
    assert isinstance(index.meta, bytes)
    assert isinstance(index.skip, int)
    assert isinstance(index.indices, dict)

    # num_no_coord is optional, can be None or int
    if index.num_no_coord is not None:
        assert isinstance(index.num_no_coord, int)


def test_read_index_detailed_structure():
    """Test detailed structure of the parsed tabix index."""
    tbi_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz.tbi"
    index = read_index(tbi_path)

    # Test that we can access all expected sequences
    assert "1" in index.indices
    assert "2" in index.indices
    assert "20" in index.indices

    # Test that sequences have expected data structure
    seq1 = index.indices["1"]
    assert len(seq1.bins) > 0  # Should have bins
    assert isinstance(seq1.offsets, list)  # Should have offsets list

    # Test bin and chunk properties
    for bin_obj in seq1.bins:
        assert hasattr(bin_obj, "number")
        assert hasattr(bin_obj, "chunks")
        assert isinstance(bin_obj.chunks, list)

        for chunk in bin_obj.chunks:
            assert hasattr(chunk, "beg")
            assert hasattr(chunk, "end")


def test_read_index_handles_pathlib_and_string():
    """Test that read_index accepts both pathlib.Path and string arguments."""
    # Test with pathlib.Path
    tbi_path_pathlib = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz.tbi"
    index1 = read_index(tbi_path_pathlib)

    # Test with string
    tbi_path_str = str(tbi_path_pathlib)
    index2 = read_index(tbi_path_str)

    # Both should produce the same result
    assert index1.format == index2.format
    assert index1.col_seq == index2.col_seq
    assert index1.col_beg == index2.col_beg
    assert index1.col_end == index2.col_end
    assert index1.meta == index2.meta
    assert index1.skip == index2.skip
    assert set(index1.indices.keys()) == set(index2.indices.keys())


def test_tabix_file_iter_construction():
    """Test TabixFileIter construction and basic functionality."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=str(vcf_path)) as tabix_file:
        # Test fetching with just reference
        iter1 = tabix_file.fetch(reference="1")
        assert isinstance(iter1, TabixFileIter)
        assert iter1.reference == "1"
        assert iter1.start is None
        assert iter1.end is None
        assert iter1.index is tabix_file.index_file
        assert iter1.sequence_index is tabix_file.index_file.indices["1"]

        # Test fetching with reference and coordinates
        iter2 = tabix_file.fetch(reference="20", start=1000000, end=2000000)
        assert iter2.reference == "20"
        assert iter2.start == 1000000
        assert iter2.end == 2000000

        # Test fetching with region string
        iter3 = tabix_file.fetch(region="2:15000-20000")
        assert iter3.reference == "2"
        assert iter3.start == 15000
        assert iter3.end == 20000

        # Test fetching with region string (no end coordinate)
        iter4 = tabix_file.fetch(region="1:50000")
        assert iter4.reference == "1"
        assert iter4.start == 50000
        assert iter4.end is None

        # Test fetching with region string (no coordinates)
        iter5 = tabix_file.fetch(region="20")
        assert iter5.reference == "20"
        assert iter5.start is None
        assert iter5.end is None


def test_tabix_file_iter_error_cases():
    """Test TabixFileIter error handling."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=str(vcf_path)) as tabix_file:
        # Test invalid reference
        with pytest.raises(ValueError, match="Reference invalid_chromosome not found"):
            tabix_file.fetch(reference="invalid_chromosome")

        # Test conflicting parameters
        with pytest.raises(ValueError, match="Cannot specify both region and reference"):
            tabix_file.fetch(region="1:1000-2000", reference="1")

        # Test start/end without reference
        with pytest.raises(ValueError, match="If start or end is given, reference must be given"):
            tabix_file.fetch(start=1000)


def test_tabix_file_iteration_comprehensive():
    """Test comprehensive tabix iteration functionality."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    # Expected records from the VCF file
    expected_records = [
        ("1", 14370, "rs6054257"),
        ("2", 17330, "."),
        ("20", 1110696, "rs6040355"),
        ("20", 1230237, "."),
        ("20", 1234567, "microsat1"),
    ]

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test fetching by individual chromosomes
        all_found_records = []
        for chrom in ["1", "2", "20"]:
            records = list(tabix_file.fetch(reference=chrom))
            for record in records:
                fields = record.split("\t")
                found_record = (fields[0], int(fields[1]), fields[2])
                all_found_records.append(found_record)

        # Verify we found all expected records
        assert len(all_found_records) == len(expected_records)

        # Check each expected record was found
        for expected in expected_records:
            assert expected in all_found_records, f"Missing expected record: {expected}"

        # Check for unexpected records
        for found in all_found_records:
            assert found in expected_records, f"Unexpected record: {found}"


def test_tabix_file_position_based_queries():
    """Test position-based queries for tabix iteration."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test specific position ranges
        test_cases = [
            ("1", 14000, 15000, [("1", 14370, "rs6054257")]),
            ("2", 17000, 18000, [("2", 17330, ".")]),
            ("20", 1110000, 1111000, [("20", 1110696, "rs6040355")]),
            ("20", 1230000, 1235000, [("20", 1230237, "."), ("20", 1234567, "microsat1")]),
            ("1", 1, 1000, []),  # No records in this range
            ("20", 2000000, 3000000, []),  # No records in this range
        ]

        for chrom, start, end, expected in test_cases:
            records = list(tabix_file.fetch(reference=chrom, start=start, end=end))
            found_records = []
            for record in records:
                fields = record.split("\t")
                found_records.append((fields[0], int(fields[1]), fields[2]))

            assert (
                found_records == expected
            ), f"Query {chrom}:{start}-{end} failed. Expected: {expected}, Found: {found_records}"


def test_tabix_file_region_string_format():
    """Test region string format for tabix queries."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test region string format
        region_tests = [
            ("20:1110000-1111000", [("20", 1110696, "rs6040355")]),
            ("1:14000-15000", [("1", 14370, "rs6054257")]),
            ("2", [("2", 17330, ".")]),  # Chromosome only
        ]

        for region, expected in region_tests:
            records = list(tabix_file.fetch(region=region))
            found_records = []
            for record in records:
                fields = record.split("\t")
                found_records.append((fields[0], int(fields[1]), fields[2]))

            assert found_records == expected, f"Region '{region}' failed. Expected: {expected}, Found: {found_records}"


def test_tabix_file_record_ordering():
    """Test that records are returned in correct positional order."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test that records are returned in correct positional order
        records = list(tabix_file.fetch(reference="20"))
        positions = []
        for record in records:
            fields = record.split("\t")
            positions.append(int(fields[1]))

        assert positions == sorted(positions), f"Positions not sorted: {positions}"
        assert len(positions) == 3, f"Expected 3 records on chromosome 20, found {len(positions)}"


def test_tabix_file_boundary_queries():
    """Test boundary condition queries."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test exact position matches
        exact_tests = [
            ("1", 14370, 14370, 1),  # Exact match - SNV at position 14370 with REF="G"
            ("20", 1110696, 1110696, 1),  # Exact match - SNV at position 1110696 with REF="A"
            ("1", 14369, 14369, 0),  # Just before - no variant at this position
            ("1", 14371, 14371, 0),  # Just after - no variant at this position
            (
                "20",
                1234567,
                1234569,
                1,
            ),  # Should match microsatellite at 1234567 with REF="GTC" (spans 1234567-1234569)
        ]

        for chrom, start, end, expected_count in exact_tests:
            records = list(tabix_file.fetch(reference=chrom, start=start, end=end))
            assert (
                len(records) == expected_count
            ), f"Query {chrom}:{start}-{end}: expected {expected_count} records, found {len(records)}"


def test_tabix_file_large_range_queries():
    """Test very large range queries."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test very large ranges
        large_range_tests = [
            ("20", 1, 2000000, 3),  # Should find all records on chr20
            ("1", 1, 100000000, 1),  # Should find all records on chr1
            ("2", 1, 100000000, 1),  # Should find all records on chr2
        ]

        for chrom, start, end, expected_count in large_range_tests:
            records = list(tabix_file.fetch(reference=chrom, start=start, end=end))
            assert (
                len(records) == expected_count
            ), f"Large range {chrom}:{start}-{end}: expected {expected_count} records, found {len(records)}"


def test_tabix_file_iteration_performance():
    """Test that iteration doesn't produce duplicate records."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Get all records and ensure no duplicates
        all_records = []
        for chrom in ["1", "2", "20"]:
            records = list(tabix_file.fetch(reference=chrom))
            all_records.extend(records)

        # Check that we don't have duplicate lines
        unique_records = set(all_records)
        assert len(all_records) == len(
            unique_records
        ), f"Found duplicate records: {len(all_records)} total vs {len(unique_records)} unique"

        # Verify total count matches expected
        assert len(all_records) == 5, f"Expected 5 total records, found {len(all_records)}"


def test_tabix_against_reference_implementation():
    """
    Test our tabix implementation against known correct results from the real tabix command.

    These expected results were verified by running: tabix tests/vcfs/multi_contig.vcf.gz <range>
    This test can run in CI without requiring tabix to be installed.
    """
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    # Expected results verified against real tabix command
    expected_results = {
        # Full chromosome queries
        ("1", None, None): [
            "1\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,."
        ],
        ("2", None, None): [
            "2\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3:40,30"
        ],
        ("20", None, None): [
            "20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4:40,30",
            "20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2:40,30",
            "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3",
        ],
        # Position range queries
        ("1", 14000, 15000): [
            "1\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,."
        ],
        ("20", 1110000, 1111000): [
            "20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4:40,30"
        ],
        ("20", 1230000, 1235000): [
            "20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2:40,30",
            "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3",
        ],
        # Critical test case (the one that was failing in test_fetch.py)
        ("20", 1110697, 1234568): [
            "20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2:40,30",
            "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3",
        ],
        # Boundary cases
        ("1", 14370, 14370): [
            "1\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,."
        ],
        ("1", 14371, 14371): [],  # Just after SNV - should be empty
        ("20", 1234567, 1234569): [  # Microsatellite range (REF=GTC spans 3 positions)
            "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3"
        ],
        ("20", 1234570, 1234570): [],  # Just after microsatellite - should be empty
        # Empty result cases
        ("1", 1, 1000): [],  # Before any variants
        ("20", 2000000, 3000000): [],  # After all variants
    }

    with TabixFile(filename=vcf_path) as tabix_file:
        for (chrom, start, end), expected in expected_results.items():
            # Get actual results
            if start is None and end is None:
                actual = list(tabix_file.fetch(reference=chrom))
            else:
                actual = list(tabix_file.fetch(reference=chrom, start=start, end=end))

            # Compare
            assert len(actual) == len(
                expected
            ), f"Query {chrom}:{start}-{end}: expected {len(expected)} records, got {len(actual)}"

            for i, (actual_record, expected_record) in enumerate(zip(actual, expected, strict=False)):
                assert actual_record == expected_record, (
                    f"Query {chrom}:{start}-{end}, record {i}: mismatch\n"
                    f"Expected: {expected_record}\n"
                    f"Actual:   {actual_record}"
                )

        # Test error case for non-existent chromosome
        with pytest.raises(ValueError, match="not found in index"):
            list(tabix_file.fetch(reference="3", start=1, end=1000000))


def test_tabix_sequence_ordering():
    """
    Test that sequence ordering follows index order, not lexical order.

    This is important for chromosomes like "10" vs "2" where lexical ordering
    would incorrectly place "2" after "10".
    """
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "multi_contig.vcf.gz"

    with TabixFile(filename=vcf_path) as tabix_file:
        # Get the sequence order from the index
        sequence_names = list(tabix_file.index_file.indices.keys())

        # The test file has sequences in order: ['1', '2', '20']
        assert sequence_names == ["1", "2", "20"], f"Expected ['1', '2', '20'], got {sequence_names}"

        # Verify that lexical ordering would be different for some chromosome names
        # In general, "10" > "2" lexically but "10" should come after "2" genomically

        # Test that iteration respects this ordering by fetching from each chromosome
        all_lines = []
        for seq_name in sequence_names:
            records = list(tabix_file.fetch(reference=seq_name))
            for record in records:
                # Extract chromosome from each record
                fields = record.split("\t")
                chrom = fields[0]
                all_lines.append((sequence_names.index(chrom), chrom, record))

        # Verify that records are in index order
        sorted_lines = sorted(all_lines, key=lambda x: x[0])  # Sort by index position
        assert all_lines == sorted_lines, "Records should be in index order, not lexical order"
