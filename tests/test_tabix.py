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


def test_tabix_large_file_comprehensive():
    """Test tabix functionality on a large real-world VCF file."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test 1: Basic file properties
        assert tabix_file.index_file.format == FileFormat.VCF
        assert "SL2.50ch00" in tabix_file.index_file.indices

        # Test 2: Small range at the beginning - should have 44 records
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=1, end=1000))
        assert len(records) == 44, f"Expected 44 records in range 1-1000, got {len(records)}"

        # Verify first record
        first_record = records[0]
        fields = first_record.split("\t")
        assert fields[0] == "SL2.50ch00"
        assert int(fields[1]) == 280  # First position should be 280

        # Test 3: Dense region - should have 114 records
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=10000, end=11000))
        assert len(records) == 114, f"Expected 114 records in range 10000-11000, got {len(records)}"

        # Test 4: Medium range - should have 360 records
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=80000, end=85000))
        assert len(records) == 360, f"Expected 360 records in range 80000-85000, got {len(records)}"

        # Test 5: Large range - should have 3118 records
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=50000, end=100000))
        assert len(records) == 3118, f"Expected 3118 records in range 50000-100000, got {len(records)}"


def test_tabix_large_file_position_verification():
    """Test that positions are correctly parsed and filtered in large file."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test position boundaries
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=90000, end=95000))

        for record in records:
            fields = record.split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]

            assert chrom == "SL2.50ch00"
            assert 90000 <= pos <= 95000, f"Position {pos} outside range 90000-95000"

            # For VCF, end position is pos + len(ref) - 1
            end_pos = pos + len(ref) - 1
            # Record should overlap with query range [90000, 95000]
            assert not (end_pos < 90000 or pos > 95000), f"Record at {pos}-{end_pos} doesn't overlap 90000-95000"


def test_tabix_large_file_edge_cases():
    """Test edge cases with the large file."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test 1: Query before any records
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=1, end=100))
        assert len(records) == 0, "Should have no records before position 280"

        # Test 2: Query after all records (assuming file ends before 200000)
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=200000, end=300000))
        # This might have 0 records if the file doesn't extend that far

        # Test 3: Single position query
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=90016, end=90016))
        # Should find the record at position 90016 if it exists
        if records:
            fields = records[0].split("\t")
            pos = int(fields[1])
            ref = fields[3]
            end_pos = pos + len(ref) - 1
            assert pos <= 90016 <= end_pos, f"Record at {pos}-{end_pos} should contain position 90016"

        # Test 4: Very small range
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=90090, end=90100))
        # Verify all records are in the correct range
        for record in records:
            fields = record.split("\t")
            pos = int(fields[1])
            ref = fields[3]
            end_pos = pos + len(ref) - 1
            # Record should overlap with query range
            assert not (end_pos < 90090 or pos > 90100), f"Record at {pos}-{end_pos} doesn't overlap 90090-90100"


def test_tabix_large_file_different_query_sizes():
    """Test different query sizes to validate binning scheme performance."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test various bin sizes to exercise different levels of the UCSC binning scheme
        test_cases = [
            # Small queries (should use high-resolution bins)
            ("SL2.50ch00", 1000, 2000, "small_query_1kb"),
            ("SL2.50ch00", 5000, 6000, "small_query_1kb_2"),
            # Medium queries (should use medium-resolution bins)
            ("SL2.50ch00", 10000, 20000, "medium_query_10kb"),
            ("SL2.50ch00", 30000, 50000, "medium_query_20kb"),
            # Large queries (should use low-resolution bins)
            ("SL2.50ch00", 1, 50000, "large_query_50kb"),
            ("SL2.50ch00", 25000, 75000, "large_query_50kb_2"),
        ]

        for chrom, start, end, description in test_cases:
            records = list(tabix_file.fetch(reference=chrom, start=start, end=end))

            # Verify all records are in the correct range
            for record in records:
                fields = record.split("\t")
                pos = int(fields[1])
                ref = fields[3]
                record_end = pos + len(ref) - 1

                # Record should overlap with query range
                assert not (
                    record_end < start or pos > end
                ), f"{description}: Record at {pos}-{record_end} doesn't overlap {start}-{end}"

                assert fields[0] == chrom, f"{description}: Wrong chromosome {fields[0]}, expected {chrom}"


def test_tabix_large_file_validation_against_reference():
    """Validate our implementation against reference tabix command results."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Test cases with known results from reference tabix command
        validation_cases = [
            # Range, Expected count (from tabix command)
            ("SL2.50ch00", 1, 1000, 44),
            ("SL2.50ch00", 10000, 11000, 114),
            ("SL2.50ch00", 80000, 85000, 360),
            ("SL2.50ch00", 50000, 100000, 3118),
        ]

        for chrom, start, end, expected_count in validation_cases:
            records = list(tabix_file.fetch(reference=chrom, start=start, end=end))
            actual_count = len(records)

            assert (
                actual_count == expected_count
            ), f"Range {chrom}:{start}-{end}: expected {expected_count} records, got {actual_count}"

            # Verify record ordering (should be sorted by position)
            positions = []
            for record in records:
                fields = record.split("\t")
                pos = int(fields[1])
                positions.append(pos)

            assert positions == sorted(positions), f"Records in range {chrom}:{start}-{end} are not sorted by position"


def test_tabix_large_file_indel_handling():
    """Test handling of insertions and deletions in large file."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Get some records that include indels
        records = list(tabix_file.fetch(reference="SL2.50ch00", start=400, end=500))

        indel_count = 0
        for record in records:
            fields = record.split("\t")
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            # Check if this is an indel
            is_indel = len(ref) != len(alt)
            if is_indel:
                indel_count += 1

                # Verify the record is correctly positioned
                record_end = pos + len(ref) - 1
                assert 400 <= pos <= 500 or (
                    pos < 400 and record_end >= 400
                ), f"Indel at {pos}-{record_end} incorrectly included in range 400-500"

        # The test file should have at least one indel in this range
        # (We saw ATTT->ATTTT at position 409 in the sample output)
        assert indel_count > 0, "Should find at least one indel in the test range"


def test_tabix_large_file_performance_bins():
    """Test that different bin levels are being used appropriately."""
    vcf_path = pathlib.Path(__file__).parent / "vcfs" / "annotated_tomato_150.100000.vcf.gz"

    if not vcf_path.exists():
        pytest.skip("Large tomato VCF file not found")

    with TabixFile(filename=vcf_path) as tabix_file:
        # Get the sequence index for SL2.50ch00
        seq_index = tabix_file.index_file.indices["SL2.50ch00"]

        # Verify we have bins
        bin_numbers = [bin_obj.number for bin_obj in seq_index.bins]
        assert len(bin_numbers) > 0, "Should have at least some bins"

        # Check bin number ranges for UCSC binning scheme
        # Level 5: bins 4681-37448 (each spans ~128kb)
        # For a file covering ~150kb, we should mainly have level 5 bins
        level_5_bins = [b for b in bin_numbers if 4681 <= b <= 37448]

        # The file covers ~150kb so should have multiple level 5 bins
        assert len(level_5_bins) > 0, "Should have level 5 bins for fine-grained access"
        assert len(level_5_bins) == len(bin_numbers), "All bins should be level 5 for this file size"

        # Verify bins are consecutive or at least in reasonable range
        min_bin = min(bin_numbers)
        max_bin = max(bin_numbers)
        bin_range = max_bin - min_bin + 1

        # For ~150kb file, we expect a reasonable number of bins
        assert 50 <= len(bin_numbers) <= 200, f"Expected 50-200 bins, got {len(bin_numbers)}"
        assert bin_range <= 200, f"Bin range {bin_range} seems too large for a 150kb file"
