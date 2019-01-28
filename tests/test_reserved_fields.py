# -*- coding: utf-8 -*-
"""Tests for reserved field handling in header
"""

from vcfpy import header

# INFO for (small) variants ---------------------------------------------------


def test_info_aa():
    assert header.RESERVED_INFO["AA"]
    assert header.RESERVED_INFO["AA"].type == "String"
    assert header.RESERVED_INFO["AA"].number == 1
    assert header.RESERVED_INFO["AA"].description


def test_info_ac():
    assert header.RESERVED_INFO["AC"]
    assert header.RESERVED_INFO["AC"].type == "Integer"
    assert header.RESERVED_INFO["AC"].number == "A"
    assert header.RESERVED_INFO["AC"].description


def test_info_ad():
    assert header.RESERVED_INFO["AD"]
    assert header.RESERVED_INFO["AD"].type == "Integer"
    assert header.RESERVED_INFO["AD"].number == "R"
    assert header.RESERVED_INFO["AD"].description


def test_info_adf():
    assert header.RESERVED_INFO["ADF"]
    assert header.RESERVED_INFO["ADF"].type == "Integer"
    assert header.RESERVED_INFO["ADF"].number == "R"
    assert header.RESERVED_INFO["ADF"].description


def test_info_adr():
    assert header.RESERVED_INFO["ADR"]
    assert header.RESERVED_INFO["ADR"].type == "Integer"
    assert header.RESERVED_INFO["ADR"].number == "R"
    assert header.RESERVED_INFO["ADR"].description


def test_info_an():
    assert header.RESERVED_INFO["AN"]
    assert header.RESERVED_INFO["AN"].type == "Integer"
    assert header.RESERVED_INFO["AN"].number == 1
    assert header.RESERVED_INFO["AN"].description


def test_info_bq():
    assert header.RESERVED_INFO["BQ"]
    assert header.RESERVED_INFO["BQ"].type == "Float"
    assert header.RESERVED_INFO["BQ"].number == 1
    assert header.RESERVED_INFO["BQ"].description


def test_info_cigar():
    assert header.RESERVED_INFO["CIGAR"]
    assert header.RESERVED_INFO["CIGAR"].type == "String"
    assert header.RESERVED_INFO["CIGAR"].number == "A"
    assert header.RESERVED_INFO["CIGAR"].description


def test_info_DP():
    assert header.RESERVED_INFO["DB"]
    assert header.RESERVED_INFO["DB"].type == "Flag"
    assert header.RESERVED_INFO["DB"].number == 0
    assert header.RESERVED_INFO["DB"].description


def test_info_dp():
    assert header.RESERVED_INFO["DP"]
    assert header.RESERVED_INFO["DP"].type == "Integer"
    assert header.RESERVED_INFO["DP"].number == 1
    assert header.RESERVED_INFO["DP"].description


def test_info_h2():
    assert header.RESERVED_INFO["H2"]
    assert header.RESERVED_INFO["H2"].type == "Flag"
    assert header.RESERVED_INFO["H2"].number == 0
    assert header.RESERVED_INFO["H2"].description


def test_info_h3():
    assert header.RESERVED_INFO["H3"]
    assert header.RESERVED_INFO["H3"].type == "Flag"
    assert header.RESERVED_INFO["H3"].number == 0
    assert header.RESERVED_INFO["H3"].description


def test_info_mq():
    assert header.RESERVED_INFO["MQ"]
    assert header.RESERVED_INFO["MQ"].type == "Integer"
    assert header.RESERVED_INFO["MQ"].number == 1
    assert header.RESERVED_INFO["MQ"].description


def test_info_mq0():
    assert header.RESERVED_INFO["MQ0"]
    assert header.RESERVED_INFO["MQ0"].type == "Integer"
    assert header.RESERVED_INFO["MQ0"].number == 1
    assert header.RESERVED_INFO["MQ0"].description


def test_info_ns():
    assert header.RESERVED_INFO["NS"]
    assert header.RESERVED_INFO["NS"].type == "Integer"
    assert header.RESERVED_INFO["NS"].number == 1
    assert header.RESERVED_INFO["NS"].description


def test_info_sb():
    assert header.RESERVED_INFO["SB"]
    assert header.RESERVED_INFO["SB"].type == "Integer"
    assert header.RESERVED_INFO["SB"].number == 4
    assert header.RESERVED_INFO["SB"].description


def test_info_somatic():
    assert header.RESERVED_INFO["SOMATIC"]
    assert header.RESERVED_INFO["SOMATIC"].type == "Flag"
    assert header.RESERVED_INFO["SOMATIC"].number == 0
    assert header.RESERVED_INFO["SOMATIC"].description


def test_info_validated():
    assert header.RESERVED_INFO["VALIDATED"]
    assert header.RESERVED_INFO["VALIDATED"].type == "Flag"
    assert header.RESERVED_INFO["VALIDATED"].number == 0
    assert header.RESERVED_INFO["VALIDATED"].description


def test_info_1000g():
    assert header.RESERVED_INFO["1000G"]
    assert header.RESERVED_INFO["1000G"].type == "Flag"
    assert header.RESERVED_INFO["1000G"].number == 0
    assert header.RESERVED_INFO["1000G"].description


# INFO for SVs ----------------------------------------------------------------


def test_info_imprecise():
    assert header.RESERVED_INFO["IMPRECISE"]
    assert header.RESERVED_INFO["IMPRECISE"].type == "Flag"
    assert header.RESERVED_INFO["IMPRECISE"].number == 0
    assert header.RESERVED_INFO["IMPRECISE"].description


def test_info_novel():
    assert header.RESERVED_INFO["NOVEL"]
    assert header.RESERVED_INFO["NOVEL"].type == "Flag"
    assert header.RESERVED_INFO["NOVEL"].number == 0
    assert header.RESERVED_INFO["NOVEL"].description


def test_info_end():
    assert header.RESERVED_INFO["END"]
    assert header.RESERVED_INFO["END"].type == "Integer"
    assert header.RESERVED_INFO["END"].number == 1
    assert header.RESERVED_INFO["END"].description


def test_info_svtype():
    assert header.RESERVED_INFO["SVTYPE"]
    assert header.RESERVED_INFO["SVTYPE"].type == "String"
    assert header.RESERVED_INFO["SVTYPE"].number == 1
    assert header.RESERVED_INFO["SVTYPE"].description


def test_info_svlen():
    assert header.RESERVED_INFO["SVLEN"]
    assert header.RESERVED_INFO["SVLEN"].type == "Integer"
    assert header.RESERVED_INFO["SVLEN"].number == 1
    assert header.RESERVED_INFO["SVLEN"].description


def test_info_cipos():
    assert header.RESERVED_INFO["CIPOS"]
    assert header.RESERVED_INFO["CIPOS"].type == "Integer"
    assert header.RESERVED_INFO["CIPOS"].number == 2
    assert header.RESERVED_INFO["CIPOS"].description


def test_info_ciend():
    assert header.RESERVED_INFO["CIEND"]
    assert header.RESERVED_INFO["CIEND"].type == "Integer"
    assert header.RESERVED_INFO["CIEND"].number == 2
    assert header.RESERVED_INFO["CIEND"].description


def test_info_homlen():
    assert header.RESERVED_INFO["HOMLEN"]
    assert header.RESERVED_INFO["HOMLEN"].type == "Integer"
    assert header.RESERVED_INFO["HOMLEN"].number == "."
    assert header.RESERVED_INFO["HOMLEN"].description


def test_info_homseq():
    assert header.RESERVED_INFO["HOMSEQ"]
    assert header.RESERVED_INFO["HOMSEQ"].type == "String"
    assert header.RESERVED_INFO["HOMSEQ"].number == "."
    assert header.RESERVED_INFO["HOMSEQ"].description


def test_info_bkptid():
    assert header.RESERVED_INFO["BKPTID"]
    assert header.RESERVED_INFO["BKPTID"].type == "String"
    assert header.RESERVED_INFO["BKPTID"].number == "."
    assert header.RESERVED_INFO["BKPTID"].description


def test_info_meinfo():
    assert header.RESERVED_INFO["MEINFO"]
    assert header.RESERVED_INFO["MEINFO"].type == "String"
    assert header.RESERVED_INFO["MEINFO"].number == 4
    assert header.RESERVED_INFO["MEINFO"].description


def test_info_metrans():
    assert header.RESERVED_INFO["METRANS"]
    assert header.RESERVED_INFO["METRANS"].type == "String"
    assert header.RESERVED_INFO["METRANS"].number == 4
    assert header.RESERVED_INFO["METRANS"].description


def test_info_dgvid():
    assert header.RESERVED_INFO["DGVID"]
    assert header.RESERVED_INFO["DGVID"].type == "String"
    assert header.RESERVED_INFO["DGVID"].number == 1
    assert header.RESERVED_INFO["DGVID"].description


def test_info_dbvarid():
    assert header.RESERVED_INFO["DBVARID"]
    assert header.RESERVED_INFO["DBVARID"].type == "String"
    assert header.RESERVED_INFO["DBVARID"].number == 1
    assert header.RESERVED_INFO["DBVARID"].description


def test_info_dbripid():
    assert header.RESERVED_INFO["DBRIPID"]
    assert header.RESERVED_INFO["DBRIPID"].type == "String"
    assert header.RESERVED_INFO["DBRIPID"].number == 1
    assert header.RESERVED_INFO["DBRIPID"].description


def test_info_mateid():
    assert header.RESERVED_INFO["MATEID"]
    assert header.RESERVED_INFO["MATEID"].type == "String"
    assert header.RESERVED_INFO["MATEID"].number == "."
    assert header.RESERVED_INFO["MATEID"].description


def test_info_parid():
    assert header.RESERVED_INFO["PARID"]
    assert header.RESERVED_INFO["PARID"].type == "String"
    assert header.RESERVED_INFO["PARID"].number == 1
    assert header.RESERVED_INFO["PARID"].description


def test_info_event():
    assert header.RESERVED_INFO["EVENT"]
    assert header.RESERVED_INFO["EVENT"].type == "String"
    assert header.RESERVED_INFO["EVENT"].number == 1
    assert header.RESERVED_INFO["EVENT"].description


def test_info_cilen():
    assert header.RESERVED_INFO["CILEN"]
    assert header.RESERVED_INFO["CILEN"].type == "Integer"
    assert header.RESERVED_INFO["CILEN"].number == 2
    assert header.RESERVED_INFO["CILEN"].description


def test_info_dpadj():
    assert header.RESERVED_INFO["DPADJ"]
    assert header.RESERVED_INFO["DPADJ"].type == "Integer"
    assert header.RESERVED_INFO["DPADJ"].number == "."
    assert header.RESERVED_INFO["DPADJ"].description


def test_info_cn():
    assert header.RESERVED_INFO["CN"]
    assert header.RESERVED_INFO["CN"].type == "Integer"
    assert header.RESERVED_INFO["CN"].number == 1
    assert header.RESERVED_INFO["CN"].description


def test_info_cnadj():
    assert header.RESERVED_INFO["CNADJ"]
    assert header.RESERVED_INFO["CNADJ"].type == "Integer"
    assert header.RESERVED_INFO["CNADJ"].number == "."
    assert header.RESERVED_INFO["CNADJ"].description


def test_info_cicn():
    assert header.RESERVED_INFO["CICN"]
    assert header.RESERVED_INFO["CICN"].type == "Integer"
    assert header.RESERVED_INFO["CICN"].number == 2
    assert header.RESERVED_INFO["CICN"].description


def test_info_cicnadj():
    assert header.RESERVED_INFO["CICNADJ"]
    assert header.RESERVED_INFO["CICNADJ"].type == "Integer"
    assert header.RESERVED_INFO["CICNADJ"].number == "."
    assert header.RESERVED_INFO["CICNADJ"].description


# FORMAT for (small) variants -------------------------------------------------


def test_format_ad():
    assert header.RESERVED_FORMAT["AD"]
    assert header.RESERVED_FORMAT["AD"].type == "Integer"
    assert header.RESERVED_FORMAT["AD"].number == "R"
    assert header.RESERVED_FORMAT["AD"].description


def test_format_adf():
    assert header.RESERVED_FORMAT["ADF"]
    assert header.RESERVED_FORMAT["ADF"].type == "Integer"
    assert header.RESERVED_FORMAT["ADF"].number == "R"
    assert header.RESERVED_FORMAT["ADF"].description


def test_format_adr():
    assert header.RESERVED_FORMAT["ADR"]
    assert header.RESERVED_FORMAT["ADR"].type == "Integer"
    assert header.RESERVED_FORMAT["ADR"].number == "R"
    assert header.RESERVED_FORMAT["ADR"].description


def test_format_dp():
    assert header.RESERVED_FORMAT["DP"]
    assert header.RESERVED_FORMAT["DP"].type == "Integer"
    assert header.RESERVED_FORMAT["DP"].number == 1
    assert header.RESERVED_FORMAT["DP"].description


def test_format_ec():
    assert header.RESERVED_FORMAT["EC"]
    assert header.RESERVED_FORMAT["EC"].type == "Integer"
    assert header.RESERVED_FORMAT["EC"].number == "A"
    assert header.RESERVED_FORMAT["EC"].description


def test_format_ft():
    assert header.RESERVED_FORMAT["FT"]
    assert header.RESERVED_FORMAT["FT"].type == "String"
    assert header.RESERVED_FORMAT["FT"].number == "1"
    assert header.RESERVED_FORMAT["FT"].description


def test_format_gq():
    assert header.RESERVED_FORMAT["GQ"]
    assert header.RESERVED_FORMAT["GQ"].type == "Integer"
    assert header.RESERVED_FORMAT["GQ"].number == "G"
    assert header.RESERVED_FORMAT["GQ"].description


def test_format_gp():
    assert header.RESERVED_FORMAT["GP"]
    assert header.RESERVED_FORMAT["GP"].type == "Float"
    assert header.RESERVED_FORMAT["GP"].number == "G"
    assert header.RESERVED_FORMAT["GP"].description


def test_format_gt():
    assert header.RESERVED_FORMAT["GT"]
    assert header.RESERVED_FORMAT["GT"].type == "String"
    assert header.RESERVED_FORMAT["GT"].number == 1
    assert header.RESERVED_FORMAT["GT"].description


def test_format_gl():
    assert header.RESERVED_FORMAT["GL"]
    assert header.RESERVED_FORMAT["GL"].type == "Float"
    assert header.RESERVED_FORMAT["GL"].number == "G"
    assert header.RESERVED_FORMAT["GL"].description


def test_format_hq():
    assert header.RESERVED_FORMAT["HQ"]
    assert header.RESERVED_FORMAT["HQ"].type == "Integer"
    assert header.RESERVED_FORMAT["HQ"].number == 2
    assert header.RESERVED_FORMAT["HQ"].description


def test_format_mq():
    assert header.RESERVED_FORMAT["MQ"]
    assert header.RESERVED_FORMAT["MQ"].type == "Integer"
    assert header.RESERVED_FORMAT["MQ"].number == 1
    assert header.RESERVED_FORMAT["MQ"].description


def test_format_pl():
    assert header.RESERVED_FORMAT["PL"]
    assert header.RESERVED_FORMAT["PL"].type == "Integer"
    assert header.RESERVED_FORMAT["PL"].number == "G"
    assert header.RESERVED_FORMAT["PL"].description


def test_format_pq():
    assert header.RESERVED_FORMAT["PQ"]
    assert header.RESERVED_FORMAT["PQ"].type == "Integer"
    assert header.RESERVED_FORMAT["PQ"].number == 1
    assert header.RESERVED_FORMAT["PQ"].description


def test_format_ps():
    assert header.RESERVED_FORMAT["PS"]
    assert header.RESERVED_FORMAT["PS"].type == "Integer"
    assert header.RESERVED_FORMAT["PS"].number == 1
    assert header.RESERVED_FORMAT["PS"].description


# FORMAT for SVs --------------------------------------------------------------


def test_format_cn():
    assert header.RESERVED_FORMAT["CN"]
    assert header.RESERVED_FORMAT["CN"].type == "Integer"
    assert header.RESERVED_FORMAT["CN"].number == 1
    assert header.RESERVED_FORMAT["CN"].description


def test_format_cnq():
    assert header.RESERVED_FORMAT["CNQ"]
    assert header.RESERVED_FORMAT["CNQ"].type == "Float"
    assert header.RESERVED_FORMAT["CNQ"].number == 1
    assert header.RESERVED_FORMAT["CNQ"].description


def test_format_cnl():
    assert header.RESERVED_FORMAT["CNL"]
    assert header.RESERVED_FORMAT["CNL"].type == "Float"
    assert header.RESERVED_FORMAT["CNL"].number == "G"
    assert header.RESERVED_FORMAT["CNL"].description


def test_format_cnp():
    assert header.RESERVED_FORMAT["CNP"]
    assert header.RESERVED_FORMAT["CNP"].type == "Float"
    assert header.RESERVED_FORMAT["CNP"].number == "G"
    assert header.RESERVED_FORMAT["CNP"].description


def test_format_nq():
    assert header.RESERVED_FORMAT["NQ"]
    assert header.RESERVED_FORMAT["NQ"].type == "Integer"
    assert header.RESERVED_FORMAT["NQ"].number == 1
    assert header.RESERVED_FORMAT["NQ"].description


def test_format_hap():
    assert header.RESERVED_FORMAT["HAP"]
    assert header.RESERVED_FORMAT["HAP"].type == "Integer"
    assert header.RESERVED_FORMAT["HAP"].number == 1
    assert header.RESERVED_FORMAT["HAP"].description


def test_format_ahap():
    assert header.RESERVED_FORMAT["AHAP"]
    assert header.RESERVED_FORMAT["AHAP"].type == "Integer"
    assert header.RESERVED_FORMAT["AHAP"].number == 1
    assert header.RESERVED_FORMAT["AHAP"].description
