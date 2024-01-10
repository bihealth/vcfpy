# Changelog

## [0.13.8](https://github.com/bihealth/vcfpy/compare/v0.13.7...v0.13.8) (2024-01-10)


### Bug Fixes

* fixing manifest for changelog ([#169](https://github.com/bihealth/vcfpy/issues/169)) ([83c5b8e](https://github.com/bihealth/vcfpy/commit/83c5b8e6cd1199245673cc0d8deb2d6f3646d183))

## [0.13.7](https://github.com/bihealth/vcfpy/compare/v0.13.6...v0.13.7) (2024-01-10)


### Bug Fixes

* remove versioneer Python 3.12 compatibility ([#160](https://github.com/bihealth/vcfpy/issues/160)) ([5e2860e](https://github.com/bihealth/vcfpy/commit/5e2860e22042aa794304c8805ca716a39c88f24e))


## [0.13.6](https://github.com/bihealth/vcfpy/compare/v0.13.5...v0.13.6) (2022-11-28)

- Fixing bug in `setup.py` that prevented `pysam` dependency to be loaded (#150).

## v0.13.5 (2022-11-13)

- Treat `.bgz` files the same as `.gz` (#145, \#149)

## v0.13.4 (2022-04-13)

- Switching to Github Actions for CI
- Fix INFO flag raises TypeError (#146)

## v0.13.3 (2020-09-14)

- Adding `Record.update_calls`.
- Making `Record.{format,calls}` use list when empty

## v0.13.2 (2020-08-20)

- Adding `Call.set_genotype()`.

## v0.13.1 (2020-08-20)

- Fixed `Call.ploidy`.
- Fixed `Call.is_variant`.

## v0.13.0 (2020-07-10)

- Fixing bug in case `GT` describes only one allele.
- Proper escaping of colon and semicolon (or the lack of escaping) in
  `INFO` and `FORMAT`.

## v0.12.2 (2020-04-29)

- Fixing bug in case `GT` describes only one allele.

## v0.12.1 (2019-03-08)

- Not warning on `PASS` filter if not defined in header.

## v0.12.0 (2019-01-29)

- Fixing tests for Python \>=3.6
- Fixing CI, improving tox integration.
- Applying `black` formatting.
- Replacing Makefile with more minimal one.
- Removing some linting errors from flake8.
- Adding support for reading VCF without `FORMAT` or any sample
  column.
- Adding support for writing headers and records without `FORMAT` and
  any sample columns.

## v0.11.2 (2018-04-16)

- Removing `pip` module from `setup.py` which is not recommended
  anyway.

## v0.11.1 (2018-03-06)

- Working around problem in HTSJDK output with incomplete `FORMAT`
  fields (#127). Writing out `.` instead of keeping trailing empty
  records empty.

## v0.11.0 (2017-11-22)

- The field `FORMAT/FT` is now expected to be a semicolon-separated
  string. Internally, we will handle it as a list.
- Switching from warning helper utility code to Python `warnings`
  module.
- Return `str` in case of problems with parsing value.

## v0.10.0 (2017-02-27)

- Extending API to allow for reading subsets of records. (Writing for
  sample subsets or reordered samples is possible through using the
  appropriate `names` list in the `SamplesInfos` for the `Writer`).
- Deep-copying header lines and samples infos on `Writer` construction
- Using `samples` attribute from `Header` in `Reader` and `Writer`
  instead of passing explicitely

## 0.9.0 (2017-02-26)

- Restructuring of requirements.txt files
- Fixing parsing of no-call `GT` fields

## 0.8.1 (2017-02-08)

- PEP8 style adjustments
- Using versioneer for versioning
- Using `requirements*.txt` files now from setup.py
- Fixing dependency on cyordereddict to be for Python \<3.6 instead of
  \<3.5
- Jumping by samtools coordinate string now also allowed

## 0.8.0 (2016-10-31)

- Adding `Header.has_header_line` for querying existence of header
  line
- `Header.add_*_line` return a `bool` no indicating any conflicts
- Construction of Writer uses samples within header and no extra
  parameter (breaks API)

## 0.7.0 (2016-09-25)

- Smaller improvements and fixes to documentation
- Adding Codacy coverage and static code analysis results to README
- Various smaller code cleanup triggered by Codacy results
- Adding `__eq__`, `__neq__` and `__hash__` to data types (where
  applicable)

## 0.6.0 (2016-09-25

- Refining implementation for breakend and symbolic allele class
- Removing `record.SV_CODES`
- Refactoring parser module a bit to make the code cleaner
- Fixing small typos and problems in documentation

## 0.5.0 (2016-09-24)

- Deactivating warnings on record parsing by default because of
  performance
- Adding validation for `INFO` and `FORMAT` fields on reading (#8)
- Adding predefined `INFO` and `FORMAT` fields to `pyvcf.header` (#32)

## 0.4.1 (2016-09-22)

- Initially enabling codeclimate

## 0.4.0 (2016-09-22)

- Exporting constants for encoding variant types
- Exporting genotype constants `HOM_REF`, `HOM_ALT`, `HET`
- Implementing `Call.is_phased`, `Call.is_het`, `Call.is_variant`,
  `Call.is_phased`, `Call.is_hom_ref`, `Call.is_hom_alt`
- Removing `Call.phased` (breaks API, next release is 0.4.0)
- Adding tests, fixing bugs for methods of `Call`

## 0.3.1 (2016-09-21)

- Work around `FORMAT/FT` being a string; this is done so in the Delly
  output

## 0.3.0 (2016-09-21)

- `Reader` and `Writer` can now be used as context manager (with
  `with`)
- Including license in documentation, including Biopython license
- Adding support for writing bgzf files (taken from Biopython)
- Adding support for parsing arrays in header lines
- Removing `example-4.1-bnd.vcf` example file because v4.1 tumor
  derival lacks `ID` field
- Adding `AltAlleleHeaderLine`, `MetaHeaderLine`,
  `PedigreeHeaderLine`, and `SampleHeaderLine`
- Renaming `SimpleHeaderFile` to `SimpleHeaderLine`
- Warn on missing `FILTER` entries on parsing
- Reordered parameters in `from_stream` and `from_file` (#18)
- Renamed `from_file` to `from_stream` (#18)
- Renamed `Reader.jump_to` to `Reader.fetch`
- Adding `header_without_lines` function
- Generally extending API to make it esier to use
- Upgrading dependencies, enabling pyup-bot
- Greatly extending documentation

## 0.2.1 (2016-09-19)

- First release on PyPI
