=======
History
=======

0.9.0 (2017-02-26)
------------------

* Restructuring of requirements.txt files
* Fixing parsing of no-call ``GT`` fields

0.8.1 (2017-02-08)
------------------

* PEP8 style adjustments
* Using versioneer for versioning
* Using ``requirements*.txt`` files now from setup.py
* Fixing dependency on cyordereddict to be for Python <3.6 instead of <3.5
* Jumping by samtools coordinate string now also allowed

0.8.0 (2016-10-31)
------------------

* Adding ``Header.has_header_line`` for querying existence of header line
* ``Header.add_*_line`` return a ``bool`` no indicating any conflicts
* Construction of Writer uses samples within header and no extra parameter (breaks API)

0.7.0 (2016-09-25)
------------------

* Smaller improvements and fixes to documentation
* Adding Codacy coverage and static code analysis results to README
* Various smaller code cleanup triggered by Codacy results
* Adding ``__eq__``, ``__neq__`` and ``__hash__`` to data types (where applicable)

0.6.0 (2016-09-25
-----------------

* Refining implementation for breakend and symbolic allele class
* Removing ``record.SV_CODES``
* Refactoring parser module a bit to make the code cleaner
* Fixing small typos and problems in documentation

0.5.0 (2016-09-24)
------------------

* Deactivating warnings on record parsing by default because of performance
* Adding validation for ``INFO`` and ``FORMAT`` fields on reading (#8)
* Adding predefined ``INFO`` and ``FORMAT`` fields to ``pyvcf.header`` (#32)

0.4.1 (2016-09-22)
------------------

* Initially enabling codeclimate

0.4.0 (2016-09-22)
------------------

* Exporting constants for encoding variant types
* Exporting genotype constants ``HOM_REF``, ``HOM_ALT``, ``HET``
* Implementing ``Call.is_phased``, ``Call.is_het``, ``Call.is_variant``, ``Call.is_phased``, ``Call.is_hom_ref``, ``Call.is_hom_alt``
* Removing ``Call.phased`` (breaks API, next release is 0.4.0)
* Adding tests, fixing bugs for methods of ``Call``

0.3.1 (2016-09-21)
------------------

* Work around ``FORMAT/FT`` being a string; this is done so in the Delly output

0.3.0 (2016-09-21)
------------------

* ``Reader`` and ``Writer`` can now be used as context manager (with ``with``)
* Including license in documentation, including Biopython license
* Adding support for writing bgzf files (taken from Biopython)
* Adding support for parsing arrays in header lines
* Removing ``example-4.1-bnd.vcf`` example file because v4.1 tumor derival lacks ``ID`` field
* Adding ``AltAlleleHeaderLine``, ``MetaHeaderLine``, ``PedigreeHeaderLine``, and ``SampleHeaderLine``
* Renaming ``SimpleHeaderFile`` to ``SimpleHeaderLine``
* Warn on missing ``FILTER`` entries on parsing
* Reordered parameters in ``from_stream`` and ``from_file`` (#18)
* Renamed ``from_file`` to ``from_stream`` (#18)
* Renamed ``Reader.jump_to`` to ``Reader.fetch``
* Adding ``header_without_lines`` function
* Generally extending API to make it esier to use
* Upgrading dependencies, enabling pyup-bot
* Greatly extending documentation

0.2.1 (2016-09-19)
------------------

* First release on PyPI
