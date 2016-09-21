=======
History
=======

HEAD
----

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
