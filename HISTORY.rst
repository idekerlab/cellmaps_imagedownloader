=======
History
=======

0.3.0 (2025-XX-XX)
-------------------

* Updated to work with 0.6 release of CM4AI dataset by
  supporting format in `manifest.csv`

* Added `r` prefix to regex strings with escape sequences to remove
  deprecation warnings that show up in python 3.11

* Updated to PEP 517 compliant build system

0.2.1 (2025-03-18)
------------------

* Add version bounds for required packages.

0.2.0 (2024-11-22)
------------------

* Added `--protein_list` and `--cell_line` flags that will fetch images from HPA for specific proteins and/or cell line.

* Bug fixes in fairscape registration.

0.1.1 (2024-08-26)
------------------

* Bug fix in adding gene node attributes. The bug was resulting in duplicate entries and
  missing gene names, related to ambiguous antibodies.

0.1.0 (2024-01-01)
------------------

* First release on PyPI.
