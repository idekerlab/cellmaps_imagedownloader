=====
Usage
=====

This script facilitates the downloading of ImmunoFluorescent (IF) labeled images from the `Human Protein Atlas`_ (HPA).
The tool requires an output directory to write results to and either a TSV_ file in CM4AI_ RO-Crate_ format,
or CSV_ file with list of IF images to download and CSV_ file of unique samples.

In a project
--------------

To use cellmaps_imagedownloader in a project::

    import cellmaps_imagedownloader

On the command line
---------------------

For information invoke :code:`cellmaps_imagedownloadercmd.py -h`

**Usage**

.. code-block::

  cellmaps_imagedownloadercmd.py OUTPUT_DIRECTORY [--provenance PROVENANCE_PATH] [OPTIONS]

**Arguments**

- ``outdir``
    The directory where the output will be written to.

*Required*

- ``--provenance PROVENANCE_PATH``
    Path to file containing provenance information about input files in JSON format.

*Optional but either `samples` and `unique` parameters, or `cm4ai_table` parameter is required*

- ``--samples SAMPLES_PATH``
    CSV file with list of IF images to download. The file follow a specific format with columns such as
    filename, if_plate_id, position, sample, locations, antibody, ensembl_ids, and gene_names.

- ``--unique UNIQUE_PATH``
    CSV file of unique samples. The file should have columns like antibody, ensembl_ids, gene_names, atlas_name, locations, and n_location.

- ``--cm4ai_table CM4AI_TABLE_PATH``
    Path to TSV file in CM4AI RO-Crate directory.

*Optional*

- ``--proteinatlasxml``: URL or path to ``proteinatlas.xml`` or ``proteinatlas.xml.gz`` file.
- ``--fake_images``: If set, the first image of each color is downloaded, and subsequent images are copies of those images. If ``--cm4ai_table`` flag is set, the ``--fake_images`` flag is ignored.
- ``--poolsize``: If using multiprocessing image downloader, this sets the number of current downloads to run.
- ``--imgsuffix``: Suffix for images to download (default is ``.jpg``).
- ``--skip_existing``: If set, skips download if the image already exists and has a size greater than 0 bytes.
- ``--skip_failed``: If set, ignores images that failed to download after retries.
- ``--logconf``: Path to the python logging configuration file.
- ``--skip_logging``: If set, certain log files will not be created.
- ``--verbose``, ``-v``: Increases verbosity of logger to standard error for log messages.
- ``--version``: Shows the current version of the tool.

**Example usage**

The example file can be downloaded from `cm4ai.org <https://cm4ai.org>`__. Go to **Products -> Data**, log in, and download file for IF images with the desired treatment,
then unpack the tar.gz (``tar -xzvf filename.tar.gz``).

.. code-block::

   cellmaps_imagedownloadercmd.py ./cellmaps_imagedownloader_outdir  --cm4ai_table path/to/downloaded/unpacked/dir --provenance examples/provenance.json


Alternatively, use the files in the example directory in the repository:

1) samples file: CSV_ file with list of IF images to download (see sample samples file in examples folder)
2) unique file: CSV_ file of unique samples (see sample unique file in examples folder)
3) provenance: file containing provenance information about input files in JSON format (see sample provenance file in examples folder)

.. code-block::

   cellmaps_imagedownloadercmd.py ./cellmaps_imagedownloader_outdir  --samples examples/samples.csv --unique examples/unique.csv --provenance examples/provenance.json

Via Docker
---------------

**Example usage**


.. code-block::

   Coming soon...

.. _RO-Crate: https://www.researchobject.org/ro-crate
.. _CSV: https://en.wikipedia.org/wiki/Comma-separated_values
.. _TSV: https://en.wikipedia.org/wiki/Tab-separated_values
.. _Human Protein Atlas: https://www.proteinatlas.org
.. _CM4AI: https://cm4ai.org


