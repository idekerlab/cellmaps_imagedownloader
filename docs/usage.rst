Usage
=====

This script facilitates the downloading of ImmunoFluorescent (IF) labeled images from the `Human Protein Atlas`_ (HPA).
The tool requires an output directory to write results to and either a TSV_ or CSV_ file in CM4AI_ RO-Crate_ format,
or CSV_ file with list of IF images to download and CSV_ file of unique samples.

In a project
*************

To use cellmaps_imagedownloader in a project::

    import cellmaps_imagedownloader

On the command line
*********************

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

*Optional but either `samples`, `cm4ai_table`, `protein_list` or `cell_line` parameter is required*

- ``--samples SAMPLES_PATH``
    CSV file with list of IF images to download. The file follow a specific format with columns such as
    filename, if_plate_id, position, sample, locations, antibody, ensembl_ids, and gene_names.

- ``--protein_list``
    List of proteins for which HPA images will be downloaded. Each protein in new line.

- ``--cell_line``
    Cell line for which HPA images will be downloaded. See available cell lines at https://www.proteinatlas.org/humanproteome/cell+line.

- ``--cm4ai_table CM4AI_TABLE_PATH``
    Path to TSV or CSV file in CM4AI RO-Crate directory.

*Optional*

- ``--unique UNIQUE_PATH``: (Deprecated: Using --samples flag only is enough) CSV file of unique samples. The file should have columns like antibody, ensembl_ids, gene_names, atlas_name, locations, and n_location.
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


Example usage
--------------

Alternatively, use the files in the example directory in the repository:

1) samples file: CSV_ file with list of IF images to download (see sample samples file in examples folder)
2) unique file: CSV_ file of unique samples (see sample unique file in examples folder)
3) provenance: file containing provenance information about input files in JSON format (see sample provenance file in examples folder)

.. code-block::

   cellmaps_imagedownloadercmd.py ./cellmaps_imagedownloader_outdir  --samples examples/samples.csv --unique examples/unique.csv --provenance examples/provenance.json

Example usage using CM4AI 0.5 Alpha Data Release
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Visit `cm4ai.org <https://cm4ai.org>`__ and go to **Products -> Data Releases**

#. Scroll down and click **CM4AI 0.5 Alpha Data Release** link (circled in red)

    .. image:: images/datarelease_0.5link.png
        :alt: Link to CM4AI 0.5 data release circled in red

#. On the newly opened page/tab, scroll down to the **cm4ai_chromatin_mda-mb-468_paclitaxel_ifimage_0.1_alpha** entry
   and click the download icon (circled in red) to bring up a pop up dialog. Click **Zip Archive** (red arrow)
   to accept the usage agreement and download the dataset

    .. image:: images/0.5imagedownload_paclitaxel.png
        :alt: CM4AI 0.5 paclitaxel image zip download link circled in red

    .. note::

        For **vorinostat** dataset, look for **cm4ai_chromatin_mda-mb-468_vorinostat_ifimage_0.1_alpha.zip** entry and perform the same
        operations above.

    .. note::

        For **untreated** dataset, this tool was already run and its results can be found in **1.cm4ai_chromatin_mda-mb-468_untreated_imageloader_initialrun0.1alpha.zip** entry

#. Unzip file

    This can be done by double clicking on the file or if on a Mac/Linux machine by running the following
    on a command line:

    .. code-block::

        unzip cm4ai_chromatin_mda-mb-468_paclitaxel_ifimage_0.1_alpha.zip


#. Running cellmaps_imagedownloader command

    .. code-block::

        # Be sure to unzip the zip file above before running this step
        cellmaps_imagedownloadercmd.py ./paclitaxel_image  \
            --cm4ai_table cm4ai_chromatin_mda-mb-468_paclitaxel_ifimage_0.1_alpha/MDA-MB-468_paclitaxel_antibody_gene_table.tsv  \
            --provenance examples/provenance.json


Example usage February 2025 Data Release (Beta)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Visit `cm4ai.org <https://cm4ai.org>`__ and go to **Products -> Data Releases**

#. Scroll down and click **February 2025 Data Release (Beta)** link (circled in red)

    .. image:: images/datarelease_0.6link.png
        :alt: Link to CM4AI February 2025 data release circled in red

#. On the newly opened page/tab, scroll down to the **cm4ai-v0.6-beta-if-images-paclitaxel.zip** entry
   and click the download icon (circled in red) to bring up a pop up dialog. Click **Zip Archive** (red arrow) to
   accept the usage agreement and download the dataset

    .. image:: images/0.6imagedownload_paclitaxel.png
        :alt: CM4AI February 2025 data release paclitaxel circled in red

    .. note::

        For **vorinostat** dataset, look for **cm4ai-v0.6-beta-if-images-vorinostat.zip** entry and perform the same
        operations above. Same goes for untreated, look for **cm4ai-v0.6-beta-if-images-untreated.zip**

#. Unzip file

    This can be done by double clicking on the file or if on a Mac/Linux machine by running the following
    on a command line:

    .. code-block::

        unzip cm4ai-v0.6-beta-if-images-paclitaxel.zip


#. Running cellmaps_imagedownloader command

    .. code-block::

        # Be sure to unzip the zip file above before running this step
        cellmaps_imagedownloadercmd.py ./paclitaxel_image  \
            --cm4ai_table paclitaxel/manifest.csv  \
            --provenance examples/provenance.json

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


