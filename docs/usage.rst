=====
Usage
=====

This page should provide information on how to use cellmaps_imagedownloader

In a project
--------------

To use cellmaps_imagedownloader in a project::

    import cellmaps_imagedownloader

On the command line
---------------------

For information invoke :code:`cellmaps_imagedownloadercmd.py -h`

**Example usage**

The cell maps image downloader requires the following input files: 

1) samples file: CSV file with list of IF images to download (see sample samples file in examples folder)
2) unique file: CSV file of unique samples (see sample unique file in examples folder)
3) provenance: file containing provenance information about input files in JSON format (see sample provenance file in examples folder)

.. code-block::

   cellmaps_imagedownloadercmd.py ./cellmaps_imagedownloader_outdir  --samples examples/samples.csv --unique examples/unique.csv --provenance examples/provenance.json

Via Docker
---------------

**Example usage**


.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` idekerlab/cellmaps_imagedownloader:0.1.0 cellmaps_imagedownloadercmd.py ./cellmaps_imagedownloader_outdir  --samples examples/samples.csv --unique examples/unique.csv --provenance examples/provenance.json


