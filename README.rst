=============================================
Cell Maps ImmunoFluorescent Image Downloader
=============================================


.. image:: https://img.shields.io/pypi/v/cellmaps_imagedownloader.svg
        :target: https://pypi.python.org/pypi/cellmaps_imagedownloader

.. image:: https://app.travis-ci.com/idekerlab/cellmaps_imagedownloader.svg?branch=main
    :target: https://app.travis-ci.com/idekerlab/cellmaps_imagedownloader

.. image:: https://readthedocs.org/projects/cellmaps-imagedownloader/badge/?version=latest
        :target: https://cellmaps-imagedownloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://zenodo.org/badge/635992073.svg
        :target: https://zenodo.org/doi/10.5281/zenodo.10607371
        :alt: Zenodo DOI badge


Downloads ImmunoFluorescent image data from `Human Protein Atlas <https://www.proteinatlas.org/>`__


* Free software: MIT license
* Documentation: https://cellmaps-imagedownloader.readthedocs.io.


Dependencies
------------

* `cellmaps_utils <https://pypi.org/project/cellmaps-utils>`__
* `requests <https://pypi.org/project/requests>`__
* `mygene <https://pypi.org/project/mygene>`__
* `tqdm <https://pypi.org/project/tqdm>`__

Compatibility
-------------

* Python 3.8+

Installation
------------

.. code-block::

    pip install cellmaps_imagedownloader

**Or directly from source:**

.. code-block::

   git clone https://github.com/idekerlab/cellmaps_imagedownloader
   cd cellmaps_imagedownloader
   make dist
   pip install dist/cellmaps_imagedownloader*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image 

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub




Needed files
------------

* samples file: CSV file with list of IF images to download (see sample samples file in examples folder)
* unique file: CSV file of unique samples (see sample unique file in examples folder)
* provenance: file containing provenance information about input files in JSON format (see sample provenance file in examples folder)

Usage
-----

For information invoke :code:`cellmaps_imagedownloadercmd.py -h`

**Example usage**


.. code-block::

    cellmaps_imagedownloadercmd.py ./cellmaps_imagedownloader_outdir  --samples examples/samples.csv --unique examples/unique.csv --provenance examples/provenance.json


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**


.. code-block::

   Coming soon...

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
