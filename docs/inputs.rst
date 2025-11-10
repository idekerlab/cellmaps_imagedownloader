=======
Inputs
=======

The tool requires one of the following inputs: a CSV file containing a list of IF images to download,  a TXT/CSV file
with a list of proteins for which IF images will be downloaded, or a single path to a TSV file located in the CM4AI
RO-Crate directory. It also requires path to file containing provenance information about input files in JSON format.

Below is the list and description of each input accepted by the tool.

- ``samples.csv``:
    CSV file with list of IF images to download. The file must contain the columns listed below. The ``z`` column
    was introduced in version 0.3.0 to keep track of the z-stack identifier when copying data from CM4AI crates.
    If the column is omitted it will default to ``z01_``. A ``linkprefix`` column is optional and only needed when
    providing pre-downloaded CM4AI assets that should be copied instead of fetched from HPA.

    Definition of columns:

    * filename - Filename of image (string)
    * if_plate_id - ID of plate for acquired image (int)
    * position - Position in plate for acquired image (string)
    * sample - Sample number identifier for acquired image (int)
    * locations - Comma delimited list of manual annotations for image (string)
    * antibody - Name of antibody used for acquired image (string)
    * ensembl_ids - Comma delimited list of Ensembl IDs (string)
    * gene_names - Comma delimited list of genes (string)
    * z - Z slice identifier, including trailing underscore (for example ``z01_``)
    * linkprefix - *(optional)* Base URL or filesystem path used by CM4AI copy workflows

**Example:**

.. code-block::

    filename,if_plate_id,position,sample,status,locations,antibody,ensembl_ids,gene_names,z
    /archive/7/7_C5_1_,7,C5,1,35,"Cytosol,Nuclear speckles",HPA005910,ENSG00000011007,ELOA,z01_
    /archive/7/7_C5_2_,7,C5,2,35,"Cytosol,Nuclear speckles",HPA005910,ENSG00000011007,ELOA,z01_
    /archive/7/7_E8_1_,7,E8,1,35,Nuclear speckles,HPA006628,ENSG00000239306,RBM14,z01_
    /archive/7/7_E8_2_,7,E8,2,35,Nuclear speckles,HPA006628,ENSG00000239306,RBM14,z01_

- ``proteins.txt``:
    Path to a plain-text file containing one gene symbol per line. The entries will be matched against the ``name``
    field in ``proteinatlas.xml`` when the Human Protein Atlas workflow is triggered.

**Example:**

.. code-block::

    ELOA
    RBM14
    SRSF11
    MCM3
    APEX1


- ``CM4AI_TABLE_PATH``:
    Path to the CM4AI antibody table inside an RO-Crate directory. It is expected the directory also contains
    ``red/`` ``blue/`` ``green/`` ``yellow/`` directories with images. Two formats are currently supported:

    *Legacy TSV format (``*.tsv``)*

    * Antibody ID - describes the antibody ID for the antibody applied to stain the protein visible in the "green" channel. The antibody ID can be looked up at proteinatlas.org to find out more information about the antibody.
    * ENSEMBL ID - indicates the ENSEMBL ID(s) of the gene(s) of the proteins visualized in the "green" channel.
    * Treatment - refers to how the cells that are depicted in the image were treated (with Paclitaxel, Vorinostat, or untreated)
    * Well - refers to the well coordinate on the 96-well plate
    * Region - is a unique identifier for the position in the well, where the cells were acquired

**Example:**

.. code-block::

    Antibody ID	ENSEMBL ID	Treatment	Well	Region
    CAB079904	ENSG00000187555	untreated	C1	R1
    CAB079904	ENSG00000187555	untreated	C1	R2
    CAB079904	ENSG00000187555	untreated	C1	R3
    CAB079904	ENSG00000187555	untreated	C1	R5

    *Manifest CSV format (``manifest.csv`` introduced in CM4AI v0.6)*

    * HPA_Antibody_ID - antibody identifier used on proteinatlas.org
    * ENSEMBL ID - gene identifiers for the antibody target(s)
    * Baselink - basename that contains the plate/region metadata encoded in the filename (for example ``B2AI_1_vorinostat_B1_R2_z01_``)
    * Plate - numeric plate identifier
    * Treatment - treatment label (untreated, paclitaxel, vorinostat, …)
    * Well - well coordinate on the plate
    * Region - region identifier within a well

    The manifest is parsed to derive the ``filename``, ``z`` and ``linkprefix`` fields that downstream components need, so do not rename the
    original columns.

- ``provenance.json``:
    Path to file containing provenance information about input files in JSON format.
    This is required and not including will output error message with example of file.

**Example:**

.. code-block:: json

    {
      "name": "Example input dataset",
      "organization-name": "CM4AI",
      "project-name": "Example",
      "edgelist": {
        "name": "sample edgelist",
        "author": "Krogan Lab",
        "version": "1.0",
        "date-published": "07-31-2023",
        "description": "AP-MS Protein interactions on HSC2 cell line, example dataset",
        "data-format": "tsv"
      },
      "baitlist": {
        "name": "sample baitlist",
        "author": "Krogan Lab",
        "version": "1.0",
        "date-published": "07-31-2023",
        "description": "AP-MS Baits used for Protein interactions on HSC2 cell line",
        "data-format": "tsv"
      },
      "samples": {
        "name": "u2os HPA IF images",
        "author": "Author of dataset",
        "version": "Version of dataset",
        "date-published": "Date dataset was published",
        "description": "Description of dataset",
        "data-format": "csv"
      },
      "unique": {
        "name": "u2os HPA IF images unique",
        "author": "Author of dataset",
        "version": "Version of dataset",
        "date-published": "Date dataset was published",
        "description": "Description of dataset",
        "data-format": "csv"
      }
    }

Automatic HPA workflow
----------------------

If neither ``--samples`` nor ``--cm4ai_table`` is supplied, the command line tool automatically downloads the
``proteinatlas.xml.gz`` file, extracts it, and builds a ``samples.csv`` in the output directory before starting the
main run. The optional ``--protein_list`` flag should point to a newline-delimited text file when you want to restrict
that download to specific genes, and ``--cell_line`` (default ``U2OS``) controls which cell line’s images are kept.
To bypass this automatic download, always provide a ``--samples`` CSV or ``--cm4ai_table`` path explicitly.
