=======
Inputs
=======

The tool requires one of the following sets of inputs: either a CSV file containing a list of IF images to download and
a separate CSV file with a list of unique samples, or a single path to a TSV file located in the CM4AI RO-Crate directory.
It also requires path to file containing provenance information about input files in JSON format.

Below is the list and description of each input accepted by the tool.

- ``samples.csv``:
    CSV file with list of IF images to download. The file follow a specific format with columns such as
    filename, if_plate_id, position, sample, locations, antibody, ensembl_ids, and gene_names.

    Definition of columns:

    * filename - Filename of image (string)
    * if_plate_id - ID of plate for acquired image (int)
    * position - Position in plate for acquired image (string)
    * sample - Sample number identifier for acquired image (int)
    * locations - Comma delimited list of manual annotations for image (string)
    * antibody - Name of antibody used for acquired image (string)
    * ensembl_ids - Comma delimited list of Ensembl IDs (string)
    * gene_names - Comma delimited list of genes (string)

**Example:**

.. code-block::

    filename,if_plate_id,position,sample,status,locations,antibody,ensembl_ids,gene_names
    /archive/7/7_C5_1_,7,C5,1,35,"Cytosol,Nuclear speckles",HPA005910,ENSG00000011007,ELOA
    /archive/7/7_C5_2_,7,C5,2,35,"Cytosol,Nuclear speckles",HPA005910,ENSG00000011007,ELOA
    /archive/7/7_E8_1_,7,E8,1,35,Nuclear speckles,HPA006628,ENSG00000239306,RBM14
    /archive/7/7_E8_2_,7,E8,2,35,Nuclear speckles,HPA006628,ENSG00000239306,RBM14

- ``unique.csv``:
    CSV file of unique samples. The file should have columns like antibody, ensembl_ids, gene_names, atlas_name, locations, and n_location.

    Definition of columns:

    * antibody - Name of antibody used for acquired image (string)
    * ensembl_ids - Comma delimited list of ensembl IDs for target protein(s) (string)
    * gene_names - Comma delimited list of gene names for target proteins(s) (string)
    * atlas_name - Cell line (string)
    * locations - Comma delimited list of subcellular locations (string)
    * n_location - Number of subcellular locations (int)

**Example:**

.. code-block::

    antibody,ensembl_ids,gene_names,atlas_name,locations,n_location
    CAB004303,ENSG00000072110,ACTN1,U-2 OS,"Focal adhesion sites,Plasma membrane",2
    HPA046964,ENSG00000182718,ANXA2,U-2 OS,"Cytosol,Plasma membrane",2
    HPA002564,ENSG00000100823,APEX1,U-2 OS,Nucleoplasm,1
    HPA006550,ENSG00000111229,ARPC3,U-2 OS,Nucleoplasm,1

- ``CM4AI_TABLE_PATH``:
    Path to TSV file in CM4AI RO-Crate directory. It is expected the directory also contains ``red/`` ``blue/`` ``green/`` ``yellow/``
    directories with images.

    The .tsv file describes each image in the data set. Each row represents one image. The columns describe the
    staining from which the image was taken. The TSV file is expected to have the following columns:

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

