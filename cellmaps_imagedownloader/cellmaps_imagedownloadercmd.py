#! /usr/bin/env python

import argparse
import glob
import sys
import logging
import logging.config
import json
import warnings
import os
import zipfile

import requests
from cellmaps_utils import logutils
from cellmaps_utils import constants
import cellmaps_imagedownloader
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError
from cellmaps_imagedownloader.runner import MultiProcessImageDownloader
from cellmaps_imagedownloader.runner import FakeImageDownloader
from cellmaps_imagedownloader.runner import CellmapsImageDownloader
from cellmaps_imagedownloader.runner import CM4AICopyDownloader
from cellmaps_imagedownloader.gene import ImageGeneNodeAttributeGenerator
from cellmaps_imagedownloader.gene import CM4AITableConverter
from cellmaps_imagedownloader.proteinatlas import ProteinAtlasReader, ProteinAtlasProcessor
from cellmaps_imagedownloader.proteinatlas import ProteinAtlasImageUrlReader
from cellmaps_imagedownloader.proteinatlas import ImageDownloadTupleGenerator
from cellmaps_imagedownloader.proteinatlas import LinkPrefixImageDownloadTupleGenerator
from cellmaps_imagedownloader.proteinatlas import CM4AIImageCopyTupleGenerator

logger = logging.getLogger(__name__)


def _parse_arguments(desc, args):
    """
    Parses command line arguments

    :param desc: description to display on command line
    :type desc: str
    :param args: command line arguments usually :py:func:`sys.argv[1:]`
    :type args: list
    :return: arguments parsed by :py:mod:`argparse`
    :rtype: :py:class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=constants.ArgParseFormatter)
    parser.add_argument('outdir',
                        help='Directory to write results to')
    parser.add_argument('--cm4ai_table',
                        help='Path to TSV file in CM4AI RO-Crate directory. '
                             ' It is expected the directory also contains '
                             'red/ blue/ green/ yellow/ directories with images. '
                             'The TSV file is expected to have the following '
                             'columns: Antibody ID     ENSEMBL ID      '
                             'Treatment      '
                             ' Well    Region')
    parser.add_argument('--dataverse_doi',
                        help='DOI to datasets within Dataverse')
    parser.add_argument('--dataverse_dataset', help='Specifies the name of dataset with images that will '
                                                    'be downloaded (required if you specify dataverse DOI) e.g. '
                                                    'cm4ai_chromatin_mda-mb-468_paclitaxel_ifimage_0.1_alpha.zip.')
    parser.add_argument('--samples',
                        help='CSV file with list of IF images to download '
                             'in format of filename,if_plate_id,position,'
                             'sample,status,locations,antibody,ensembl_ids,'
                             'gene_names\n/archive/1/1_A1_1_,1,A1,1,35,'
                             'Golgi apparatus,HPA000992,ENSG00000066455,GOLGA5')
    parser.add_argument('--unique',
                        help='(Deprecated: Using --samples flag only is enough) '
                             'CSV file of unique samples '
                             'in format of:\n'
                             'antibody,ensembl_ids,gene_names,atlas_name,'
                             'locations,n_location\n'
                             'HPA040086,ENSG00000094914,AAAS,U-2 OS,'
                             'Nuclear membrane,1')
    parser.add_argument('--protein_list', help='List of proteins for which HPA images will be downloaded')
    parser.add_argument('--cell_line',
                        help='Cell line for which HPA images will be downloaded. See available cell lines at '
                             'https://www.proteinatlas.org/humanproteome/cell+line.', default='U2OS')
    parser.add_argument('--provenance',
                        help='Path to file containing provenance '
                             'information about input files in JSON format. '
                             'This is required and not including will output '
                             'and error message with example of file')
    parser.add_argument('--proteinatlasxml',
                        default=ProteinAtlasReader.DEFAULT_PROTEINATLAS_URL,
                        help='URL or path to proteinatlas.xml or proteinatlas.xml.gz file '
                             'used to look for images not found in the standard location '
                             'on HPA')
    parser.add_argument('--fake_images', action='store_true',
                        help='If set, 1st image of each color is downloaded '
                             'and subsequent images are just copies of those '
                             'images')
    parser.add_argument('--poolsize', type=int,
                        default=MultiProcessImageDownloader.POOL_SIZE,
                        help='If using multiprocessing image downloader, '
                             'this sets number of current downloads to run. '
                             'Note: Going above the default overloads the server')
    parser.add_argument('--imgsuffix', default=CellmapsImageDownloader.IMG_SUFFIX,
                        help='Suffix for images to download')
    parser.add_argument('--skip_existing', action='store_true',
                        help='If set, skips download if image already exists and '
                             'has size greater then 0 bytes')
    parser.add_argument('--skip_failed', action='store_true',
                        help='If set, ignores images that failed to download after retries')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/'
                             'logging.config.html#logging-config-fileformat '
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')
    parser.add_argument('--skip_logging', action='store_true',
                        help='If set, output.log, error.log '
                             'files will not be created')
    parser.add_argument('--verbose', '-v', action='count', default=1,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module. Messages are '
                             'output at these python logging levels '
                             '-v = WARNING, -vv = INFO, '
                             '-vvv = DEBUG, -vvvv = NOTSET (default ERROR '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 cellmaps_imagedownloader.__version__))

    return parser.parse_args(args)


def _get_dataset_id(dataverse_api, dataverse_dataset):
    response = requests.get(dataverse_api)
    if response.status_code != 200:
        raise CellMapsImageDownloaderError(f"Failed to fetch dataverse JSON. Status code: {response.status_code}")

    dataset_json = response.json()

    data_id = None

    for file_entry in dataset_json["datasetVersion"]["files"]:
        if file_entry["label"] == dataverse_dataset:
            data_id = file_entry["dataFile"]["id"]
            break

    if data_id is None:
        raise CellMapsImageDownloaderError(f"Dataset '{dataverse_dataset}' not found in the dataverse JSON response.")

    return data_id


def download_zip_file(data_id, outdir, dataverse_dataset):
    download_url = f"https://dataverse.lib.virginia.edu/api/access/datafile/{data_id}"
    zip_file_path = os.path.join(outdir, dataverse_dataset)

    response = requests.get(download_url, stream=True)
    if response.status_code == 200:
        with open(zip_file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    else:
        raise CellMapsImageDownloaderError(f"Failed to download file. Status code: {response.status_code}")
    return zip_file_path


def main(args):
    """
    Main entry point for program

    :param args: arguments passed to command line usually :py:func:`sys.argv[1:]`
    :type args: list

    :return: return value of :py:meth:`cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`
             or ``2`` if an exception is raised
    :rtype: int
    """
    withguids_json = json.dumps(CellmapsImageDownloader.get_example_provenance(with_ids=True), indent=2)
    register_json = json.dumps(CellmapsImageDownloader.get_example_provenance(), indent=2)

    desc = """
Version {version}

Downloads immunofluorescent labeled images from the Human Protein Atlas
(https://www.proteinatlas.org/)

To use pass in a CSV file containing links to the images to download
from HPA via --samples flag

Format of CSV file:

filename,if_plate_id,position,sample,locations,antibody,ensembl_ids,gene_names
/archive/1/1_A1_1_,1,A1,1,Golgi apparatus,HPA000992,ENSG00000066455,GOLGA5

Definition of columns:

* filename - Filename of image (string)
* if_plate_id - ID of plate for acquired image (int)
* position - Position in plate for acquired image (string)
* sample - Sample number identifier for acquired image (int)
* locations - Comma delimited list of manual annotations for image (string)
* antibody - Name of antibody used for acquired image (string)
* ensembl_ids - Comma delimited list of Ensembl IDs (string)
* gene_names - Comma delimited list of genes (string)

The downloaded images are stored under the output directory
specified on the command line in color specific directories
(red, blue, green, yellow) with name format of:
<NAME>_color.jpg

Example: 1_A1_1_blue.jpg

The --unique flag should be given a CSF file containing best or desired antibodies to
use.

Format of CSV file:

antibody,ensembl_ids,gene_names,atlas_name,locations,n_location
HPA040086,ENSG00000094914,AAAS,U-2 OS,Nuclear membrane,1

Definition of columns:

* antibody - Name of antibody used for acquired image (string)
* ensembl_ids - Comma delimited list of ensembl IDs for target protein(s) (string)
* gene_names - Comma delimited list of gene names for target proteins(s) (string)
* atlas_name - Cell line (string)
* locations - Comma delimited list of subcellular locations (string)
* n_location - Number of subcellular locations (int)

In addition, the --provenance flag is required and must be set to a path
to a JSON file.

If datasets are already registered with FAIRSCAPE then the following is sufficient:

{withguids}

If datasets are NOT registered, then the following is required:

{register}

Additional optional fields for registering datasets include
'url', 'used-by', 'associated-publication', and 'additional-documentation'


    """.format(version=cellmaps_imagedownloader.__version__,
               withguids=withguids_json,
               register=register_json)
    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = cellmaps_imagedownloader.__version__

    try:
        logutils.setup_cmd_logging(theargs)
        if theargs.provenance is None:
            sys.stderr.write('\n\n--provenance flag is required to run this tool. '
                             'Please pass '
                             'a path to a JSON file with the following data:\n\n')
            sys.stderr.write('If datasets are already registered with '
                             'FAIRSCAPE then the following is sufficient:\n\n')
            sys.stderr.write(withguids_json + '\n\n')
            sys.stderr.write('If datasets are NOT registered, then the following is required:\n\n')
            sys.stderr.write(register_json + '\n\n')
            return 1

        # load the provenance as a dict
        with open(theargs.provenance, 'r') as f:
            json_prov = json.load(f)

        created_outdir = False
        if all(arg is None for arg in [theargs.cm4ai_table, theargs.samples, theargs.protein_list, theargs.cell_line,
                                       theargs.dataverse_doi]):
            raise CellMapsImageDownloaderError("Either protein list, cell line, samples, cm4ai table, or "
                                               "dataverse_doi parameter should be specified.")
        if theargs.cm4ai_table is None and theargs.dataverse_doi:
            if theargs.dataverse_dataset is None:
                raise CellMapsImageDownloaderError("Name of dataset within the DOI must be specified.")
            os.makedirs(theargs.outdir, mode=0o755)
            created_outdir = True
            doi_number = theargs.dataverse_doi.split(".org/")[
                1] if '.org/' in theargs.dataverse_doi else theargs.dataverse_doi
            dataverse_api = (f'https://dataverse.lib.virginia.edu/api/datasets/export?exporter=dataverse_json'
                             f'&persistentId=doi%3A{doi_number}')

            data_id = _get_dataset_id(dataverse_api, theargs.dataverse_dataset)
            zip_file_path = download_zip_file(data_id, theargs.outdir, theargs.dataverse_dataset)

            with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
                zip_ref.extractall(theargs.outdir)

            dataset_dir = os.path.join(theargs.outdir, theargs.dataverse_dataset.replace(".zip", ""))
            tsv_files = glob.glob(os.path.join(dataset_dir, "*.tsv"))
            if len(tsv_files) == 1:
                theargs.cm4ai_table = os.path.join(dataset_dir, os.path.basename(tsv_files[0]))
            else:
                raise FileNotFoundError("No .tsv file found or multiple .tsv files exist in the directory.")

        elif theargs.samples is None and (theargs.protein_list is not None or theargs.cell_line is not None):
            hpa_processor = ProteinAtlasProcessor(theargs.outdir, theargs.proteinatlasxml, theargs.protein_list,
                                                  theargs.cell_line)
            theargs.samples, theargs.proteinatlasxml = hpa_processor.get_sample_list_from_hpa()
            created_outdir = True

        if theargs.cm4ai_table is not None:
            converter = CM4AITableConverter(cm4ai=theargs.cm4ai_table)
            samples_list, unique_list = converter.get_samples_and_unique_lists()
            dloader = CM4AICopyDownloader()
        else:
            samples_list = ImageGeneNodeAttributeGenerator.get_samples_from_csvfile(theargs.samples)
            unique_list = None
            if theargs.unique is not None:
                logger.debug('Found --unique file passed in. loading')
                unique_list = ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile(theargs.unique)

            if theargs.fake_images is True:
                warnings.warn('FAKE IMAGES ARE BEING DOWNLOADED!!!!!')
                dloader = FakeImageDownloader()
            else:
                dloader = MultiProcessImageDownloader(poolsize=theargs.poolsize,
                                                      skip_existing=theargs.skip_existing)

        imagegen = ImageGeneNodeAttributeGenerator(unique_list=unique_list,
                                                   samples_list=samples_list)

        if theargs.cm4ai_table is not None:
            imageurlgen = CM4AIImageCopyTupleGenerator(samples_list=imagegen.get_samples_list())
        elif 'linkprefix' in imagegen.get_samples_list()[0]:
            imageurlgen = LinkPrefixImageDownloadTupleGenerator(samples_list=imagegen.get_samples_list())
        else:
            proteinatlas_reader = ProteinAtlasReader(theargs.outdir, proteinatlas=theargs.proteinatlasxml)
            proteinatlas_urlreader = ProteinAtlasImageUrlReader(reader=proteinatlas_reader)
            imageurlgen = ImageDownloadTupleGenerator(reader=proteinatlas_urlreader,
                                                      samples_list=imagegen.get_samples_list(),
                                                      valid_image_ids=imagegen.get_samples_list_image_ids())
        return CellmapsImageDownloader(outdir=theargs.outdir,
                                       imagedownloader=dloader,
                                       imgsuffix=theargs.imgsuffix,
                                       imagegen=imagegen,
                                       imageurlgen=imageurlgen,
                                       skip_logging=theargs.skip_logging,
                                       input_data_dict=theargs.__dict__,
                                       provenance=json_prov,
                                       skip_failed=theargs.skip_failed,
                                       existing_outdir=created_outdir).run()
    except Exception as e:
        logger.exception('Caught exception: ' + str(e))
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
