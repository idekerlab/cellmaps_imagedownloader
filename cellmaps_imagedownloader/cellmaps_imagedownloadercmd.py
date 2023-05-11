#! /usr/bin/env python

import argparse
import sys
import logging
import logging.config
import json
import warnings

from cellmaps_utils import logutils
from cellmaps_utils import constants
import cellmaps_imagedownloader
from cellmaps_imagedownloader.runner import MultiProcessImageDownloader
from cellmaps_imagedownloader.runner import FakeImageDownloader
from cellmaps_imagedownloader.runner import CellmapsImageDownloader
from cellmaps_imagedownloader.gene import ImageGeneNodeAttributeGenerator

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
    parser.add_argument('--samples',
                        help='CSV file with list of IF images to download '
                             'in format of filename,if_plate_id,position,'
                             'sample,status,locations,antibody,ensembl_ids,'
                             'gene_names\n/archive/1/1_A1_1_,1,A1,1,35,'
                             'Golgi apparatus,HPA000992,ENSG00000066455,GOLGA5')
    parser.add_argument('--unique',
                        help='CSV file of unique samples '
                             'in format of:\n'
                             'antibody,ensembl_ids,gene_names,atlas_name,'
                             'locations,n_location\n'
                             'HPA040086,ENSG00000094914,AAAS,U-2 OS,'
                             'Nuclear membrane,1')
    parser.add_argument('--provenance',
                        help='Path to file containing provenance '
                             'information about input files in JSON format. '
                             'This is required and not including will output '
                             'and error message with example of file')
    parser.add_argument('--image_url', default='https://images.proteinatlas.org',
                        help='Base URL for downloading IF images')
    parser.add_argument('--fake_images', action='store_true',
                        help='If set, 1st image of each color is downloaded '
                             'and subsequent images are just copies of those '
                             'images')
    parser.add_argument('--poolsize', type=int,
                        default=4,
                        help='If using multiprocessing image downloader, '
                             'this sets number of current downloads to run. '
                             'Note: Going above the default overloads the server')
    parser.add_argument('--imgsuffix', default='.jpg',
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
                        help='If set, output.log, error.log and '
                             'task_#_start/finish.json '
                             'files will not be created')
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module. Messages are '
                             'output at these python logging levels '
                             '-v = ERROR, -vv = WARNING, -vvv = INFO, '
                             '-vvvv = DEBUG, -vvvvv = NOTSET (default no '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 cellmaps_imagedownloader.__version__))

    return parser.parse_args(args)


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

filename,if_plate_id,position,sample,status,locations,antibody,ensembl_ids,gene_names
/archive/1/1_A1_1_,1,A1,1,35,Golgi apparatus,HPA000992,ENSG00000066455,GOLGA5

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

        imagegen = ImageGeneNodeAttributeGenerator(unique_list=ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile(theargs.unique),
                                                   samples_list=ImageGeneNodeAttributeGenerator.get_samples_from_csvfile(theargs.samples))

        if theargs.fake_images is True:
            warnings.warn('FAKE IMAGES ARE BEING DOWNLOADED!!!!!')
            dloader = FakeImageDownloader()
        else:
            dloader = MultiProcessImageDownloader(poolsize=theargs.poolsize,
                                                  skip_existing=theargs.skip_existing)
        return CellmapsImageDownloader(outdir=theargs.outdir,
                                       imagedownloader=dloader,
                                       imgsuffix=theargs.imgsuffix,
                                       imagegen=imagegen,
                                       image_url=theargs.image_url,
                                       skip_logging=theargs.skip_logging,
                                       input_data_dict=theargs.__dict__,
                                       provenance=json_prov,
                                       skip_failed=theargs.skip_failed).run()
    except Exception as e:
        logger.exception('Caught exception: ' + str(e))
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
