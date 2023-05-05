#! /usr/bin/env python

import os
from multiprocessing import Pool
import re
import csv
import subprocess
import logging
import logging.config
import requests
import time
from datetime import date
import json
from tqdm import tqdm
from cellmaps_utils import logutils
import cellmaps_imagedownloader
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError

logger = logging.getLogger(__name__)


def download_file_skip_existing(downloadtuple):
    """
    Downloads file in **downloadtuple** unless the file already exists
    with a size greater then 0 bytes, in which case function
    just returns

    :param downloadtuple: (download link, dest file path)
    :type downloadtuple: tuple
    :return: None upon success otherwise:
             (requests status code, text from request, downloadtuple)
    :rtype: tuple
    """
    if os.path.isfile(downloadtuple[1]) and os.path.getsize(downloadtuple[1]) > 0:
        return None
    return download_file(downloadtuple)


def download_file(downloadtuple):
    """
    Downloads file pointed to by 'download_url' to
    'destfile'

    :param downloadtuple: (download link, dest file path)
    :type downloadtuple: tuple
    :raises Exception: from requests library if there is an error or non 200 status
    :return: None upon success otherwise:
             (requests status code, text from request, downloadtuple)
    :rtype: tuple
    """
    logger.debug('Downloading ' + downloadtuple[0] + ' to ' + downloadtuple[1])
    with requests.get(downloadtuple[0], stream=True) as r:
        if r.status_code != 200:
            return r.status_code, r.text, downloadtuple
        with open(downloadtuple[1], 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
    return None


class ImageDownloader(object):
    """

    """
    def __init__(self):
        """

        """
        pass

    def download_images(self, download_list=None):
        """
        Subclasses should implement
        :param download_list: list of tuples where first element is
                              full URL of image to download and 2nd
                              element is destination path
        :type download_list: list
        :return: 
        """
        raise CellMapsImageDownloaderError('Subclasses should implement this')


class MultiProcessImageDownloader(ImageDownloader):
    """
    Uses multiprocess package to download images in parallel
    """

    def __init__(self, poolsize=1, skip_existing=False,
                 override_dfunc=None):
        """
        Constructor
        """
        super().__init__()
        self._poolsize = poolsize
        if override_dfunc is not None:
            self._dfunc = override_dfunc
        else:
            self._dfunc = download_file
            if skip_existing is True:
                self._dfunc = download_file_skip_existing

    def download_images(self, download_list=None):
        """
        Downloads images returning a list of failed downloads

        :param download_list:
        :return: of tuples (`http status code`, `text of error`, (`link`, `destfile`))
        :rtype: list
        """
        failed_downloads = []
        logger.debug('Poolsize for image downloader set to: ' +
                     str(self._poolsize))
        with Pool(processes=self._poolsize) as pool:
            num_to_download = len(download_list)
            logger.info(str(num_to_download) + ' images to download')
            t = tqdm(total=num_to_download, desc='Download',
                     unit='images')
            for i in pool.imap_unordered(self._dfunc,
                                         download_list):
                t.update()
                if i is not None:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug('Failed download: ' + str(i))
                    failed_downloads.append(i)
        return failed_downloads


class CellmapsImageDownloader(object):
    """
    Class to run algorithm
    """

    RED = 'red'
    """
    Red color directory name and color name
    in red color files
    """

    BLUE = 'blue'
    """
    Blue color directory name and color name
    in blue color files
    """
    GREEN = 'green'
    """
    Green color directory name and color name
    in green color files
    """

    YELLOW = 'yellow'
    """
    Yellow color directory name and color name
    in yellow color files
    """

    COLORS = [RED, BLUE, GREEN, YELLOW]
    """
    List of colors
    """

    SAMPLES_CSVFILE = 'samples.csv'
    """
    Copy of input csv file that is stored in output
    directory by the :py:meth:`~cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`
    """

    UNIQUE_CSVFILE = 'unique.csv'
    """
    Copy of input csv file that is stored in output
    directory by the :py:meth:`~cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`
    """

    SAMPLES_FILEKEY = 'samples'
    UNIQUE_FILEKEY = 'unique'

    IMAGE_GENE_NODE_ATTR_FILE = 'image_gene_node_attributes.tsv'
    IMAGE_GENE_NODE_ERRORS_FILE = 'image_gene_node_attributes.errors'
    IMAGE_GENE_NODE_COLS = ['name', 'represents', 'ambiguous',
                            'antibody', 'filename']

    def __init__(self, outdir=None,
                 imgsuffix='.jpg',
                 imagedownloader=MultiProcessImageDownloader(),
                 imagegen=None,
                 image_url=None,
                 skip_logging=False,
                 provenance=None,
                 input_data_dict=None):
        """
        Constructor

        :param outdir: directory where images will be downloaded to
        :type outdir: str
        :param imgsuffix: suffix to append to image file names
        :type imgsuffix: str
        :param imagedownloader: object that will perform image downloads
        :type imagedownloader: :py:class:`~cellmaps_downloader.runner.ImageDownloader`
        :param imagegen: gene node attribute generator for IF image data
        :type imagegen: :py:class:`~cellmaps_imagedownloader.gene.ImageGeneNodeAttributeGenerator`
        :param image_url: Base URL for image download
        :type image_url: str
        :param skip_logging:
        :type skip_logging: bool
        :param provenance:
        :type provenance: dict
        :param input_data_dict:
        :type input_data_dict: dict
        """
        if outdir is None:
            raise CellMapsImageDownloaderError('outdir is None')
        self._outdir = os.path.abspath(outdir)
        self._imagedownloader = imagedownloader
        self._imgsuffix = imgsuffix
        self._start_time = int(time.time())
        self._end_time = -1
        self._imagegen = imagegen
        self._image_url = image_url
        self._provenance = provenance
        self._input_data_dict = input_data_dict
        if skip_logging is None:
            self._skip_logging = False
        else:
            self._skip_logging = skip_logging
        self._samples_datasetid = None
        self._unique_datasetid = None
        self._softwareid = None
        self._image_gene_attrid = None

    @staticmethod
    def get_example_provenance(requiredonly=True,
                               with_ids=False):
        """
        Gets a dict of provenance parameters needed to add/register
        a dataset with FAIRSCAPE

        :param requiredonly: If ``True`` only output required fields,
                             otherwise output all fields. This value
                             is ignored if **with_ids** is ``True``
        :type requiredonly: bool
        :param with_ids: If ``True`` only output the fields
                         to set dataset guids and ignore value of
                         **requiredonly** parameter.
        :type with_ids: bool
        :return:
        """
        base_dict = {'organization-name': 'Name of organization',
                     'project-name': 'Name of project'}
        if with_ids is not None and with_ids is True:
            base_dict.update({'samples_file': {'guid': 'ID of dataset'},
                              'unique_file': {'guid': 'ID of dataset'}})
            return base_dict

        field_dict = {'name': 'Name of dataset',
                      'author': 'Author of dataset',
                      'version': 'Version of dataset',
                      'date-published': 'Date dataset was published',
                      'description': 'Description of dataset',
                      'data-format': 'Format of data'}

        if requiredonly is None or requiredonly is False:
            field_dict.update({'url': 'URL of datset',
                               'used-by': '?',
                               'derived-from': '?',
                               'associated-publication': '?',
                               'additional-documentation': '?'})

        base_dict.update({'samples_file': field_dict,
                          'unique_file': field_dict})
        return base_dict

    def _create_output_directory(self):
        """
        Creates output directory if it does not already exist

        :raises CellmapsDownloaderError: If output directory is None or if directory already exists
        """
        if os.path.isdir(self._outdir):
            raise CellMapsImageDownloaderError(self._outdir + ' already exists')

        for cur_color in CellmapsImageDownloader.COLORS:
            cdir = os.path.join(self._outdir, cur_color)
            if not os.path.isdir(cdir):
                logger.debug('Creating directory: ' + cdir)
                os.makedirs(cdir,
                            mode=0o755)
            else:
                logger.debug(cdir + ' already exists')

    def _run_cmd(self, cmd):
        """
        Runs hidef command as a command line process
        :param cmd_to_run: command to run as list
        :type cmd_to_run: list
        :return: (return code, standard out, standard error)
        :rtype: tuple
        """
        logger.debug('Running command: ' + str(cmd))
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        out, err = p.communicate()

        return p.returncode, out, err

    def _register_software(self):
        """
        Registers this tool

        :return:
        """
        cmd = ['fairscape-cli', 'rocrate', 'add', 'software',
               '--name', cellmaps_imagedownloader.__name__,
               '--description', cellmaps_imagedownloader.__description__,
               '--author', cellmaps_imagedownloader.__author__,
               '--version', cellmaps_imagedownloader.__version__,
               '--file-format', '.py',
               '--url', cellmaps_imagedownloader.__repo_url__,
               self._outdir]
        exit_code, out_str, err_str = self._run_cmd(cmd)

        logger.debug('add dataset exit code: ' + str(exit_code))
        if exit_code != 0:
            raise CellMapsImageDownloaderError('Error adding dataset: ' +
                                               str(out_str) + ' : ' + str(err_str))
        logger.debug('add data set out_str: ' + str(out_str))
        logger.debug('add data set err_str: ' + str(err_str))
        self._softwareid = out_str

    def _register_image_gene_node_attrs(self):
        """
        Registers image_gene_node_attributes.tsv file with create as a dataset

        """

        self._image_gene_attrid = self._add_dataset_to_crate(crate_path=self._outdir,
                                   source_file=self.get_image_gene_node_attributes_file(),
                                   data_dict={'name': cellmaps_imagedownloader.__name__ + ' output file',
                                              'description': 'Image gene node attributes file',
                                              'data-format': 'tsv',
                                              'author': '',
                                              'version': cellmaps_imagedownloader.__version__,
                                              'date-published': date.today().strftime('%m-%d-%Y')})

    def _add_dataset_to_crate(self, crate_path=None, data_dict=None,
                              source_file=None, skip_copy=True):
        """

        :param crate_path:
        :param data_dict:
        :return:
        """
        cmd = ['fairscape-cli', 'rocrate', 'add', 'dataset',
               '--name', data_dict['name'],
               '--version', data_dict['version'],
               '--data-format', data_dict['data-format'],
               '--description', data_dict['description'],
               '--date-published', data_dict['date-published'],
               '--author', data_dict['author'],
               '--source-filepath', source_file]
        if skip_copy is not None and skip_copy is False:
            cmd.append('--destination-filepath')
            cmd.append(os.path.join(crate_path,
                                    os.path.basename(source_file)))
        cmd.append(crate_path)
        exit_code, out_str, err_str = self._run_cmd(cmd)
        logger.debug('add dataset exit code: ' + str(exit_code))
        if exit_code != 0:
            raise CellMapsImageDownloaderError('Error adding dataset: ' +
                                               str(out_str) + ' : ' + str(err_str))
        logger.debug('add data set out_str: ' + str(out_str))
        logger.debug('add data set err_str: ' + str(err_str))
        return out_str

    def _register_computation(self):
        """

        :return:
        """
        cmd = ['fairscape-cli', 'rocrate', 'add', 'computation',
               '--name', cellmaps_imagedownloader.__name__ + ' computation',
               '--author', '',
               '--run-by', '',
               '--command', '',
               '--description', 'run of ' + cellmaps_imagedownloader.__name__,
               '--used-software', self._softwareid,
               '--used-dataset', self._unique_datasetid,
               '--used-dataset', self._samples_datasetid,
               '--generated', self._image_gene_attrid,
               self._outdir]
        exit_code, out_str, err_str = self._run_cmd(cmd)
        logger.debug('add dataset exit code: ' + str(exit_code))
        if exit_code != 0:
            raise CellMapsImageDownloaderError('Error adding dataset: ' +
                                               str(out_str) + ' : ' + str(err_str))
        logger.debug('add data set out_str: ' + str(out_str))
        logger.debug('add data set err_str: ' + str(err_str))

    def _create_run_crate(self):
        """

        :return:
        """
        # create rocrate
        try:
            cmd = ['fairscape-cli', 'rocrate', 'create',
                   '--name', 'cellmaps_imagedownloader',
                   '--organization-name', self._provenance['organization-name'],
                   '--project-name', self._provenance['project-name'],
                   os.path.abspath(self._outdir)]
        except TypeError as te:
            raise CellMapsImageDownloaderError('Invalid provenance: ' + str(te))
        except KeyError as ke:
            raise CellMapsImageDownloaderError('Provenance data malformed: ' + str(ke))
        exit_code, out_str, err_str = self._run_cmd(cmd)
        logger.debug('creation of crate stdout: ' + str(out_str))
        logger.debug('creation of crate stdout: ' + str(err_str))
        logger.debug('creation of crate exit code: ' + str(exit_code))
        if exit_code != 0:
            raise CellMapsImageDownloaderError('Error creating crate: ' +
                                               str(out_str) + ' : ' + str(err_str))

    def _register_input_datasets(self):
        """
        Registers samples and unique input datasets with FAIRSCAPE
        setting **self._samples_datasetid** and **self._unique_datasetid**
        values.

        :return:
        """

        if 'guid' in self._provenance[CellmapsImageDownloader.SAMPLES_FILEKEY]:
            self._samples_datasetid = self._provenance[CellmapsImageDownloader.SAMPLES_FILEKEY]['guid']
        if 'guid' in self._provenance[CellmapsImageDownloader.UNIQUE_FILEKEY]:
            self._unique_datasetid = self._provenance[CellmapsImageDownloader.UNIQUE_FILEKEY]['guid']

        if self._samples_datasetid is not None and self._unique_datasetid is not None:
            logger.debug('Both samples and unique have dataset ids. Just returning')
            return

        if self._samples_datasetid is None:
            # write file and add samples dataset
            self._samples_datasetid = self._add_dataset_to_crate(crate_path=self._outdir,
                                                                 data_dict=self._provenance[CellmapsImageDownloader.SAMPLES_FILEKEY],
                                                                 source_file=self._input_data_dict[CellmapsImageDownloader.SAMPLES_FILEKEY],
                                                                 skip_copy=False)
            logger.debug('Samples dataset id: ' + str(self._samples_datasetid))

        if self._unique_datasetid is None:
            # write file and add unique dataset
            self._unique_datasetid = self._add_dataset_to_crate(crate_path=self._outdir,
                                                                data_dict=self._provenance[CellmapsImageDownloader.UNIQUE_FILEKEY],
                                                                source_file=self._input_data_dict[CellmapsImageDownloader.UNIQUE_FILEKEY],
                                                                skip_copy=False)
            logger.debug('Unique dataset id: ' + str(self._unique_datasetid))

    def _get_color_download_map(self):
        """
        Creates a dict where key is color name and value is directory
        path for files for that color

        ``{'red': '/tmp/foo/red'}``

        :return: map of colors to directory paths
        :rtype: dict
        """
        color_d_map = {}
        for c in CellmapsImageDownloader.COLORS:
            color_d_map[c] = os.path.join(self._outdir, c)
        return color_d_map

    def _get_sample_url_and_filename(self, sample=None, color=None):
        """

        :param sample:
        :return:
        """
        file_name = sample['if_plate_id'] + '_' + sample['position'] + '_' + sample['sample'] + '_' + color + self._imgsuffix
        return self._image_url + '/' + re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + file_name, file_name

    def _get_download_tuples_from_csv(self):
        """
        Gets download list from CSV file for the 4 colors

        :return: list of (image download URL prefix,
                          file path where image should be written)
        :rtype: list
        """
        dtuples = []

        color_d_map = self._get_color_download_map()
        for row in self._imagegen.get_samples_list():
            for c in CellmapsImageDownloader.COLORS:
                image_url, file_name = self._get_sample_url_and_filename(sample=row, color=c)
                dtuples.append((image_url,
                                os.path.join(color_d_map[c], file_name)))
        return dtuples

    def _write_task_start_json(self):
        """
        Writes task_start.json file with information about
        what is to be run

        """
        data = {'image_downloader': str(self._imagedownloader),
                'image_suffix': self._imgsuffix}

        if self._input_data_dict is not None:
            data.update({'commandlineargs': self._input_data_dict})

        logutils.write_task_start_json(outdir=self._outdir,
                                       start_time=self._start_time,
                                       version=cellmaps_imagedownloader.__version__,
                                       data=data)

    def _retry_failed_images(self, failed_downloads=None):
        """

        :param failed_downloads:
        :return:
        """
        downloads_to_retry = []
        error_code_map = {}
        for entry in failed_downloads:
            if entry[0] not in error_code_map:
                error_code_map[entry[0]] = 0
            error_code_map[entry[0]] += 1
            downloads_to_retry.append(entry[2])
        logger.debug('Failed download counts by http error code: ' + str(error_code_map))
        return self._imagedownloader.download_images(downloads_to_retry)

    def _download_images(self, max_retry=5):
        """
        Uses downloader specified in constructor to download images noted in
        tsvfile file also specified in constructor
        :raises CellMapsImageDownloaderError: if image downloader is ``None`` or
                                         if there are failed downloads
        :return: 0 upon success otherwise, failure
        :rtype: int
        """
        if self._imagedownloader is None:
            raise CellMapsImageDownloaderError('Image downloader is None')

        downloadtuples = self._get_download_tuples_from_csv()

        failed_downloads = self._imagedownloader.download_images(downloadtuples)
        retry_count = 0
        while len(failed_downloads) > 0 and retry_count < max_retry:
            retry_count += 1
            logger.error(str(len(failed_downloads)) +
                         ' images failed to download. Retrying #' + str(retry_count))

            # try one more time with files that failed
            failed_downloads = self._retry_failed_images(failed_downloads=failed_downloads)

        if len(failed_downloads) > 0:
            raise CellMapsImageDownloaderError('Failed to download: ' +
                                               str(len(failed_downloads)) + ' images')
        return 0

    def get_image_gene_node_attributes_file(self):
        """
        Gets full path to image gene node attribute file under output directory
        created when invoking :py:meth:`~cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`

        :return: Path to file
        :rtype: str
        """
        return os.path.join(self._outdir,
                            CellmapsImageDownloader.IMAGE_GENE_NODE_ATTR_FILE)

    def get_image_gene_node_errors_file(self):
        """
        Gets full path to image gene node attribute errors file under output directory
        created when invoking :py:meth:`~cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`

        :return: Path to file
        :rtype: str
        """
        return os.path.join(self._outdir,
                            CellmapsImageDownloader.IMAGE_GENE_NODE_ERRORS_FILE)

    def _write_image_gene_node_attrs(self, gene_node_attrs=None,
                                     errors=None):
        """

        :param gene_node_attrs:
        :param errors:
        :return:
        """
        with open(self.get_image_gene_node_attributes_file(), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=CellmapsImageDownloader.IMAGE_GENE_NODE_COLS, delimiter='\t')
            writer.writeheader()
            for key in gene_node_attrs:
                writer.writerow(gene_node_attrs[key])
        if errors is not None:
            with open(self.get_image_gene_node_errors_file(), 'w') as f:
                for e in errors:
                    f.write(str(e) + '\n')

    def run(self):
        """
        Downloads images to output directory specified in constructor
        using tsvfile for list of images to download

        :raises CellMapsImageDownloaderError: If there is an error
        :return: 0 upon success, otherwise failure
        """
        try:
            exitcode = 99
            self._create_output_directory()
            if self._skip_logging is False:
                logutils.setup_filelogger(outdir=self._outdir,
                                          handlerprefix='cellmaps_imagedownloader')
                self._write_task_start_json()

            self._create_run_crate()
            self._register_input_datasets()

            # Todo: uncomment when fixed
            # register software fails due to this bug:
            # https://github.com/fairscape/fairscape-cli/issues/7
            # self._register_software()

            # write image attribute data
            image_gene_node_attrs, errors = self._imagegen.get_gene_node_attributes()

            # write image attribute data
            self._write_image_gene_node_attrs(image_gene_node_attrs, errors)

            # Todo: uncomment when fixed
            # register image_gene_node_attrs fails due to this bug:
            # https://github.com/fairscape/fairscape-cli/issues/6
            # self._register_image_gene_node_attrs()

            # Todo: Need to put this below download images once we know
            #       how to register all 100k+ images

            # Todo: uncomment when above work
            # Above registrations need to work for this to work
            # register computation
            # self._register_computation()

            exitcode = self._download_images()
            # todo need to validate downloaded image data

            return exitcode
        finally:
            self._end_time = int(time.time())
            if self._skip_logging is False:
                # write a task finish file
                logutils.write_task_finish_json(outdir=self._outdir,
                                                start_time=self._start_time,
                                                end_time=self._end_time,
                                                status=exitcode)
