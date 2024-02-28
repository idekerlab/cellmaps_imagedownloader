#! /usr/bin/env python

import os
from multiprocessing import Pool
import re
import csv
import shutil
import logging
import logging.config
import requests
import time
from datetime import date
import warnings
from tqdm import tqdm
from cellmaps_utils import logutils
from cellmaps_utils.provenance import ProvenanceUtil
from cellmaps_utils import constants
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

    .. note::

        Default download function used by :py:class:`~MultiProcessImageDownloader`

    :param downloadtuple: `(download link, dest file path)`
    :type downloadtuple: tuple
    :return: None upon success otherwise:
             `(requests status code, text from request, downloadtuple)`
    :rtype: tuple
    """
    logger.debug('Downloading ' + downloadtuple[0] + ' to ' + downloadtuple[1])
    try:
        with requests.get(downloadtuple[0], stream=True) as r:
            if r.status_code != 200:
                return r.status_code, r.text, downloadtuple
            with open(downloadtuple[1], 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:  # filter out keep-alive new chunks
                        f.write(chunk)
        return None
    except requests.exceptions.HTTPError as e:
        return -1, str(e), downloadtuple
    except requests.exceptions.ConnectionError as e:
        return -2, str(e), downloadtuple
    except requests.exceptions.Timeout as e:
        return -3, str(e), downloadtuple
    except requests.exceptions.RequestException as e:
        return -4, str(e), downloadtuple
    except Exception as e:
        return -5, str(e), downloadtuple


class ImageDownloader(object):
    """
    Abstract class that defines interface for classes that download images

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


class CM4AICopyDownloader(ImageDownloader):
    """
    Copies over images from CM4AI RO-Crate
    """

    def __init__(self):
        """
        Constructor
        """
        super().__init__()

    def download_images(self, download_list=None):
        """
        Subclasses should implement

        :param download_list: list of tuples where first element is
                              full URL of image to download and 2nd
                              element is destination path
        :type download_list: list
        :return:
        """
        num_to_copy = len(download_list)
        logger.info(str(num_to_copy) + ' images to copy')
        t = tqdm(total=num_to_copy, desc='Copy',
                 unit='images')
        for entry in download_list:
            t.update()
            shutil.copy(entry[0], entry[1])

        return []


class FakeImageDownloader(ImageDownloader):
    """
    Creates fake download by downloading
    the first image in each color from
    `Human Protein Atlas <https://www.proteinatlas.org/>`__
    and making renamed copies. The :py:func:`download_file` function
    is used to download the first image of each color

    """

    def __init__(self):
        """
        Constructor

        """
        super().__init__()
        warnings.warn('This downloader generates FAKE images\n'
                      'You have been warned!!!\n'
                      'Have a nice day')

    def download_images(self, download_list=None):
        """
        Downloads 1st image from server and then
        and makes renamed copies for subsequent images

        :param download_list:
        :type download_list: list of tuple
        :return:
        """
        num_to_download = len(download_list)
        logger.info(str(num_to_download) + ' images to download')
        t = tqdm(total=num_to_download, desc='Download',
                 unit='images')

        src_image_dict = {}
        # assume 1st four images are the colors for the first image
        for entry in download_list[0:4]:
            t.update()
            if download_file(entry) is not None:
                raise CellMapsImageDownloaderError('Unable to download ' +
                                                   str(entry))
            fname = os.path.basename(entry[1])
            color = re.sub('\..*$', '', re.sub('^.*_', '', fname))
            src_image_dict[color] = entry[1]

        for entry in download_list[5:]:
            t.update()
            fname = os.path.basename(entry[1])
            color = re.sub('\..*$', '', re.sub('^.*_', '', fname))
            shutil.copy(src_image_dict[color], entry[1])
        return []


class MultiProcessImageDownloader(ImageDownloader):
    """
    Uses multiprocess package to download images in parallel

    """

    def __init__(self, poolsize=1, skip_existing=False,
                 override_dfunc=None):
        """
        Constructor

        .. warning::

            Exceeding **poolsize** of ``4`` causes errors from Human Protein Atlas site

        :param poolsize: Number of concurrent downloaders to use.
        :type poolsize: int
        :param skip_existing: If ``True`` skip download if image file exists and has size
                              greater then ``0``
        :type skip_existing: bool
        :param override_dfunc: Function that takes a tuple `(image URL, download str path)`
                               and downloads the image. If ``None`` :py:func:`download_file`
                               function is used
        :type override_dfunc: :py:class:`function`
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

        .. code-block::

            from cellmaps_imagedownloader.runner import MultiProcessImageDownloader

            dloader = MultiProcessImageDownloader(poolsize=2)

            d_list = [('https://images.proteinatlas.org/992/1_A1_1_red.jpg',
                       '/tmp/1_A1_1_red.jpg')]
            failed = dloader.download_images(download_list=d_list)

        :param download_list: Each tuple of format `(image URL, dest file path)`
        :type download_list: list of tuple
        :return: Failed downloads, format of tuple
                 (`http status code`, `text of error`, (`link`, `destfile`))
        :rtype: list of tuple
        """
        failed_downloads = []
        logger.debug('Poolsize for image downloader set to: ' +
                     str(self._poolsize))
        num_to_download = len(download_list)
        logger.info(str(num_to_download) + ' images to download')
        t = tqdm(total=num_to_download, desc='Download',
                 unit='images')
        if self._poolsize <= 1:
            for entry in download_list:
                t.update()
                res = self._dfunc(entry)
                if res is not None:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug('Failed download: ' + str(res))
                    failed_downloads.append(res)
        else:
            with Pool(processes=self._poolsize) as pool:
                for i in pool.imap_unordered(self._dfunc,
                                             download_list):
                    t.update()
                    if i is not None:
                        if logger.isEnabledFor(logging.DEBUG):
                            logger.debug('Failed download: ' + str(i))
                        failed_downloads.append(i)
        t.close()
        return failed_downloads


class CellmapsImageDownloader(object):
    """
    Downloads Immunofluorescent images from
    `Human Protein Atlas <https://www.proteinatlas.org>`__
    storing them in an output directory that is locally
    registered as an `RO-Crate <https://www.researchobject.org/ro-crate>`__

    """

    SAMPLES_FILEKEY = 'samples'
    UNIQUE_FILEKEY = 'unique'

    def __init__(self, outdir=None,
                 imgsuffix='.jpg',
                 imagedownloader=MultiProcessImageDownloader(),
                 imagegen=None,
                 imageurlgen=None,
                 skip_logging=True,
                 provenance=None,
                 input_data_dict=None,
                 provenance_utils=ProvenanceUtil(),
                 skip_failed=False):
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
        :param image_url: Base URL for image download from Human Protein Atlas
        :type image_url: str
        :param skip_logging: If ``True`` skip logging, if ``None`` or ``False`` do NOT skip logging
        :type skip_logging: bool
        :param provenance:
        :type provenance: dict
        :param input_data_dict:
        :type input_data_dict: dict
        :param provenance_utils: Wrapper for `fairscape-cli <https://pypi.org/project/fairscape-cli>`__
                                 which is used for
                                 `RO-Crate <https://www.researchobject.org/ro-crate>`__ creation and population
        :type provenance_utils: :py:class:`~cellmaps_utils.provenance.ProvenanceUtil`
        """
        if outdir is None:
            raise CellMapsImageDownloaderError('outdir is None')
        self._outdir = os.path.abspath(outdir)
        self._imagedownloader = imagedownloader
        self._imgsuffix = imgsuffix
        self._start_time = int(time.time())
        self._end_time = -1
        self._imagegen = imagegen
        self._imageurlgen = imageurlgen
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
        self._provenance_utils = provenance_utils
        self._skip_failed = skip_failed
        self._image_dataset_ids = None

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
        base_dict = {'name': 'Name for pipeline run',
                     'organization-name': 'Name of lab or group. Ex: Ideker',
                     'project-name': 'Name of funding source or project',
                     'cell-line': 'Name of cell line. Ex: U2OS',
                     'treatment': 'Name of treatment, Ex: untreated',
                     'release': 'Name of release. Example: 0.1 alpha',
                     'gene-set': 'Name of gene set. Example chromatin'}
        if with_ids is not None and with_ids is True:
            guid_dict = ProvenanceUtil.example_dataset_provenance(with_ids=with_ids)
            base_dict.update({CellmapsImageDownloader.SAMPLES_FILEKEY: guid_dict,
                              CellmapsImageDownloader.UNIQUE_FILEKEY: guid_dict})
            return base_dict

        field_dict = ProvenanceUtil.example_dataset_provenance(requiredonly=requiredonly)

        base_dict.update({CellmapsImageDownloader.SAMPLES_FILEKEY: field_dict,
                          CellmapsImageDownloader.UNIQUE_FILEKEY: field_dict})
        return base_dict

    def _update_provenance_with_keywords(self):
        """
        Generates appropriate keywords from provenance data set in constructor

        :return: keywords as str values
        :rtype: list
        """
        if self._provenance is None:
            logger.warning('Provenance is None')
            return
        keywords = []
        for key in ['organization-name', 'project-name', 'release',
                    'cell-line', 'treatment', 'gene-set', 'name']:
            if key in self._provenance:
                keywords.append(self._provenance[key])
        keywords.extend(['IF microscopy', 'images'])
        self._provenance['keywords'] = keywords

    def _update_provenance_with_description(self):
        """
        Gets description from provenance
        :return:
        """
        if self._provenance is None:
            logger.warning('Provenance is None')
            return
        desc = ''
        for key in ['organization-name', 'project-name', 'release',
                    'cell-line', 'treatment', 'gene-set', 'name']:
            if key in self._provenance:
                if desc != '':
                    desc += ' '
                desc += self._provenance[key]
        self._provenance['description'] = desc + ' IF microscopy images'

    def _create_output_directory(self):
        """
        Creates output directory if it does not already exist

        :raises CellmapsDownloaderError: If output directory is None or if directory already exists
        """
        if os.path.isdir(self._outdir):
            raise CellMapsImageDownloaderError(self._outdir + ' already exists')

        os.makedirs(self._outdir, mode=0o755)
        for cur_color in constants.COLORS:
            cdir = os.path.join(self._outdir, cur_color)
            if not os.path.isdir(cdir):
                logger.debug('Creating directory: ' + cdir)
                os.makedirs(cdir,
                            mode=0o755)

    def _register_software(self):
        """
        Registers this tool, self._provenance must exist and have 'keywords' and
        'description' set

        :raises CellMapsProvenanceError: If fairscape call fails
        """
        software_keywords = self._provenance['keywords']
        software_keywords.extend(['tools', cellmaps_imagedownloader.__name__])
        software_description = self._provenance['description'] + \
                               ' ' + \
                               cellmaps_imagedownloader.__description__
        self._softwareid = self._provenance_utils.register_software(self._outdir,
                                                                    name=cellmaps_imagedownloader.__name__,
                                                                    description=software_description,
                                                                    author=cellmaps_imagedownloader.__author__,
                                                                    version=cellmaps_imagedownloader.__version__,
                                                                    file_format='py',
                                                                    keywords=software_keywords,
                                                                    url=cellmaps_imagedownloader.__repo_url__)

    def _register_image_gene_node_attrs(self, fold=1):
        """
        Registers <fold>_image_gene_node_attributes.tsv file with
        ro-crate as a dataset

        :param fold: name of fold
        :type fold: int
        """
        if self._image_gene_attrid is None:
            self._image_gene_attrid = []

        keywords = self._provenance['keywords']
        keywords.extend(['gene', 'attributes', 'file'])
        description = self._provenance['description'] + ' Fold ' + str(fold) + ' Image gene node attributes file'

        data_dict = {'name': cellmaps_imagedownloader.__name__ +
                             ' output file',
                     'description': description,
                     'data-format': 'tsv',
                     'author': cellmaps_imagedownloader.__name__,
                     'version': cellmaps_imagedownloader.__version__,
                     'schema': 'https://raw.githubusercontent.com/fairscape/cm4ai-schemas/main/v0.1.0/cm4ai_schema_imageloader_gene_node_attributes.json',
                     'keywords': keywords,
                     'date-published': date.today().strftime(self._provenance_utils.get_default_date_format_str())}
        src_file = self.get_image_gene_node_attributes_file(fold)
        self._image_gene_attrid.append(self._provenance_utils.register_dataset(self._outdir,
                                                                               source_file=src_file,
                                                                               data_dict=data_dict))

    def _add_dataset_to_crate(self, data_dict=None,
                              source_file=None, skip_copy=True):
        """
        Adds a dataset to rocrate

        :param crate_path: Path to rocrate
        :type crate_path: str
        :param data_dict: information needed to add dataset
        :type data_dict: dict
        :return:
        """
        return self._provenance_utils.register_dataset(self._outdir,
                                                       source_file=source_file,
                                                       data_dict=data_dict,
                                                       skip_copy=skip_copy)

    def _register_computation(self):
        """
        Registers computation with rocrate.
        Current implementation only registers the 1st 2000 images

        """
        if self._image_gene_attrid is not None:
            generated = self._image_gene_attrid
        else:
            generated = []
        if self._image_dataset_ids is not None:
            generated.extend(self._image_dataset_ids)

        keywords = self._provenance['keywords']
        keywords.extend(['computation', 'download'])
        description = self._provenance['description'] + ' run of ' + cellmaps_imagedownloader.__name__

        self._provenance_utils.register_computation(self._outdir,
                                                    name=cellmaps_imagedownloader.__computation_name__,
                                                    run_by=str(self._provenance_utils.get_login()),
                                                    command=str(self._input_data_dict),
                                                    description=description,
                                                    keywords=keywords,
                                                    used_software=[self._softwareid],
                                                    used_dataset=[self._unique_datasetid, self._samples_datasetid],
                                                    generated=generated)

    def _create_run_crate(self):
        """
        Creates rocrate for output directory

        :raises CellMapsProvenanceError: If there is an error
        """
        try:
            self._provenance_utils.register_rocrate(self._outdir,
                                                    name=self._provenance['name'],
                                                    organization_name=self._provenance['organization-name'],
                                                    project_name=self._provenance['project-name'],
                                                    description=self._provenance['description'],
                                                    keywords=self._provenance['keywords'])
        except TypeError as te:
            raise CellMapsImageDownloaderError('Invalid provenance: ' + str(te))
        except KeyError as ke:
            raise CellMapsImageDownloaderError('Key missing in provenance: ' + str(ke))

    def _register_samples_dataset(self):
        """

        :return:
        """
        if 'guid' in self._provenance[CellmapsImageDownloader.SAMPLES_FILEKEY]:
            self._samples_datasetid = self._provenance[CellmapsImageDownloader.SAMPLES_FILEKEY]['guid']
            return

        # if input file for samples was not set then write the samples we
        # have to the output directory and use that path as dataset to register
        if self._input_data_dict is None or \
            CellmapsImageDownloader.SAMPLES_FILEKEY not in self._input_data_dict or \
            self._input_data_dict[CellmapsImageDownloader.SAMPLES_FILEKEY] is None:
            logger.debug('no samples passed in, just write out copy to output directory')
            samples_file = os.path.join(self._outdir, 'samplescopy.csv')
            self._imagegen.write_samples_to_csvfile(csvfile=samples_file)
            skip_samples_copy = True
        else:
            samples_file = self._input_data_dict[CellmapsImageDownloader.SAMPLES_FILEKEY]
            skip_samples_copy = False

        # add samples dataset
        self._samples_datasetid = self._add_dataset_to_crate(
            data_dict=self._provenance[CellmapsImageDownloader.SAMPLES_FILEKEY],
            source_file=samples_file,
            skip_copy=skip_samples_copy)
        logger.debug('Samples dataset id: ' + str(self._samples_datasetid))

    def _register_unique_dataset(self):
        """

        :return:
        """
        if 'guid' in self._provenance[CellmapsImageDownloader.UNIQUE_FILEKEY]:
            self._unique_datasetid = self._provenance[CellmapsImageDownloader.UNIQUE_FILEKEY]['guid']
            return

        # if input file for unique list was not set then write the unique list we
        # have to the output directory and use that path as dataset to register
        if self._input_data_dict is None or \
            CellmapsImageDownloader.UNIQUE_FILEKEY not in self._input_data_dict or \
            self._input_data_dict[CellmapsImageDownloader.UNIQUE_FILEKEY] is None:
            unique_file = os.path.join(self._outdir, 'uniquecopy.csv')
            self._imagegen.write_unique_list_to_csvfile(csvfile=unique_file)
            skip_unique_copy = True
        else:
            unique_file = self._input_data_dict[CellmapsImageDownloader.UNIQUE_FILEKEY]
            skip_unique_copy = False

        # add unique dataset
        self._unique_datasetid = self._add_dataset_to_crate(
            data_dict=self._provenance[CellmapsImageDownloader.UNIQUE_FILEKEY],
            source_file=unique_file,
            skip_copy=skip_unique_copy)
        logger.debug('Unique dataset id: ' + str(self._unique_datasetid))

    def _register_downloaded_images(self):
        """
        Registers all the downloaded images
        :return:
        """
        data_dict = {'name': cellmaps_imagedownloader.__name__ + ' downloaded image file',
                     'description': self._provenance['description'] + ' IF image file',
                     'data-format': self._imgsuffix,
                     'author': '???',
                     'version': '???',
                     'date-published': date.today().strftime(self._provenance_utils.get_default_date_format_str())}

        self._image_dataset_ids = []

        for c in constants.COLORS:
            cntr = 0
            for entry in tqdm(os.listdir(os.path.join(self._outdir, c)),
                              desc='FAIRSCAPE ' + c + ' images registration'):
                if not entry.endswith(self._imgsuffix):
                    continue
                fullpath = os.path.join(self._outdir, c, entry)
                data_dict['name'] = entry + ' ' + c + \
                                    ' channel image'
                if len(data_dict['name']) >= 64:
                    data_dict['name'] = data_dict['name'][:63]

                data_dict['keywords'] = [c, 'IF', 'image']
                self._image_dataset_ids.append(self._add_dataset_to_crate(data_dict=data_dict,
                                                                          source_file=fullpath,
                                                                          skip_copy=True))
                cntr += 1
                if cntr > 25:
                    # Todo: https://github.com/fairscape/fairscape-cli/issues/9
                    logger.error('FAIRSCAPE cannot handle too many images, skipping rest')
                    break

    def _get_color_download_map(self):
        """
        Creates a dict where key is color name and value is directory
        path for files for that color

        ``{'red': '/tmp/foo/red'}``

        :return: map of colors to directory paths
        :rtype: dict
        """
        color_d_map = {}
        for c in constants.COLORS:
            color_d_map[c] = os.path.join(self._outdir, c)
        return color_d_map

    def _get_download_tuples(self):
        """
        Gets download list from **imageurlgen** object set via constructor

        :return: list of (image download URL prefix,
                          file path where image should be written)
        :rtype: list
        """
        dtuples = []
        color_map = self._get_color_download_map()
        for image_url, image_dest in self._imageurlgen.get_next_image_url(color_map):
            dtuples.append((image_url, image_dest))
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
        :return: (int with value of 0 upon success otherwise failure, list of failed downloads)
        :rtype: tuple
        """
        if self._imagedownloader is None:
            raise CellMapsImageDownloaderError('Image downloader is None')

        downloadtuples = self._get_download_tuples()

        failed_downloads = self._imagedownloader.download_images(downloadtuples)
        retry_count = 0
        while len(failed_downloads) > 0 and retry_count < max_retry:
            retry_count += 1
            logger.error(str(len(failed_downloads)) +
                         ' images failed to download. Retrying #' + str(retry_count))

            # try one more time with files that failed
            failed_downloads = self._retry_failed_images(failed_downloads=failed_downloads)

        if len(failed_downloads) > 0 and (self._skip_failed is None or self._skip_failed is False):
            raise CellMapsImageDownloaderError('Failed to download: ' +
                                               str(len(failed_downloads)) + ' images')
        return 0, failed_downloads

    def get_image_gene_node_attributes_file(self, fold):
        """
        Gets full path to image gene node attribute file under output directory
        created when invoking :py:meth:`~cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`

        :return: Path to file
        :rtype: str
        """
        return os.path.join(self._outdir, str(fold) + '_' +
                            constants.IMAGE_GENE_NODE_ATTR_FILE)

    def get_image_gene_node_errors_file(self):
        """
        Gets full path to image gene node attribute errors file under output directory
        created when invoking :py:meth:`~cellmaps_imagedownloader.runner.CellmapsImageDownloader.run`

        :return: Path to file
        :rtype: str
        """
        return os.path.join(self._outdir,
                            constants.IMAGE_GENE_NODE_ERRORS_FILE)

    def _write_image_gene_node_attrs(self, gene_node_attrs=None, fold=1,
                                     errors=None):
        """

        :param gene_node_attrs:
        :param errors:
        :return:
        """
        with open(self.get_image_gene_node_attributes_file(fold), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=constants.IMAGE_GENE_NODE_COLS, delimiter='\t')
            writer.writeheader()
            for key in gene_node_attrs:
                writer.writerow(gene_node_attrs[key])
        if errors is not None:
            with open(self.get_image_gene_node_errors_file(), 'w') as f:
                for e in errors:
                    f.write(str(e) + '\n')

    def _add_imageurl_to_gene_node_attrs(self, gene_node_attrs=None):
        """
        Adds imageurl to **gene_node_attrs** passed in

        :param gene_node_attrs:
        :type gene_node_attrs: dict
        """
        sample_urlmap = self._imageurlgen.get_sample_urlmap()
        for key in gene_node_attrs:
            sample = gene_node_attrs[key]
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + \
                       sample['filename']
            if image_id in sample_urlmap:
                if sample_urlmap[image_id].startswith('http'):
                    gene_node_attrs[key][constants.IMAGE_GENE_NODE_IMAGEURL_COL] = sample_urlmap[image_id]
                else:
                    gene_node_attrs[key][constants.IMAGE_GENE_NODE_IMAGEURL_COL] = 'no image url found'
            else:
                # this should NOT happen, but just in case
                logger.error(image_id + ' not in sample urlmap. setting to no image url found')
                gene_node_attrs[key][constants.IMAGE_GENE_NODE_IMAGEURL_COL] = 'no image url found'

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

            self._update_provenance_with_description()
            self._update_provenance_with_keywords()
            self._create_run_crate()
            self._register_samples_dataset()
            self._register_unique_dataset()

            self._register_software()

            exitcode, failed_downloads = self._download_images()
            # todo need to validate downloaded image data

            # Remove entries from samples that lack a download URL
            self._imagegen.filter_samples_by_sample_urlmap(self._imageurlgen.get_sample_urlmap())

            # write image attribute data
            for fold in [1, 2]:
                image_gene_node_attrs, errors = self._imagegen.get_gene_node_attributes(fold)
                self._add_imageurl_to_gene_node_attrs(gene_node_attrs=image_gene_node_attrs)
                # write image attribute data
                self._write_image_gene_node_attrs(image_gene_node_attrs, fold, errors)

                self._register_image_gene_node_attrs(fold)

            self._register_downloaded_images()

            self._register_computation()

            return exitcode
        finally:
            self._end_time = int(time.time())
            # write a task finish file
            logutils.write_task_finish_json(outdir=self._outdir,
                                            start_time=self._start_time,
                                            end_time=self._end_time,
                                            status=exitcode)
