import os
import re
import gzip
import logging
import shutil
import time
import xml.etree.ElementTree as ET
import pandas as pd

import requests
from requests import RequestException

from tqdm import tqdm
from cellmaps_utils import constants
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError

logger = logging.getLogger(__name__)


def download_proteinalas_file(outdir, proteinatlas, max_retries=3, retry_wait=10):
    # use python requests to download the file and then get its results
    local_file = os.path.join(outdir,
                              proteinatlas.split('/')[-1])

    retry_num = 0
    download_successful = False
    while retry_num < max_retries and not download_successful:
        try:
            with requests.get(proteinatlas, stream=True) as r:
                content_size = int(r.headers.get('content-length', 0))
                tqdm_bar = tqdm(desc='Downloading ' + os.path.basename(local_file),
                                total=content_size,
                                unit='B', unit_scale=True,
                                unit_divisor=1024)
                logger.debug('Downloading ' + str(proteinatlas) +
                             ' of size ' + str(content_size) +
                             'b to ' + local_file)
                try:
                    r.raise_for_status()
                    with open(local_file, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
                            tqdm_bar.update(len(chunk))
                    download_successful = True
                finally:
                    tqdm_bar.close()

            return local_file

        except RequestException as he:
            logger.debug(str(he.response.text))
            retry_num += 1
            time.sleep(retry_wait)

    if not download_successful:
        raise CellMapsImageDownloaderError(f'{max_retries} attempts to download proteinatlas file failed')


class ProteinAtlasProcessor(object):

    def __init__(self, outdir=None,
                 proteinatlas=None,
                 proteinlist_file=None):

        self._outdir = outdir
        if proteinatlas is None:
            self._proteinatlas = ProteinAtlasReader.DEFAULT_PROTEINATLAS_URL
        else:
            self._proteinatlas = proteinatlas
        self._proteinlist_file = proteinlist_file

    def _create_output_directory(self):
        if os.path.isdir(self._outdir):
            raise CellMapsImageDownloaderError(self._outdir + ' already exists')
        else:
            os.makedirs(self._outdir, mode=0o755)

    def _read_protein_list(self):
        protein_list = list()
        with open(self._proteinlist_file, 'r') as f:
            for line in f:
                protein_list.append(line.strip())

        self.protein_set = set(protein_list)

    @staticmethod
    def _decompress_gz_file(gz_file, output_file):
        """
        Decompress a .gz file.
        """
        with gzip.open(gz_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return output_file

    def get_sample_list_from_hpa(self):
        self._create_output_directory()
        self._read_protein_list()

        if not os.path.isfile(self._proteinatlas):
            proteinatlas_gz_file = download_proteinalas_file(self._outdir, self._proteinatlas)
            proteinatlas_xml_file = proteinatlas_gz_file.replace('.gz', '')
            self._proteinatlas = self._decompress_gz_file(proteinatlas_gz_file, proteinatlas_xml_file)

        data = []

        context = ET.iterparse(self._proteinatlas, events=('end',))
        for event, elem in context:
            if elem.tag == 'entry':
                gene_name = elem.find('name').text if elem.find('name') is not None else None
                if gene_name not in self.protein_set:
                    elem.clear()
                    continue

                antibody_elem = elem.find('antibody')
                if antibody_elem is None:
                    elem.clear()
                    continue

                antibody = antibody_elem.attrib['id'] if elem.find('antibody') is not None else None

                image_data = antibody_elem.find('cellExpression')
                if image_data is None or image_data.attrib['source'] != "HPA":
                    elem.clear()
                    continue

                data_points = image_data.findall('.//data')
                for dp in data_points:
                    cellline = dp.find('cellLine').text
                    locations = [location.text for location in dp.findall('.//location')]
                    image_urls = [image.find('imageUrl').text for image in dp.findall('.//image') if
                                  'blue' in image.find('imageUrl').text]
                    if len(image_urls) == 0:
                        elem.clear()
                        continue

                    for im_url in image_urls:
                        filename = im_url.split("/")[-1].replace('blue_red_green.jpg', '')
                        values_file = filename.split("_")
                        if_plate, position, sample = values_file[0], values_file[1], values_file[2]
                        data.append({
                            'filename': filename,
                            'if_plate_id': if_plate,
                            'position': position,
                            'sample': sample,
                            'status': 35,
                            'locations': ', '.join(locations),
                            'antibody': antibody,
                            'ensembl_ids': elem.attrib['url'].split('/')[-1],
                            'gene_names': gene_name,
                            'atlas_name': cellline,
                            'image_urls': im_url
                        })

                elem.clear()

        df = pd.DataFrame(data)
        samples_file = os.path.join(self._outdir, 'samples.csv')
        df.to_csv(samples_file)

        return samples_file, self._proteinatlas


class ProteinAtlasReader(object):
    """
    Returns contents of proteinatlas.xml file one
    line at a time
    """

    DEFAULT_PROTEINATLAS_URL = 'https://www.proteinatlas.org/download/proteinatlas.xml.gz'

    def __init__(self, outdir=None,
                 proteinatlas=None):
        """
        Constructor

        :param outdir: Path to directory where results can be written to
        :type outdir: str
        :param proteinatlas: URL or path to proteinatlas.xml| proteinatlas.xml.gz file
        :type proteinatlas: str
        """
        self._outdir = outdir
        if proteinatlas is None:
            self._proteinatlas = ProteinAtlasReader.DEFAULT_PROTEINATLAS_URL
        else:
            self._proteinatlas = proteinatlas

    def readline(self):
        """
        Generator that returns next line of proteinatlas data
        set via constructor

        :return: next line of file
        :rtype: str
        """
        for line in self._readline(self._proteinatlas):
            yield line

    def _readline(self, proteinatlas, max_retries=3, retry_wait=10):
        """
        Generator that returns next line of **proteinatlas**

        :param proteinatlas: Path to xml or xml.gz file or URL to download
                             xml or xml.gz file
        :type proteinatlas: str
        :param max_retries:
        :type max_retries: int
        :param retry_wait:
        :type retry_wait: int
        :return: next line of file
        :rtype: str
        """
        if os.path.isfile(proteinatlas):
            if proteinatlas.endswith('.gz'):
                with gzip.open(proteinatlas, mode='rt') as f:
                    for line in tqdm(f, desc='Processing proteinatlas.xml.gz', unit='bytes',
                                     total=os.path.getsize(proteinatlas)):
                        yield line
                return
            with open(proteinatlas, 'r') as f:
                for line in tqdm(f, desc='Processing proteinatlas.xml', unit='bytes',
                                 total=os.path.getsize(proteinatlas)):
                    yield line
            return

        local_file = download_proteinalas_file(self._outdir, proteinatlas, max_retries, retry_wait)

        for line in self._readline(local_file):
            yield line


class ProteinAtlasImageUrlReader(object):
    """
    Takes a proteinatlas generator to get
    value between <imageUrl>XXX</imageUrl> lines
    with the keyword _blue in them
    """

    def __init__(self, reader=None):
        """
        Constructor

        :param reader:
        :type reader: :py:class:`~cellmaps_imagedownloader.proteinatlas.ProteinAtlasReader`
        """
        self._reader = reader

    def _get_url_from_line(self, line):
        """
        Gets value between rightmost >< in text file

        :param line:
        :type line: str
        :return: content between > and < characters
        :rtype: str
        """
        m = re.search('^.*>(.*)<.*$', line)
        return m.group(1)

    def _get_image_id(self, image_url):
        """

        :param image_url:
        :return:
        """
        antibody_and_id = '/'.join(image_url.split('/')[-2:])
        return '_'.join(antibody_and_id.split('_')[0:3]) + '_'

    def get_next_image_id_and_url(self):
        """

        :return: (image id, image_url)
        :rtype: tuple
        """
        for line in self._reader.readline():
            line = line.rstrip()
            if '<imageUrl>' not in line:
                continue
            if 'blue' not in line:
                continue
            image_url = self._get_url_from_line(line)
            yield self._get_image_id(image_url), image_url


class ImageDownloadTupleGenerator(object):
    """
    Gets URL to download images for given samples
    """

    def __init__(self, samples_list=None,
                 reader=None,
                 valid_image_ids=None):
        """
        Constructor

        :param samples_list:
        :type samples_list: list
        :param reader: Used to get download URLs for images
        :type reader: :py:class:`~cellmaps_imagedownloader.proteinatlas.ProteinAtlasImageUrlReader`
        :param valid_image_ids: Image ids that need a download URL in format of
                                ``<ANTIBODY ID minus HPA or CAB prefix>/<IMAGE ID>``
        :type valid_image_ids: set
        """
        self._samples_list = samples_list
        self._reader = reader
        self._sample_urlmap = None
        self._valid_image_ids = valid_image_ids

    def _populate_sample_urlmap(self):
        """
        Iterates over reader and builds a
        map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg
        for all images in ``self._valid_image_ids`` list or for everything if that is ``None``

        The map is stored in ``self._sample_urlmap``
        """
        self._sample_urlmap = {}
        for image_id, image_url in self._reader.get_next_image_id_and_url():
            if self._valid_image_ids is not None and image_id not in self._valid_image_ids:
                continue
            self._sample_urlmap[image_id] = image_url

    def get_sample_urlmap(self):
        """
        Gets map of ``ANTIBODY/PLATE_ID_POSITION_SAMPLE_`` => download url of _blue_red_green.jpg

        :return: map or ``None``
        :rtype: dict
        """
        if self._sample_urlmap is None:
            self._populate_sample_urlmap()
        return self._sample_urlmap

    def _get_image_prefix_suffix(self, image_url):
        """
        Extracts URL prefix and filename suffix from **image_url**

        :param image_url: URL to download image
        :type image_url: str
        :return: (image url prefix, suffix ie .jpg)
        :rtype: tuple
        """
        prefix = image_url[:image_url.index('_blue') + 1]
        suffix = image_url[image_url.rindex('.'):]
        return prefix, suffix

    def get_next_image_url(self, color_download_map=None):
        """
        `Generator function <https://docs.python.org/3/glossary.html#index-19>`__
        that gets the next image URL to download

        :param color_download_map: dict of colors to location on filesystem
                                   ``{'red': '/tmp/foo/red'}``
        :return: list of tuples (image download URL, destination file path)
        :rtype: list
        """
        if self._sample_urlmap is None:
            self._populate_sample_urlmap()

        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + \
                       sample['if_plate_id'] + \
                       '_' + sample['position'] + \
                       '_' + sample['sample'] + '_'
            if image_id not in self._sample_urlmap:
                logger.error(image_id + ' not in sample map which means '
                                        'no URL was found to acquire '
                                        'said image')
                continue

            image_filename = sample['if_plate_id'] + \
                             '_' + sample['position'] + \
                             '_' + sample['sample'] + '_'
            for c in constants.COLORS:
                sample_url = self._sample_urlmap[image_id]
                (image_url_prefix,
                 image_suffix) = self._get_image_prefix_suffix(sample_url)

                yield (image_url_prefix + c + image_suffix,
                       os.path.join(color_download_map[c],
                                    image_filename + c + image_suffix))


class LinkPrefixImageDownloadTupleGenerator(object):
    """
    Gets URL to download images for given samples
    """

    def __init__(self, samples_list=None):
        """

        :param samples_list:
        """
        self._samples_list = samples_list
        self._sample_urlmap = None

    def _populate_sample_urlmap(self):
        """
        Iterates over reader and builds a
        map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg

        :return:
        """
        self._sample_urlmap = {}
        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + \
                       sample['filename']

            self._sample_urlmap[image_id] = sample['linkprefix'] + 'blue_red_green.jpg'

        logger.debug(self._sample_urlmap)

    def get_sample_urlmap(self):
        """
        Gets map of ``ANTIBODY/PLATE_ID_POSITION_SAMPLE_`` => download url of _blue_red_green.jpg

        :return: map or ``None``
        :rtype: dict
        """
        if self._sample_urlmap is None:
            self._populate_sample_urlmap()
        return self._sample_urlmap

    def _get_image_prefix_suffix(self, image_url):
        """
        Extracts URL prefix and filename suffix from **image_url**
        :param image_url:
        :type image_url: str
        :return: (image url prefix, suffix ie .jpg)
        :rtype: tuple
        """
        prefix = image_url[:image_url.index('_blue') + 1]
        suffix = image_url[image_url.rindex('.'):]
        return prefix, suffix

    def get_next_image_url(self, color_download_map=None):
        """

        :param color_download_map: dict of colors to location on filesystem
                                   ``{'red': '/tmp/foo/red'}``
        :return: list of tuples (image download URL, destination file path)
        :rtype: list
        """
        if self._sample_urlmap is None:
            self._populate_sample_urlmap()

        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + \
                       sample['if_plate_id'] + \
                       '_' + sample['position'] + \
                       '_' + sample['sample'] + '_'
            if image_id not in self._sample_urlmap:
                logger.error(image_id + ' not in sample map')
                continue

            image_filename = sample['if_plate_id'] + \
                             '_' + sample['position'] + \
                             '_' + sample['sample'] + '_'
            for c in constants.COLORS:
                sample_url = self._sample_urlmap[image_id]
                (image_url_prefix,
                 image_suffix) = self._get_image_prefix_suffix(sample_url)

                yield (image_url_prefix + c + image_suffix,
                       os.path.join(color_download_map[c],
                                    image_filename + c + image_suffix))


class CM4AIImageCopyTupleGenerator(object):
    """
    Gets URL to download images for given samples
    """

    def __init__(self, samples_list=None):
        """

        :param samples_list:
        """
        self._samples_list = samples_list
        self._sample_urlmap = None

    def _populate_sample_urlmap(self):
        """
        Iterates over reader and builds a
        map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg

        :return:
        """
        self._sample_urlmap = {}
        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + \
                       sample['filename']

            self._sample_urlmap[image_id] = os.path.join(sample['linkprefix'],
                                                         sample['filename'] + 'z01_blue.jpg')

        logger.debug(self._sample_urlmap)

    def get_sample_urlmap(self):
        """
        Gets map of ``ANTIBODY/PLATE_ID_POSITION_SAMPLE_`` => download url of _blue_red_green.jpg

        :return: map or ``None``
        :rtype: dict
        """
        if self._sample_urlmap is None:
            self._populate_sample_urlmap()
        return self._sample_urlmap

    def _get_image_prefix_suffix(self, image_url):
        """
        Extracts URL prefix and filename suffix from **image_url**
        :param image_url:
        :type image_url: str
        :return: (image url prefix, suffix ie .jpg)
        :rtype: tuple
        """
        prefix = image_url[:image_url.index('_blue') + 1]
        suffix = image_url[image_url.rindex('.'):]
        return prefix, suffix

    def get_next_image_url(self, color_download_map=None):
        """

        :param color_download_map: dict of colors to location on filesystem
                                   ``{'red': '/tmp/foo/red'}``
        :return: list of tuples (image download URL, destination file path)
        :rtype: list
        """
        if self._sample_urlmap is None:
            self._populate_sample_urlmap()

        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' + \
                       sample['if_plate_id'] + \
                       '_' + sample['position'] + \
                       '_' + sample['sample'] + '_'
            if image_id not in self._sample_urlmap:
                logger.error(image_id + ' not in sample map')
                continue

            image_filename = sample['if_plate_id'] + \
                             '_' + sample['position'] + \
                             '_' + sample['sample'] + '_'
            for c in constants.COLORS:
                sample_url = self._sample_urlmap[image_id]
                (image_url_prefix,
                 image_suffix) = self._get_image_prefix_suffix(sample_url)
                basedir = os.path.dirname(image_url_prefix)

                yield (os.path.join(basedir, c, image_filename + 'z01_' + c + image_suffix),
                       os.path.join(color_download_map[c],
                                    image_filename + c + image_suffix))
