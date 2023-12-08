
import os
import re
import gzip
import logging
import requests

from tqdm import tqdm
from cellmaps_utils import constants
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError


logger = logging.getLogger(__name__)


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

    def _readline(self, proteinatlas):
        """
        Generator that returns next line of **proteinatlas**

        :param proteinatlas: Path to xml or xml.gz file or URL to download
                             xml or xml.gz file
        :type proteinatlas: str
        :return: next line of file
        :rtype: str
        """
        if os.path.isfile(proteinatlas):
            if proteinatlas.endswith('.gz'):
                with gzip.open(proteinatlas, mode='rt') as f:
                    for line in f:
                        yield line
                return
            with open(proteinatlas, 'r') as f:
                for line in f:
                    yield line
            return
        # use python requests to download the file and then get its results
        local_file = os.path.join(self._outdir,
                                  proteinatlas.split('/')[-1])

        with requests.get(proteinatlas,
                          stream=True) as r:
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
            finally:
                tqdm_bar.close()

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
        Gets `image id` from **image_url** where `image id` is basically
        ``<ANTIBODY_ID>/<IMAGE_ID>_`` and the **image_url** is expected to
        be a URL in this format: ``http://xxxx/<ANTIBODY_ID>/<IMAGE_ID>_blue_red_green.jpg``
        or ``http://xxxx/<ANTIBODY_ID>/<IMAGE_ID>_YYYY_blue_red_green.jpg`` where ``xxxx``
        is url of site usually ``images.proteinatlas.org`` and if set ``YYYY`` is
        something like ``cr5asdjd3983``

        **Example 1:

        ``http://images.proteinatlas.org/4109/1832_C1_2_blue_red_green.jpg`` will
        return as ``4109/1832_C2_``

        **Example 2:

        ``http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue_red_green.jpg``
        will return as ``4109/1843_B2_17_``

        :param image_url:
        :type image_url: str
        :return: id of image
        :rtype: str
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


class ProteinAtlasImageIdToURLMapper(object):
    """
    Given a :py:class:`ProteinAtlasImageUrlReader` instances
    of this class create a dict of image id to image URL
    keeping only entries with image ids set in constructor
    """
    def __init__(self, reader=None,
                 valid_image_ids=None):
        """
        Constructor

        :param reader:
        :param valid_image_ids:
        """
        self._reader = reader
        self._valid_image_ids = valid_image_ids

    def get_image_id_to_url_map(self):
        """

        :return:
        """
        id_url_map = {}
        for image_id, image_url in self._reader.get_next_image_id_and_url():
            if self._valid_image_ids is not None and image_id not in self._valid_image_ids:
                continue
            id_url_map[image_id] = image_url
        return id_url_map


class ImageDownloadTupleGenerator(object):
    """
    Gets URL to download images for given samples
    """
    def __init__(self, samples_list=None,
                 mapper=None):
        """

        :param samples_list:
        """
        self._samples_list = samples_list
        self._mapper = mapper
        self._sample_urlmap = None

    def _populate_sample_urlmap(self):
        """
        Iterates over reader and builds a
        map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg

        :return:
        """
        self._sample_urlmap = self._mapper.get_image_id_to_url_map()

    def get_sample_urlmap(self):
        """
        Gets map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg


        .. note::

            This only returns a map if get_next_image_url() has already been called

        :return: map or ``None``
        :rtype: dict
        """
        return self._sample_urlmap

    def _get_image_prefix_suffix(self, image_url):
        """
        Extracts URL prefix and filename suffix from **image_url**

        :param image_url:
        :type image_url: str
        :return: (image url prefix, suffix ie .jpg)
        :rtype: tuple
        """
        prefix = image_url[:image_url.index('_blue')+1]
        suffix = image_url[image_url.rindex('.'):]
        return prefix, suffix

    def get_next_image_url(self, color_download_map=None):
        """

        :param color_download_map: dict of colors to location on filesystem
                                   ``{'red': '/tmp/foo/red'}``
        :return: list of tuples (image download URL, destination file path)
        :rtype: list
        """
        self._populate_sample_urlmap()

        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' +\
                       sample['if_plate_id'] +\
                       '_' + sample['position'] +\
                       '_' + sample['sample'] + '_'
            if image_id not in self._sample_urlmap:
                logger.error(image_id + ' not in sample map which means '
                                        'no URL was found to acquire '
                                        'said image')
                continue

            image_filename = sample['if_plate_id'] +\
                             '_' + sample['position'] +\
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
        Gets map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg


        .. note::

            This only returns a map if get_next_image_url() has already been called

        :return: map or ``None``
        :rtype: dict
        """
        return self._sample_urlmap

    def _get_image_prefix_suffix(self, image_url):
        """
        Extracts URL prefix and filename suffix from **image_url**
        :param image_url:
        :type image_url: str
        :return: (image url prefix, suffix ie .jpg)
        :rtype: tuple
        """
        prefix = image_url[:image_url.index('_blue')+1]
        suffix = image_url[image_url.rindex('.'):]
        return prefix, suffix

    def get_next_image_url(self, color_download_map=None):
        """

        :param color_download_map: dict of colors to location on filesystem
                                   ``{'red': '/tmp/foo/red'}``
        :return: list of tuples (image download URL, destination file path)
        :rtype: list
        """
        self._populate_sample_urlmap()

        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' +\
                       sample['if_plate_id'] +\
                       '_' + sample['position'] +\
                       '_' + sample['sample'] + '_'
            if image_id not in self._sample_urlmap:
                logger.error(image_id + ' not in sample map')
                continue

            image_filename = sample['if_plate_id'] +\
                             '_' + sample['position'] +\
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
        Gets map of ANTIBODY/PLATE_ID_POSITION_SAMPLE_ => download url of _blue_red_green.jpg


        .. note::

            This only returns a map if get_next_image_url() has already been called

        :return: map or ``None``
        :rtype: dict
        """
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
        self._populate_sample_urlmap()

        for sample in self._samples_list:
            image_id = re.sub('^HPA0*|^CAB0*', '', sample['antibody']) + '/' +\
                       sample['if_plate_id'] +\
                       '_' + sample['position'] +\
                       '_' + sample['sample'] + '_'
            if image_id not in self._sample_urlmap:
                logger.error(image_id + ' not in sample map')
                continue

            image_filename = sample['if_plate_id'] +\
                             '_' + sample['position'] +\
                             '_' + sample['sample'] + '_'
            for c in constants.COLORS:
                sample_url = self._sample_urlmap[image_id]
                (image_url_prefix,
                 image_suffix) = self._get_image_prefix_suffix(sample_url)
                basedir = os.path.dirname(image_url_prefix)

                yield (os.path.join(basedir, c, image_filename + 'z01_' + c + image_suffix),
                       os.path.join(color_download_map[c],
                                    image_filename + c + image_suffix))
