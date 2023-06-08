
import os
import gzip
import requests
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError

"""
Todo: 1) Write proteinatlas streaming class that takes URL, proteinatlas.xml or
         proteinatlas.xml.gz file and returns data in file, one line at a time
      2) Write parser that gets url from <imageUrl> </imageUrl>lines that contain word 'blue' inside
         Take that URL
         http://images.proteinatlas.org/15021/1736_F10_19_cr5805e721b05c3_blue_red_green.jpg
         and build a dict of
         15021/1736_F10_19_ => http://images.proteinatlas.org/15021/1736_F10_19_cr5805e721b05c3_blue_red_green.jpg
      3) For any failed images in the downloader, code should do a similar URL parse and use that image
         in dict for download. Need to denote this alternate download
         (probably add URL for image in image_gene_node_attributes.tsv)


Todo: 1) Downloader should just grab single image ###/####_X##_##_blue_red_green_yellow.jpg and split into
         separate files or store as a single file in images/ directory (requires changing, image_embedding)


Some code frags:

# built image url list from grep "<imageUrl>" of proteinatlas.xml
with open(sys.argv[1], 'r') as f:
    for line in f:

        if 'blue' not in line:
            continue
        line = line.rstrip()
        m = re.search('^.*>(.*)<.*$', line)
        full_url = m.group(1)
        print(full_url)


finding failed image in image url list pulled out of proteinatlas.xml

image_url_dict = {}
# list of URLs
with open('image_url_list.list', 'r') as f:
    for line in f:
        line = line.rstrip()
        end_of_url = re.sub('^.*org\/','', line)
        image_url_dict['_'.join(end_of_url.split('_')[:3]) + '_'] = line

#print(image_url_dict)

image_id_set = set()
# list of failed URLs
with open('failed.list', 'r') as f:
    for line in f:
        line = line.rstrip()
        end_of_url = re.sub('^.*org\/','', line)
        image_id = '_'.join(end_of_url.split('_')[:3]) + '_'
        image_id_set.add(image_id)

for image_id in image_id_set:
    print(image_url_dict[image_id])



"""


class ProteinAtlasReader(object):
    """
    Returns contents of proteinatlas.xml file one
    line at a time
    """
    def __init__(self, outdir):
        """
        Constructor
        """
        self._outdir = outdir

    def readline(self, proteinatlas):
        """
        Generator that returns next line of proteinatlas.xml file

        :param proteinatlas: URL or path to proteinatlas.xml| proteinatlas.xml.gz file
        :type proteinatlas: str
        :return: next line of file
        :rtype: str
        """
        if proteinatlas is None:
            raise CellMapsImageDownloaderError('proteinatlas is None')

        if os.path.isfile(proteinatlas):
            if proteinatlas.endswith('.gz'):
                with gzip.open(proteinatlas, mode='rt') as f:
                    for line in f.readline():
                        yield line
                return
            with open(proteinatlas, 'r') as f:
                for line in f.readline():
                    yield line
            return
        # use python requests to download the file and then get its results
        local_file = os.path.join(self._outdir, proteinatlas.split('/')[-1])
        with requests.get(proteinatlas, stream=True) as r:
            r.raise_for_status()
            with open(local_file, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        self.readline(local_file)
