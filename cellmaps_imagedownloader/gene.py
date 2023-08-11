
import re
import csv
import mygene
import logging
from tqdm import tqdm

from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError

logger = logging.getLogger(__name__)


class GeneQuery(object):
    """
    Gets information about genes from mygene
    """
    def __init__(self, mygeneinfo=mygene.MyGeneInfo()):
        """
        Constructor
        """
        self._mg = mygeneinfo

    def querymany(self, queries, species=None,
                  scopes=None,
                  fields=None):

        """
        Simple wrapper that calls MyGene querymany
        returning the results

        :param queries: list of gene ids/symbols to query
        :type queries: list
        :param species:
        :type species: str
        :param scopes:
        :type scopes: str
        :param fields:
        :type fields: list
        :return: dict from MyGene usually in format of
        :rtype: list
        """
        mygene_out = self._mg.querymany(queries,
                                        scopes=scopes,
                                        fields=fields,
                                        species=species)
        return mygene_out

    def get_symbols_for_genes(self, genelist=None,
                              scopes='_id'):
        """
        Queries for genes via GeneQuery() object passed in via
        constructor

        :param genelist: genes to query for valid symbols and ensembl ids
        :type genelist: list
        :param scopes: field to query on _id for gene id, ensemble.gene
                       for ENSEMBLE IDs
        :type scopes: str
        :return: result from mygene which is a list of dict objects where
                 each dict is of format:

                 .. code-block::

                     { 'query': 'ID',
                       '_id': 'ID', '_score': #.##,
                       'ensembl': { 'gene': 'ENSEMBLEID' },
                       'symbol': 'GENESYMBOL' }
        :rtype: list
        """
        res = self.querymany(genelist,
                             species='human',
                             scopes=scopes,
                             fields=['ensembl.gene', 'symbol'])
        return res


class GeneNodeAttributeGenerator(object):
    """
    Base class for GeneNodeAttribute Generator
    """
    def __init__(self):
        """
        Constructor
        """
        pass

    def get_gene_node_attributes(self):
        """
        Should be implemented by subclasses

        :raises NotImplementedError: Always
        """
        raise NotImplementedError('Subclasses should implement')


class ImageGeneNodeAttributeGenerator(GeneNodeAttributeGenerator):
    """
    Creates Image Gene Node Attributes table
    """

    SAMPLES_HEADER_COLS = ['filename', 'if_plate_id',
                           'position', 'sample', #'status',
                           'locations', 'antibody',
                           'ensembl_ids', 'gene_names']
    LINKPREFIX_HEADER = 'linkprefix'
    """
    Column labels for samples file
    """

    UNIQUE_HEADER_COLS = ['antibody', 'ensembl_ids',
                          'gene_names', 'atlas_name',
                          'locations',
                          'n_location']
    """
    Column labels for unique file
    """

    def __init__(self, samples_list=None,
                 unique_list=None,
                 genequery=GeneQuery()):
        """
        Constructor

        **samples_list** is expected to be a list of :py:class:`dict`
        objects with this format:

        # TODO: Move this to a separate data document

        .. code-block::

            {
             'filename': HPA FILENAME,
             'if_plate_id': HPA PLATE ID,
             'position': POSITION,
             'sample': SAMPLE,
             'status': STATUS,
             'locations': COMMA DELIMITED LOCATIONS,
             'antibody': ANTIBODY_ID,
             'ensembl_ids': COMMA DELIMITED ENSEMBLID IDS,
             'gene_names': COMMA DELIMITED GENE SYMBOLS
            }

        **Example:**

        .. code-block::

            {
             'filename': '/archive/1/1_A1_1_',
             'if_plate_id': '1',
             'position': 'A1',
             'sample': '1',
             'status': '35',
             'locations': 'Golgi apparatus',
             'antibody': 'HPA000992',
             'ensembl_ids': 'ENSG00000066455',
             'gene_names': 'GOLGA5'
            }

        **unique_list** is expected to be a list of :py:class:`dict`
        objects with this format:

        .. code-block::

            {
             'antibody': ANTIBODY,
             'ensembl_ids': COMMA DELIMITED ENSEMBL IDS,
             'gene_names': COMMA DELIMITED GENE SYMBOLS,
             'atlas_name': ATLAS NAME?,
             'locations': COMMA DELIMITED LOCATIONS IN CELL,
             'n_location': NUMBER OF LOCATIONS IN CELL,
             }

        **Example:**

        .. code-block::

            {
             'antibody': 'HPA040086',
             'ensembl_ids': 'ENSG00000094914',
             'gene_names': 'AAAS',
             'atlas_name': 'U-2',
             'locations': 'OS,Nuclear membrane',
             'n_location': '2',
             }


        :param samples_list: List of samples
        :type samples_list: list
        :param unique_list: List of unique samples
        :type unique_list: list
        :param genequery: Object to query for updated gene symbols
        :type genequery: :py:class:`~cellmaps_imagedownloader.gene.GeneQuery`
        """
        super().__init__()
        self._samples_list = samples_list
        self._unique_list = unique_list
        self._genequery = genequery

    def get_samples_list(self):
        """
        Gets **samples_list** passed in via the constructor


        :return: list of samples set via constructor
        :rtype: list
        """
        return self._samples_list

    @staticmethod
    def get_samples_from_csvfile(csvfile=None):
        """

        :param tsvfile:
        :return:
        """
        if csvfile is None:
            raise CellMapsImageDownloaderError('csvfile is None')

        samples = []
        with open(csvfile, 'r') as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                sample_entry = {}
                for key in ImageGeneNodeAttributeGenerator.SAMPLES_HEADER_COLS:
                    sample_entry[key] = row[key]
                if ImageGeneNodeAttributeGenerator.LINKPREFIX_HEADER in row:
                    sample_entry[ImageGeneNodeAttributeGenerator.LINKPREFIX_HEADER] = row[ImageGeneNodeAttributeGenerator.LINKPREFIX_HEADER]
                samples.append(sample_entry)
        return samples

    def get_unique_list(self):
        """
        Gets antibodies_list passed in via the constructor

        :return:
        """
        return self._unique_list

    @staticmethod
    def get_unique_list_from_csvfile(csvfile=None):
        """

        :param csvfile:
        :return:
        """
        if csvfile is None:
            raise CellMapsImageDownloaderError('csvfile is None')

        u_list = []
        with open(csvfile, 'r') as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                unique_entry = {}
                for key in ImageGeneNodeAttributeGenerator.UNIQUE_HEADER_COLS:
                    unique_entry[key] = row[key]
                u_list.append(unique_entry)
        return u_list

    def write_unique_list_to_csvfile(self, csvfile=None):
        """
        Writes unique list to file

        :param csvfile: path to file to write
        :type csvfile: str
        """
        if csvfile is None:
            raise CellMapsImageDownloaderError('csvfile is None')
        with open(csvfile, 'w', newline='') as f:
            writer = csv.DictWriter(f,
                                    fieldnames=ImageGeneNodeAttributeGenerator.UNIQUE_HEADER_COLS,
                                    delimiter=',')
            writer.writeheader()
            for unique_entry in self._unique_list:
                writer.writerow(unique_entry)

    def write_samples_to_csvfile(self, csvfile=None):
        """
        Writes samples to file

        :param csvfile: path to file to write
        :type csvfile: str
        """
        if csvfile is None:
            raise CellMapsImageDownloaderError('csvfile is None')
        with open(csvfile, 'w', newline='') as f:
            writer = csv.DictWriter(f,
                                    fieldnames=ImageGeneNodeAttributeGenerator.SAMPLES_HEADER_COLS,
                                    delimiter=',')
            writer.writeheader()
            for sample in self._samples_list:
                writer.writerow(sample)

    def _get_unique_ids_from_samplelist(self, column='ensembl_ids'):
        """
        Gets a unique list of ids split by comma from the samples
        under **column**.

        For example for a sample with these values and column set to ``ensembl_ids``:

        .. code-block:: python

            {'ensembl_ids': 'ENSG00000240682,ENSG00000261796'}

        The resulting tuple would be:

        .. code-block:: python

            ['ENSG00000240682', 'ENSG00000261796']

        :return: split values from corresponding **column** in samples list
        :rtype: list
        """
        id_set = set()

        for row in self._samples_list:
            geneid=row[column]
            split_str = re.split('\W*,\W*', geneid)
            id_set.update(split_str)

        return list(id_set)

    def _get_set_of_antibodies_from_unique_list(self):
        """
        Extract a unique set of antibodies from antibodies list
        passed in via constructor

        :return: unique antibodies
        :rtype: set
        """
        if self._unique_list is None:
            raise CellMapsImageDownloaderError('unique list is None')

        antibody_set = set()
        for a in self._unique_list:
            if 'antibody' not in a:
                logger.warning('Skipping because antibody not found '
                               'in unique entry: ' + str(a))
                continue
            antibody_set.add(a['antibody'])
        return antibody_set

    def get_dicts_of_gene_to_antibody_filename(self, allowed_antibodies=None):
        """
        Gets a tuple of dictionaries from the sample list passed in via
        the constructor.


        :param allowed_antibodies: Skip samples whose antibody is NOT in this list.
                                   If ``None`` then all samples are included
        :type allowed_antibodies: list or set
        :return: (:py:class:`dict` of ensembl_id => antibody,
                  :py:class:`dict` of antibody => filename,
                  :py:class:`dict` of antibody => comma delimited ambiguous ensembl_ids)

        :rtype: tuple
        """
        if self._samples_list is None:
            raise CellMapsImageDownloaderError('samples list is None')

        g_antibody_dict = {}
        antibody_filename_dict = {}
        ambiguous_antibody_dict = {}

        for sample in self._samples_list:
            antibody = sample['antibody']
            if allowed_antibodies is not None and antibody not in allowed_antibodies:
                # skipping cause antibody is not in allowed set
                continue

            ensembl_ids = sample['ensembl_ids'].split(',')
            if len(ensembl_ids) > 1:
                ambiguous_antibody_dict[antibody] = sample['ensembl_ids']

            if antibody not in antibody_filename_dict:
                antibody_filename_dict[antibody] = set()
            antibody_filename_dict[antibody].add(sample['if_plate_id'] + '_' +
                                   sample['position'] + '_' +
                                   sample['sample'] + '_')
            for g in ensembl_ids:
                #if gene already has nonambgiuous antibody, use that one
                if g in g_antibody_dict:
                    if g in ambiguous_antibody_dict:
                        continue
                g_antibody_dict[g] = sample['antibody']

        return g_antibody_dict, antibody_filename_dict, ambiguous_antibody_dict

    def get_gene_node_attributes(self, fold=1):
        """
        Using **samples_list** and **unique_list**, builds
        a list of :py:class:`dict` objects with updated Gene Symbols.

        Format of each resulting :py:class:`dict`:

        .. code-block::

            {'name': GENE_SYMBOL,
             'represents': ENSEMBL_ID,
             'ambiguous': AMBIGUOUS_GENES,
             'antibody': ANTIBODY,
             'filename': FILENAME}

        **Example**

        .. code-block::

            {'ENSG00000066455': {'name': 'GOLGA5',
                                 'represents': 'ensembl:ENSG00000066455',
                                 'ambiguous': '',
                                 'antibody': 'HPA000992',
                                 'filename': '1_A1_2_,1_A1_1_'}}

        :return: (list of dict, list of errors)
        :rtype: tuple
        """
        # Todo: Refactor because this method is doing waaay too much

        t = tqdm(total=5, desc='Get updated gene symbols',
                 unit='steps')

        t.update()
        # get the unique set of ensembl_ids for mygene query
        ensembl_id_list = self._get_unique_ids_from_samplelist()

        t.update()

        # queries mygene and gets a list of dicts that look like this:
        # {'query': 'ENSG00000066455',
        #  '_id': '9950',
        #  '_score': 25.046944,
        #  'ensembl': {'gene':'ENSG00000066455'},
        #  'symbol': 'GOLGA5'
        # }
        query_res = self._genequery.get_symbols_for_genes(genelist=ensembl_id_list,
                                                          scopes='ensembl.gene')

        t.update()

        # get the unique or best antibodies to use
        unique_antibodies = self._get_set_of_antibodies_from_unique_list()

        t.update()
        # create a mapping of ensembl id to antibody and ensembl_id to filenames
        # where entries NOT in unique_antibodies are filtered out
        # get mapping of ambiguous genes
        g_antibody_dict, antibody_filename_dict, ambiguous_antibody_dict = self.get_dicts_of_gene_to_antibody_filename(allowed_antibodies=unique_antibodies)

        errors = []
        gene_node_attrs = {}
        for x in query_res:

            # skips item that lacks a symbol like this one:
            # {'query': 'ENSG00000282988',
            #  '_id': 'ENSG00000282988',
            #  '_score': 25.04868,
            #  'ensembl': {'gene': 'ENSG00000282988'}}
            if 'symbol' not in x:
                errors.append('Skipping ' + str(x) +
                              ' no symbol in query result: ' + str(x))
                logger.error(errors[-1])
                continue

            ensemblstr = 'ensembl:'

            # skips item that lacks anything like this one:
            # {'query': 'ENSG000001', 'notfound': True}
            if 'ensembl' not in x:
                errors.append('Skipping ' + str(x) +
                              ' no ensembl in query result: ' + str(x))
                logger.error(errors[-1])
                continue

            ensembl_id = None

            # check if item 'ensembl' has more then 1 element
            #
            # {'query': 'ENSG00000273706',
            #  '_id': '3975', '_score': 24.515644,
            #  'ensembl': [{'gene': 'ENSG00000273706'},
            #              {'gene': 'ENSG00000274577'}],
            #  'symbol': 'LHX1'}
            if len(x['ensembl']) > 1:
                for g in x['ensembl']:
                    if g['gene'] in g_antibody_dict:
                        # we need an ensembl id for filename and antibody
                        # lookup so just grab the 1st one
                        ensembl_id = g['gene']
                        break
                # concatenate ensembl ids and delimit with ;
                ensemblstr += ';'.join([g['gene'] for g in x['ensembl']])
            else:
                ensemblstr += x['ensembl']['gene']
                ensembl_id = x['ensembl']['gene']

            if ensembl_id not in g_antibody_dict:
                continue
            antibody_str = g_antibody_dict[ensembl_id]

            filenames = list(antibody_filename_dict[antibody_str])
            if len(filenames) < fold:
                filename_str = filenames[0]
            else:
                filename_str = filenames[fold - 1]

            ambiguous_str = ''
            if antibody_str in ambiguous_antibody_dict:
                ambiguous_str = ambiguous_antibody_dict[antibody_str]

            gene_node_attrs[x['query']] = {'name': x['symbol'],
                                           'represents': ensemblstr,
                                           'ambiguous': ambiguous_str,
                                           'antibody': antibody_str,
                                           'filename': filename_str}

        return gene_node_attrs, errors
