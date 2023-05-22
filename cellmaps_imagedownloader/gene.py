
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

    @staticmethod
    def add_geneids_to_set(gene_set=None,
                           ambiguous_gene_dict=None,
                           geneid=None):
        """
        Examines **geneid** passed in and if a comma exists
        in value split by comma and assume multiple genes.
        Adds those genes into **gene_set** and add entry
        to **ambiguous_gene_dict** with key set to each gene
        name and value set to original **geneid** value

        :param gene_set: unique set of genes
        :type gene_set: set
        :param geneid: name of gene or comma delimited string of genes
        :type geneid: str
        :return: genes found in **geneid** or None if **gene_set**
                 or **geneid** is ``None``
        :rtype: list
        """
        if gene_set is None:
            return None
        if geneid is None:
            return None

        split_str = re.split('\W*,\W*', geneid)
        gene_set.update(split_str)
        if ambiguous_gene_dict is not None:
            if len(split_str) > 1:
                for entry in split_str:
                    ambiguous_gene_dict[entry] = geneid
        return split_str

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
                           'position', 'sample', 'status',
                           'locations', 'antibody',
                           'ensembl_ids', 'gene_names']
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

    def _get_set_of_antibodies_from_unique_list(self):
        """
        Extract a unique set of antibodies from antibodies list
        passed in via constructor

        :return:
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
        Gets a dictionary where key is ensembl id and value is
        the file_name value

        :return:
        :rtype: dict
        """
        if self._samples_list is None:
            raise CellMapsImageDownloaderError('samples list is None')

        g_antibody_dict = {}
        g_filename_dict = {}

        for sample in self._samples_list:
            if allowed_antibodies is not None and sample['antibody'] not in allowed_antibodies:
                # skipping cause antibody is not in allowed set
                continue

            ensembl_ids = sample['ensembl_ids'].split(',')
            for g in ensembl_ids:
                if g not in g_antibody_dict:
                    g_antibody_dict[g] = set()
                if g not in g_filename_dict:
                    g_filename_dict[g] = set()
                g_antibody_dict[g].add(sample['antibody'])
                g_filename_dict[g].add(sample['if_plate_id'] + '_' +
                                       sample['position'] + '_' +
                                       sample['sample'] + '_')

        return g_antibody_dict, g_filename_dict

    def _get_unique_ids_from_samplelist(self, column='ensembl_ids'):
        """
        Gets a unique list of ids split by comma from the samples
        under **column**. In addition a dict is also created where
        key is split id and value is original unsplit values

        For example for a sample with these values and column set to ``ensembl_ids``:

        .. code-block:: python

            {'ensembl_ids': 'ENSG00000240682,ENSG00000261796'}

        The resulting tuple would be:

        .. code-block:: python

            (['ENSG00000240682', 'ENSG00000261796'],
             {'ENSG00000240682': 'ENSG00000240682,ENSG00000261796',
              'ENSG00000261796': 'ENSG00000240682,ENSG00000261796'})

        :return: (list of ids, dict where key is id and value is original unsplit value)
        :rtype: tuple
        """
        id_set = set()
        ambiguous_id_dict = {}

        for row in self._samples_list:
            GeneNodeAttributeGenerator.add_geneids_to_set(gene_set=id_set,
                                                          ambiguous_gene_dict=ambiguous_id_dict,
                                                          geneid=row[column])
        return list(id_set), ambiguous_id_dict

    def get_gene_node_attributes(self):
        """
        TODO: need to implement this

        :return:
        """
        t = tqdm(total=5, desc='Get updated gene symbols',
                 unit='steps')

        t.update()
        # get the unique set of ensembl_ids for mygene query
        ensembl_id_list, _ = self._get_unique_ids_from_samplelist()

        t.update()
        query_res = self._genequery.get_symbols_for_genes(genelist=ensembl_id_list,
                                                          scopes='ensembl.gene')
        t.update()
        # get mapping of ambiguous genes
        _, ambiguous_gene_dict = self._get_unique_ids_from_samplelist(column='gene_names')

        t.update()
        # get the unique or best antibodies to use
        unique_antibodies = self._get_set_of_antibodies_from_unique_list()

        t.update()
        # create a mapping of ensembl id to antibody and ensembl_id to filenames
        # where entries NOT in unique_antibodies are filtered out
        g_antibody_dict, g_filename_dict = self.get_dicts_of_gene_to_antibody_filename(allowed_antibodies=unique_antibodies)

        errors = []
        gene_node_attrs = {}
        for x in query_res:
            if 'symbol' not in x:
                errors.append('Skipping ' + str(x) +
                              ' no symbol in query result: ' + str(x))
                logger.error(errors[-1])
                continue
            ensemblstr = 'ensembl:'
            if 'ensembl' not in x:
                errors.append('Skipping ' + str(x) +
                              ' no ensembl in query result: ' + str(x))
                logger.error(errors[-1])
                continue

            ensembl_id = None

            if len(x['ensembl']) > 1:
                for g in x['ensembl']:
                    if g['gene'] in g_antibody_dict:
                        ensembl_id = g['gene']
                        break
                ensemblstr += ';'.join([g['gene'] for g in x['ensembl']])
            else:
                ensemblstr += x['ensembl']['gene']
                ensembl_id = x['ensembl']['gene']

            filename_str = ','.join(list(g_filename_dict[ensembl_id]))
            antibody_str = ','.join(list(g_antibody_dict[ensembl_id]))

            ambiguous_str = ''
            if x['symbol'] in ambiguous_gene_dict:
                ambiguous_str = ambiguous_gene_dict[x['symbol']]

            gene_node_attrs[x['query']] = {'name': x['symbol'],
                                           'represents': ensemblstr,
                                           'ambiguous': ambiguous_str,
                                           'antibody': antibody_str,
                                           'filename': filename_str}

        return gene_node_attrs, errors
