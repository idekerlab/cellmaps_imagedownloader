#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import os
import unittest
import tempfile
import shutil
import json
from unittest.mock import patch, mock_open
from unittest.mock import MagicMock
import cellmaps_imagedownloader
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError
from cellmaps_imagedownloader.gene import GeneQuery
from cellmaps_imagedownloader.gene import ImageGeneNodeAttributeGenerator


class TestImageGeneNodeAttributeGenerator(unittest.TestCase):
    """Tests for `cellmaps_imagedownloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_unique_list(self):
        gen = ImageGeneNodeAttributeGenerator(unique_list='foo')
        self.assertEqual('foo', gen.get_unique_list())

    def test_filter_samples_by_unique_list_unique_list_is_none(self):
        gen = ImageGeneNodeAttributeGenerator(samples_list=[])
        gen._filter_samples_by_unique_list()

    def test_filter_samples_by_unique_list_samples_list_is_none(self):
        gen = ImageGeneNodeAttributeGenerator(unique_list=[])
        gen._filter_samples_by_unique_list()

    def test_filter_samples_by_unique_list(self):
        gen = ImageGeneNodeAttributeGenerator(samples_list=[{'antibody': 'xxx'},
                                                            {'antibody': 'HPA123'}],
                                              unique_list=[{'no_antibody': 1},
                                                           {'antibody': 'HPA123'},
                                                           {'antibody': 'CAB456'}])
        self.assertEqual(2, len(gen.get_samples_list()))
        gen._filter_samples_by_unique_list()

        self.assertEqual([{'antibody': 'HPA123'}], gen.get_samples_list())

    def test_get_image_id_for_sample_none(self):
        try:
            ImageGeneNodeAttributeGenerator.get_image_id_for_sample(None)
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('sample is None', str(ce))

    def test_get_image_id_for_sample_not_dict(self):
        try:
            ImageGeneNodeAttributeGenerator.get_image_id_for_sample('foo')
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('sample is not a dict', str(ce))

    def test_get_image_id_for_sample_invalidsample(self):
        for key_to_del in ['antibody', 'position', 'sample', 'if_plate_id']:
            sample = {'antibody': 'a',
                      'position': 'b',
                      'sample': 'c',
                      'if_plate_id': 'd'}
            del sample[key_to_del]
            try:
                ImageGeneNodeAttributeGenerator.get_image_id_for_sample(sample)
                self.fail('Expected exception')
            except CellMapsImageDownloaderError as ce:
                self.assertEqual(key_to_del + ' not in sample', str(ce))

    def test_get_image_id_for_sample_valid(self):
        sample = {'antibody': 'HPA012345',
                  'position': 2,
                  'sample': 'C1',
                  'if_plate_id': 1}
        res = ImageGeneNodeAttributeGenerator.get_image_id_for_sample(sample)
        self.assertEqual('12345/1_2_C1_', res)

    def test_get_samples_list_image_ids_samples_list_none(self):
        gen = ImageGeneNodeAttributeGenerator(samples_list=None)
        try:
            gen.get_samples_list_image_ids()
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('samples list is None', str(ce))

    def test_get_samples_list_image_ids(self):
        samples = [{'antibody': 'HPA012345',
                    'position': 2,
                    'sample': 'C1',
                    'if_plate_id': 1},
                   {'antibody': 'CAB000123',
                    'position': 'A',
                    'sample': 'A1',
                    'if_plate_id': '2'}
                   ]
        gen = ImageGeneNodeAttributeGenerator(samples_list=samples)
        res = gen.get_samples_list_image_ids()
        self.assertEqual(2, len(res))
        self.assertTrue('12345/1_2_C1_' in res)
        self.assertTrue('123/2_A_A1_' in res)

    def test_filter_samples_by_sample_urlmap_samples_none(self):
        gen = ImageGeneNodeAttributeGenerator(samples_list=None)
        try:
            gen.filter_samples_by_sample_urlmap(sample_url_map={})
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('samples list is None', str(ce))

    def test_filter_samples_by_sample_urlmap_sample_url_map_none(self):
        gen = ImageGeneNodeAttributeGenerator(samples_list=[])
        gen.filter_samples_by_sample_urlmap(None)

    def test_filter_samples_by_sample_urlmap(self):
        samples = [{'antibody': 'HPA012345',
                    'position': 2,
                    'sample': 'C1',
                    'if_plate_id': 1},
                   {'antibody': 'CAB000123',
                    'position': 'A',
                    'sample': 'A1',
                    'if_plate_id': '2'}
                   ]
        gen = ImageGeneNodeAttributeGenerator(samples_list=samples)
        gen.filter_samples_by_sample_urlmap(None)
        self.assertEqual(samples, gen.get_samples_list())
        gen.filter_samples_by_sample_urlmap({'12345/1_2_C1_': 'http://foo'})
        self.assertEqual([samples[0]], gen.get_samples_list())

    def test_get_samples_from_csv_file(self):
        try:
            # Test case when csvfile is None
            result = ImageGeneNodeAttributeGenerator.get_samples_from_csvfile()
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('csvfile is None', str(ce))

        # Test case when file is not found
        with self.assertRaises(FileNotFoundError):
            ImageGeneNodeAttributeGenerator.get_samples_from_csvfile('non_existent_file.csv')

        # Test case when csvfile is empty
        with patch('builtins.open', mock_open(read_data='')):
            result = ImageGeneNodeAttributeGenerator.get_samples_from_csvfile('test.csv')
            self.assertEqual(result, [])

        # Test case when csvfile has data
        csv_data = 'filename,if_plate_id,position,sample,' \
                   'status,locations,antibody,ensembl_ids,' \
                   'gene_names\n/archive/1/1_A1_1_,1,A1,1,35,' \
                   'Golgi apparatus,HPA000992,ENSG00000066455,' \
                   'GOLGA5\n/archive/1/1_A1_2_,1,A1,2,35,' \
                   'Golgi apparatus,HPA000992,ENSG00000066455,' \
                   'GOLGA5\n/archive/1/1_A3_1_,1,A3,1,35,' \
                   'Cytosol,Nucleoplasm,HPA002899,' \
                   'ENSG00000183092,BEGAIN\n'
        with patch('builtins.open', mock_open(read_data=csv_data)):
            result = ImageGeneNodeAttributeGenerator.get_samples_from_csvfile('test.csv')
            expected_result = [{'filename': '/archive/1/1_A1_1_',
                                'if_plate_id': '1', 'position': 'A1',
                                'sample': '1',
                                'locations': 'Golgi apparatus',
                                'antibody': 'HPA000992',
                                'ensembl_ids': 'ENSG00000066455',
                                'gene_names': 'GOLGA5'},
                               {'filename': '/archive/1/1_A1_2_',
                                'if_plate_id': '1', 'position': 'A1',
                                'sample': '2',
                                'locations': 'Golgi apparatus',
                                'antibody': 'HPA000992',
                                'ensembl_ids': 'ENSG00000066455',
                                'gene_names': 'GOLGA5'},
                               {'filename': '/archive/1/1_A3_1_',
                                'if_plate_id': '1', 'position': 'A3',
                                'sample': '1',
                                'locations': 'Cytosol',
                                'antibody': 'Nucleoplasm',
                                'ensembl_ids': 'HPA002899',
                                'gene_names': 'ENSG00000183092'}]

            self.assertEqual(result, expected_result)

    def test_get_image_antibodies_from_csvfile(self):
        try:
            # Test case when csvfile is None
            result = ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile()
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('csvfile is None', str(ce))

        # Test case when file is not found
        with self.assertRaises(FileNotFoundError):
            ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile('non_existent_file.csv')

        # Test case when csvfile is empty
        with patch('builtins.open', mock_open(read_data='')):
            result = ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile('test.csv')
            self.assertEqual(result, [])

        # Test case when csvfile has data
        csv_data = 'antibody,ensembl_ids,gene_names,atlas_name,' \
                   'locations,n_location\nABC,ENSG00000171921,CDK5,' \
                   'Brain,CA1,1\nDEF,ENSG00000173672,GAD2,Brain,CA2,' \
                   '2\nGHI,ENSG00000172137,DLG4,Brain,CA3,3\n'
        with patch('builtins.open', mock_open(read_data=csv_data)):
            result = ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile('test.csv')
            expected_result = [{'antibody': 'ABC',
                                'ensembl_ids': 'ENSG00000171921',
                                'gene_names': 'CDK5', 'atlas_name': 'Brain',
                                'locations': 'CA1', 'n_location': '1'},
                               {'antibody': 'DEF',
                                'ensembl_ids': 'ENSG00000173672',
                                'gene_names': 'GAD2', 'atlas_name': 'Brain',
                                'locations': 'CA2', 'n_location': '2'},
                               {'antibody': 'GHI',
                                'ensembl_ids': 'ENSG00000172137',
                                'gene_names': 'DLG4', 'atlas_name': 'Brain',
                                'locations': 'CA3', 'n_location': '3'}]
            self.assertEqual(result, expected_result)

    def test_get_set_of_antibodies_from_unique_list(self):
        # try with None for unique list
        try:
            imagegen = ImageGeneNodeAttributeGenerator()
            imagegen._get_set_of_antibodies_from_unique_list()
            self.fail('Expected Exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('unique list is None', str(ce))

        # try with empty unique list
        imagegen = ImageGeneNodeAttributeGenerator(unique_list=[])
        self.assertEqual(0, len(imagegen._get_set_of_antibodies_from_unique_list()))

        # try with two unique entries and one entry that is invalid
        unique_list = [{'antibody': 'one'}, {'foo': 'invalid'},
                       {'antibody': 'one'},
                       {'antibody': 'two'}]
        imagegen = ImageGeneNodeAttributeGenerator(unique_list=unique_list)
        res = imagegen._get_set_of_antibodies_from_unique_list()
        self.assertEqual(2, len(res))
        self.assertTrue('one' in res)
        self.assertTrue('two' in res)

    def test_get_dicts_of_gene_to_antibody_filename(self):
        imagegen = ImageGeneNodeAttributeGenerator()

        # test where samples list is None
        try:
            imagegen.get_dicts_of_gene_to_antibody_filename()
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('samples list is None', str(ce))

        # test where samples list is empty
        imagegen = ImageGeneNodeAttributeGenerator(samples_list=[])
        antibody_dict, filename_dict, ambig_dict = imagegen.get_dicts_of_gene_to_antibody_filename()
        self.assertEqual({}, antibody_dict)
        self.assertEqual({}, filename_dict)
        self.assertEqual({}, ambig_dict)

        # test two samples
        samples = [{'ensembl_ids': 'ensemble_one',
                    'antibody': 'antibody_one',
                    'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '2'},
                   {'ensembl_ids': 'ensemble_two,ensemble_three',
                    'antibody': 'antibody_two',
                    'if_plate_id': '3',
                    'position': 'B1',
                    'sample': '4'}]
        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples)
        antibody_dict, filename_dict, ambig_dict = imagegen.get_dicts_of_gene_to_antibody_filename()
        self.assertEqual(3, len(antibody_dict.keys()))

        self.assertEqual('antibody_one', antibody_dict['ensemble_one'])
        self.assertEqual('antibody_two', antibody_dict['ensemble_two'])
        self.assertEqual('antibody_two', antibody_dict['ensemble_three'])

        self.assertEqual({'antibody_one': {'1_A1_2_'},
                          'antibody_two': {'3_B1_4_'}}, filename_dict)
        self.assertEqual({'antibody_two': ['ensemble_two', 'ensemble_three']}, ambig_dict)

    def test_get_unique_ids_from_samplelist(self):

        samples = [{'ensembl_ids': 'ENSG01,ENSG02'}]

        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples)
        id_set = imagegen._get_unique_ids_from_samplelist()
        self.assertTrue(2, len(id_set))
        self.assertTrue('ENSG02' in id_set)
        self.assertTrue('ENSG01' in id_set)

    def test_get_gene_node_attributes_simple(self):
        samples = [{'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1',
                    'antibody': 'HPA000992',
                    'ensembl_ids': 'ENSG00000066455',
                    'gene_names': 'GOLGA5'},
                   {'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '2',
                    'antibody': 'HPA000992',
                    'ensembl_ids': 'ENSG00000066455',
                    'gene_names': 'GOLGA5'},
                   {'if_plate_id': '1',
                    'position': 'A3',
                    'sample': '1',
                    'antibody': 'HPA002899',
                    'ensembl_ids': 'ENSG00000183092',
                    'gene_names': 'BEGAIN'},
                   {'if_plate_id': '1',
                    'position': 'A3',
                    'sample': '2',
                    'antibody': 'HPA002899',
                    'ensembl_ids': 'ENSG00000183092',
                    'gene_names': 'BEGAIN'}]
        unique = [{'antibody': 'HPA000992',
                   'ensembl_ids': 'ENSG00000066455',
                   'gene_names': 'GOLGA5'},
                  {'antibody': 'HPA002899',
                   'ensembl_ids': 'ENSG00000183092',
                   'gene_names': 'BEGAIN'}
                  ]

        mockgenequery = MagicMock()
        mockgenequery.get_symbols_for_genes = MagicMock(return_value=[{'query': 'ENSG00000066455',
                                                                       '_id': '9950',
                                                                       '_score': 25.046944,
                                                                       'ensembl': {'gene':'ENSG00000066455'},
                                                                       'symbol': 'GOLGA5'
                                                                       },
                                                                       {'query': 'ENSG00000183092',
                                                                        '_id': '57596', '_score': 24.51739,
                                                                        'ensembl': {'gene': 'ENSG00000183092'},
                                                                        'symbol': 'BEGAIN'}])

        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples,
                                                   unique_list=unique,
                                                   genequery=mockgenequery)

        res = imagegen.get_gene_node_attributes()
        gene_node_attrs = res[0]
        print(gene_node_attrs)
        print('\n\n\n\n')
        self.assertTrue('GOLGA5' in gene_node_attrs)
        self.assertEqual('GOLGA5', gene_node_attrs['GOLGA5']['name'])
        self.assertEqual('ENSG00000066455', gene_node_attrs['GOLGA5']['represents'])
        self.assertTrue(gene_node_attrs['GOLGA5']['filename'] == '1_A1_2_' or
                        gene_node_attrs['GOLGA5']['filename'] == '1_A1_1_')

        self.assertTrue('GOLGA5' in gene_node_attrs)

        # check errors is empty
        self.assertEqual([], res[1])

    def test_get_gene_node_attributes_no_symbol(self):
        samples = [{'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1',
                    'antibody': 'HPA000992',
                    'ensembl_ids': 'ENSG1',
                    'gene_names': 'XXX'}]
        unique = [{'antibody': 'HPA000992',
                   'ensembl_ids': 'ENSG1',
                   'gene_names': 'XXX'}]

        mockgenequery = MagicMock()
        mockgenequery.get_symbols_for_genes = MagicMock(return_value=[{'query': 'ENSG1',
                                                                       '_id': '9950',
                                                                       '_score': 25.046944,
                                                                       'ensembl': {'gene':'ENSG1'},
                                                                       }])

        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples,
                                                   unique_list=unique,
                                                   genequery=mockgenequery)

        res = imagegen.get_gene_node_attributes()
        gene_node_attrs = res[0]

        self.assertTrue('ENSG1' in gene_node_attrs)
        self.assertEqual('ENSG1', gene_node_attrs['ENSG1']['name'])
        self.assertEqual('ENSG1', gene_node_attrs['ENSG1']['represents'])
        self.assertTrue(gene_node_attrs['ENSG1']['filename'] == '1_A1_2_' or
                        gene_node_attrs['ENSG1']['filename'] == '1_A1_1_')

    def test_get_gene_node_attributes_no_ensembl(self):
        samples = [{'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1',
                    'antibody': 'HPA000992',
                    'ensembl_ids': 'ENSG1',
                    'gene_names': 'XXX'}]
        unique = [{'antibody': 'HPA000992',
                   'ensembl_ids': 'ENSG1',
                   'gene_names': 'XXX'}]

        mockgenequery = MagicMock()
        mockgenequery.get_symbols_for_genes = MagicMock(return_value=[{'query': 'ENSG1',
                                                                       'notfound': True,
                                                                       }])

        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples,
                                                   unique_list=unique,
                                                   genequery=mockgenequery)

        res = imagegen.get_gene_node_attributes()
        self.assertEqual({}, res[0])

        # check we got an error
        self.assertTrue(1, len(res[1]))
        self.assertTrue('no ensembl in query result' in res[1][0])

    def test_get_gene_node_attributes_no_ensembl(self):
        samples = [{'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1',
                    'antibody': 'HPA000992',
                    'ensembl_ids': 'ENSG1',
                    'gene_names': 'XXX'}]
        unique = [{'antibody': 'HPA000992',
                   'ensembl_ids': 'ENSG1',
                   'gene_names': 'XXX'}]

        mockgenequery = MagicMock()
        mockgenequery.get_symbols_for_genes = MagicMock(return_value=[{'query': 'ENSG1',
                                                                       'symbol': 'XXX',
                                                                       }])

        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples,
                                                   unique_list=unique,
                                                   genequery=mockgenequery)

        res = imagegen.get_gene_node_attributes()
        self.assertEqual({}, res[0])

        # check we got an error
        self.assertTrue(1, len(res[1]))
        self.assertTrue('no ensembl in query result' in res[1][0])

    def test_get_gene_node_attributes_multiple_ensembl(self):
        samples = [{'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1',
                    'antibody': 'HPA000992',
                    'ensembl_ids': 'ENSG1,ENSG2',
                    'gene_names': 'geneA'}]
        unique = [{'antibody': 'HPA000992',
                   'ensembl_ids': 'ENSG1,ENSG2',
                   'gene_names': 'geneA'}]

        mockgenequery = MagicMock()
        mockgenequery.get_symbols_for_genes = MagicMock(return_value=[{'query': 'ENSG1',
                                                                       '_id': '3975',
                                                                       '_score': 24.515644,
                                                                       'ensembl': [{'gene': 'ENSG1'},
                                                                                   {'gene': 'ENSG2'}],
                                                                       'symbol': 'geneA'}])

        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples,
                                                   unique_list=unique,
                                                   genequery=mockgenequery)

        res = imagegen.get_gene_node_attributes()
        self.assertEqual('ENSG1,ENSG2', res[0]['geneA']['represents'])

        # check we got no error
        self.assertEqual(0, len(res[1]))

    def test_write_samples_to_csvfile_with_no_file(self):
        mockgenequery = MagicMock()
        imagegen = ImageGeneNodeAttributeGenerator(genequery=mockgenequery)
        try:
            imagegen.write_samples_to_csvfile(None)
            self.fail('expected exception')
        except CellMapsImageDownloaderError as e:
            self.assertEqual('csvfile is None', str(e))

    def test_write_samples_to_csvfile(self):
        mockgenequery = MagicMock()
        temp_dir = tempfile.mkdtemp()
        try:
            csvfile = os.path.join(temp_dir, 'foo.csv')
            samples = [{'filename': '/archive/1/1_A1_1_',
                        'if_plate_id': '1', 'position': 'A1',
                        'sample': '1',
                        'locations': 'Golgi apparatus',
                        'antibody': 'HPA000992',
                        'ensembl_ids': 'ENSG00000066455',
                        'gene_names': 'GOLGA5'}]
            imagegen = ImageGeneNodeAttributeGenerator(genequery=mockgenequery,
                                                       samples_list=samples)
            imagegen.write_samples_to_csvfile(csvfile)

            res = ImageGeneNodeAttributeGenerator.get_samples_from_csvfile(csvfile)
            self.assertEqual(samples[0], res[0])
        finally:
            shutil.rmtree(temp_dir)

    def test_write_unique_list_to_csvfile_with_no_file(self):
        mockgenequery = MagicMock()
        imagegen = ImageGeneNodeAttributeGenerator(genequery=mockgenequery)
        try:
            imagegen.write_unique_list_to_csvfile(None)
            self.fail('expected exception')
        except CellMapsImageDownloaderError as e:
            self.assertEqual('csvfile is None', str(e))

    def test_write_unique_list_to_csvfile(self):
        mockgenequery = MagicMock()
        temp_dir = tempfile.mkdtemp()
        try:
            csvfile = os.path.join(temp_dir, 'foo.csv')
            unique = [{'antibody': 'ABC',
                                'ensembl_ids': 'ENSG00000171921',
                                'gene_names': 'CDK5', 'atlas_name': 'Brain',
                                'locations': 'CA1', 'n_location': '1'}]
            imagegen = ImageGeneNodeAttributeGenerator(genequery=mockgenequery,
                                                       unique_list=unique)
            imagegen.write_unique_list_to_csvfile(csvfile)

            res = ImageGeneNodeAttributeGenerator.get_unique_list_from_csvfile(csvfile)
            self.assertEqual(unique[0], res[0])
        finally:
            shutil.rmtree(temp_dir)





