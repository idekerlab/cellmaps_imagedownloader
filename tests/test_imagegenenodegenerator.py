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
                                'sample': '1', 'status': '35',
                                'locations': 'Golgi apparatus',
                                'antibody': 'HPA000992',
                                'ensembl_ids': 'ENSG00000066455',
                                'gene_names': 'GOLGA5'},
                               {'filename': '/archive/1/1_A1_2_',
                                'if_plate_id': '1', 'position': 'A1',
                                'sample': '2', 'status': '35',
                                'locations': 'Golgi apparatus',
                                'antibody': 'HPA000992',
                                'ensembl_ids': 'ENSG00000066455',
                                'gene_names': 'GOLGA5'},
                               {'filename': '/archive/1/1_A3_1_',
                                'if_plate_id': '1', 'position': 'A3',
                                'sample': '1', 'status': '35',
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

    def test_write_samples_as_csvfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test write when samples is None
            try:
                imagegen = ImageGeneNodeAttributeGenerator()
                imagegen.write_samples_as_csvfile(os.path.join(temp_dir, 'foo.csv'))
                self.fail('Expected exception')
            except CellMapsImageDownloaderError as ce:
                self.assertEqual('samples list is None', str(ce))

                # test with empty samples
                imagegen = ImageGeneNodeAttributeGenerator(samples_list=[])
                out_file = os.path.join(temp_dir, 'emptyfile.csv')
                imagegen.write_samples_as_csvfile(outfile=os.path.join(temp_dir, out_file))

                with open(out_file, 'r') as f:
                    data = f.read()
                    self.assertEqual(','.join(ImageGeneNodeAttributeGenerator.SAMPLES_HEADER_COLS) + '\n', data)

                # test with a row
                sample = {}
                for key in ImageGeneNodeAttributeGenerator.SAMPLES_HEADER_COLS:
                    sample[key] = key + '_' + str(len(key))

                imagegen = ImageGeneNodeAttributeGenerator(samples_list=[sample])
                one_line_file = os.path.join(temp_dir, 'oneline.csv')
                imagegen.write_samples_as_csvfile(outfile=os.path.join(temp_dir, one_line_file))

                with open(one_line_file, 'r') as f:
                    data = f.readline()
                    self.assertEqual(','.join(ImageGeneNodeAttributeGenerator.SAMPLES_HEADER_COLS) + '\n', data)
                    data = f.readline()
                    self.assertEqual('filename_8,if_plate_id_11,position_8,'
                                     'sample_6,status_6,locations_9,'
                                     'antibody_8,ensembl_ids_11,'
                                     'gene_names_10\n', data)
        finally:
            shutil.rmtree(temp_dir)

    def test_write_unique_list_as_csvfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # test write when samples is None
            try:
                imagegen = ImageGeneNodeAttributeGenerator()
                imagegen.write_unique_list_as_csvfile(os.path.join(temp_dir, 'foo.csv'))
                self.fail('Expected exception')
            except CellMapsImageDownloaderError as ce:
                self.assertEqual('unique list is None', str(ce))

                # test with empty unique
                imagegen = ImageGeneNodeAttributeGenerator(unique_list=[])
                out_file = os.path.join(temp_dir, 'emptyfile.csv')
                imagegen.write_unique_list_as_csvfile(outfile=os.path.join(temp_dir, out_file))

                with open(out_file, 'r') as f:
                    data = f.read()
                    self.assertEqual(','.join(ImageGeneNodeAttributeGenerator.UNIQUE_HEADER_COLS) + '\n', data)

                # test with a row
                unique = {}
                for key in ImageGeneNodeAttributeGenerator.UNIQUE_HEADER_COLS:
                    unique[key] = key + '_' + str(len(key))

                imagegen = ImageGeneNodeAttributeGenerator(unique_list=[unique])
                one_line_file = os.path.join(temp_dir, 'oneline.csv')
                imagegen.write_unique_list_as_csvfile(outfile=os.path.join(temp_dir, one_line_file))

                with open(one_line_file, 'r') as f:
                    data = f.readline()
                    self.assertEqual(','.join(ImageGeneNodeAttributeGenerator.UNIQUE_HEADER_COLS) + '\n', data)
                    data = f.readline()
                    self.assertEqual('antibody_8,ensembl_ids_11,'
                                     'gene_names_10,atlas_name_10,'
                                     'locations_9,n_location_10\n', data)
        finally:
            shutil.rmtree(temp_dir)

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
        antibody_dict, filename_dict = imagegen.get_dicts_of_gene_to_antibody_filename()
        self.assertEqual({}, antibody_dict)
        self.assertEqual({}, filename_dict)

        # test two samples
        samples = [{'ensembl_ids': 'ensemble_one',
                    'antibody': 'antibody_one',
                    'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '2'},
                   {'ensembl_ids': 'ensemble_two',
                    'antibody': 'antibody_two',
                    'if_plate_id': '3',
                    'position': 'B1',
                    'sample': '4'}]
        imagegen = ImageGeneNodeAttributeGenerator(samples_list=samples)
        antibody_dict, filename_dict = imagegen.get_dicts_of_gene_to_antibody_filename()
        self.assertEqual({'ensemble_one': {'antibody_one'},
                          'ensemble_two': {'antibody_two'}}, antibody_dict)
        self.assertEqual({'ensemble_one': {'1_A1_2_'},
                          'ensemble_two': {'3_B1_4_'}}, filename_dict)

        # run again this time limit two antibody_two
        antibody_dict, filename_dict = imagegen.get_dicts_of_gene_to_antibody_filename(allowed_antibodies={'antibody_two'})
        self.assertEqual({'ensemble_two': {'antibody_two'}}, antibody_dict)
        self.assertEqual({'ensemble_two': {'3_B1_4_'}}, filename_dict)
