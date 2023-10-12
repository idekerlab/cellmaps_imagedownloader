#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import os
import unittest
import tempfile
import shutil
import requests
import requests_mock
from unittest.mock import MagicMock
import json
from cellmaps_utils import constants
import cellmaps_imagedownloader
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError
from cellmaps_imagedownloader.runner import CellmapsImageDownloader
from cellmaps_imagedownloader import runner


class TestCellmapsdownloaderrunner(unittest.TestCase):
    """Tests for `cellmaps_imagedownloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_constructor(self):
        """Tests constructor"""
        myobj = CellmapsImageDownloader(outdir='foo')
        self.assertIsNotNone(myobj)

    def test_constructor_no_outdir(self):
        try:
            CellmapsImageDownloader()
            self.fail('Expected Exception')
        except CellMapsImageDownloaderError as e:
            self.assertEqual('outdir is None', str(e))

    def test_register_software(self):
        prov = MagicMock()
        prov.register_software = MagicMock(return_value='12345')
        prov_json = {'keywords': ['hi'],
                     'description': 'some desc'}
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_json,
                                        provenance_utils=prov)
        myobj._register_software()
        self.assertEqual('12345', myobj._softwareid)
        self.assertEqual(1, prov.register_software.call_count)
        self.assertEqual('/foo',
                         prov.register_software.call_args_list[0][0][0])
        self.assertEqual(cellmaps_imagedownloader.__name__,
                         prov.register_software.call_args_list[0][1]['name'])
        self.assertEqual(cellmaps_imagedownloader.__version__,
                         prov.register_software.call_args_list[0][1]['version'])
        self.assertEqual(cellmaps_imagedownloader.__author__,
                         prov.register_software.call_args_list[0][1]['author'])
        self.assertEqual('some desc ' + cellmaps_imagedownloader.__description__,
                         prov.register_software.call_args_list[0][1]['description'])
        self.assertEqual(['hi', 'tools', cellmaps_imagedownloader.__name__],
                         prov.register_software.call_args_list[0][1]['keywords'])
        self.assertEqual(cellmaps_imagedownloader.__repo_url__,
                         prov.register_software.call_args_list[0][1]['url'])
        self.assertEqual('py',
                         prov.register_software.call_args_list[0][1]['file_format'])

    def test_register_image_gene_node_attrs(self):
        prov = MagicMock()
        prov.get_default_date_format_str = MagicMock(return_value='%Y-%m-%d')
        prov.register_dataset = MagicMock()
        prov.register_dataset.side_effect =['1', '2']
        prov_json = {'keywords': ['hi'],
                     'description': 'some desc'}
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_json,
                                        provenance_utils=prov)
        myobj._register_image_gene_node_attrs()
        myobj._register_image_gene_node_attrs(fold=2)
        self.assertEqual(['1', '2'], myobj._image_gene_attrid)
        self.assertEqual(2, prov.register_dataset.call_count)
        self.assertEqual('/foo/1_image_gene_node_attributes.tsv',
                         prov.register_dataset.call_args_list[0][1]['source_file'])
        self.assertEqual('/foo/2_image_gene_node_attributes.tsv',
                         prov.register_dataset.call_args_list[1][1]['source_file'])

        # todo verify all fields in data_dict are correct

    def test_add_dataset_to_crate(self):
        prov = MagicMock()
        prov.register_dataset = MagicMock(return_value='1')
        prov_json = {'keywords': ['hi'],
                     'description': 'some desc'}
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_json,
                                        provenance_utils=prov)
        self.assertEqual('1',
                         myobj._add_dataset_to_crate(data_dict={},
                                                     source_file='/srcfile'))

    def test_register_computation_no_datasets(self):
        prov = MagicMock()
        prov.register_computation = MagicMock(return_value='1')
        prov.get_login = MagicMock(return_value='smith')
        prov_json = {'keywords': ['hi'],
                     'description': 'some desc'}
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_json,
                                        provenance_utils=prov)
        myobj._softwareid = 'softwareid'
        myobj._unique_datasetid = 'uniqueid'
        myobj._samples_datasetid = 'samples'
        myobj._register_computation()
        self.assertEqual(1, prov.register_computation.call_count)
        self.assertEqual('/foo',
                         prov.register_computation.call_args_list[0][0][0])
        self.assertEqual('IF Image Loader',
                         prov.register_computation.call_args_list[0][1]['name'])
        self.assertEqual('smith',
                         prov.register_computation.call_args_list[0][1]['run_by'])
        self.assertEqual('None',
                         prov.register_computation.call_args_list[0][1]['command'])
        self.assertEqual('some desc run of ' + cellmaps_imagedownloader.__name__,
                         prov.register_computation.call_args_list[0][1]['description'])
        self.assertEqual(['hi', 'computation', 'download'],
                         prov.register_computation.call_args_list[0][1]['keywords'])
        self.assertEqual(['softwareid'],
                         prov.register_computation.call_args_list[0][1]['used_software'])
        self.assertEqual(['uniqueid', 'samples'],
                         prov.register_computation.call_args_list[0][1]['used_dataset'])
        self.assertEqual([],
                         prov.register_computation.call_args_list[0][1]['generated'])

    def test_register_computation_no_one_gene_attr_and_one_image(self):
        prov = MagicMock()
        prov.register_computation = MagicMock(return_value='1')
        prov.get_login = MagicMock(return_value='smith')
        prov_json = {'keywords': ['hi'],
                     'description': 'some desc'}
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_json,
                                        provenance_utils=prov)
        myobj._image_gene_attrid = ['image1']
        myobj._image_dataset_ids = ['1']
        myobj._softwareid = 'softwareid'
        myobj._unique_datasetid = 'uniqueid'
        myobj._samples_datasetid = 'samples'
        myobj._register_computation()
        self.assertEqual(1, prov.register_computation.call_count)
        self.assertEqual('/foo',
                         prov.register_computation.call_args_list[0][0][0])
        self.assertEqual('IF Image Loader',
                         prov.register_computation.call_args_list[0][1]['name'])
        self.assertEqual('smith',
                         prov.register_computation.call_args_list[0][1]['run_by'])
        self.assertEqual('None',
                         prov.register_computation.call_args_list[0][1]['command'])
        self.assertEqual('some desc run of ' + cellmaps_imagedownloader.__name__,
                         prov.register_computation.call_args_list[0][1]['description'])
        self.assertEqual(['hi', 'computation', 'download'],
                         prov.register_computation.call_args_list[0][1]['keywords'])
        self.assertEqual(['softwareid'],
                         prov.register_computation.call_args_list[0][1]['used_software'])
        self.assertEqual(['uniqueid', 'samples'],
                         prov.register_computation.call_args_list[0][1]['used_dataset'])
        self.assertEqual(['image1', '1'],
                         prov.register_computation.call_args_list[0][1]['generated'])

    def test_register_computation_no_two_gene_attr_and_2001_images(self):
        prov = MagicMock()
        prov.register_computation = MagicMock(return_value='1')
        prov.get_login = MagicMock(return_value='smith')
        prov_json = {'keywords': ['hi'],
                     'description': 'some desc'}
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_json,
                                        provenance_utils=prov)
        myobj._image_gene_attrid = ['image1', 'image2']
        myobj._image_dataset_ids = []
        for x in range(2010):
            myobj._image_dataset_ids.append(str(x))
        myobj._softwareid = 'softwareid'
        myobj._unique_datasetid = 'uniqueid'
        myobj._samples_datasetid = 'samples'
        myobj._register_computation()
        self.assertEqual(1, prov.register_computation.call_count)
        self.assertEqual('/foo',
                         prov.register_computation.call_args_list[0][0][0])
        self.assertEqual('IF Image Loader',
                         prov.register_computation.call_args_list[0][1]['name'])
        self.assertEqual('smith',
                         prov.register_computation.call_args_list[0][1]['run_by'])
        self.assertEqual('None',
                         prov.register_computation.call_args_list[0][1]['command'])
        self.assertEqual('some desc run of ' + cellmaps_imagedownloader.__name__,
                         prov.register_computation.call_args_list[0][1]['description'])
        self.assertEqual(['hi', 'computation', 'download'],
                         prov.register_computation.call_args_list[0][1]['keywords'])

        self.assertEqual(['softwareid'],
                         prov.register_computation.call_args_list[0][1]['used_software'])
        self.assertEqual(['uniqueid', 'samples'],
                         prov.register_computation.call_args_list[0][1]['used_dataset'])
        gen_list = prov.register_computation.call_args_list[0][1]['generated']
        self.assertEqual('image1', gen_list[0])
        self.assertEqual('image2', gen_list[1])
        self.assertEqual('0', gen_list[2])
        self.assertEqual('1999', gen_list[2001])

    def test_create_run_crate_missing_value_in_provenance(self):
        prov = MagicMock()
        prov.register_rocrate = MagicMock()
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance_utils=prov,
                                        provenance={})
        try:
            myobj._create_run_crate()
            self.fail('expected exception')
        except CellMapsImageDownloaderError as e:
            self.assertEqual('Key missing in provenance: \'name\'', str(e))

    def test_create_run_crate_invalidprovenancetype(self):
        prov = MagicMock()
        prov.register_rocrate = MagicMock()
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance_utils=prov,
                                        provenance='')
        try:
            myobj._create_run_crate()
            self.fail('expected exception')
        except CellMapsImageDownloaderError as e:
            self.assertTrue('Invalid provenance: ' in str(e))

    def test_create_run_crate(self):
        prov = MagicMock()
        prov.register_rocrate = MagicMock()

        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance_utils=prov,
                                        provenance={'name': 'foo',
                                                    'organization-name': 'icorp',
                                                    'project-name': 'myproj',
                                                    'cell-line': 'U2OS',
                                                    'release': '0.1 alpha',
                                                    'treatment': 'untreated'})
        myobj._update_provenance_with_keywords()
        myobj._update_provenance_with_description()
        myobj._create_run_crate()
        self.assertEqual(1, prov.register_rocrate.call_count)
        self.assertEqual('/foo',
                         prov.register_rocrate.call_args_list[0][0][0])
        self.assertEqual('foo',
                         prov.register_rocrate.call_args_list[0][1]['name'])
        self.assertEqual('icorp',
                         prov.register_rocrate.call_args_list[0][1]['organization_name'])
        self.assertEqual('myproj',
                         prov.register_rocrate.call_args_list[0][1]['project_name'])
        self.assertEqual('icorp myproj 0.1 alpha U2OS untreated foo IF microscopy images',
                         prov.register_rocrate.call_args_list[0][1]['description'])
        self.assertEqual(['icorp', 'myproj', '0.1 alpha', 'U2OS', 'untreated',
                          'foo', 'IF microscopy', 'images'],
                         prov.register_rocrate.call_args_list[0][1]['keywords'])

    def test_register_samples_dataset_guid_already_set(self):
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance={CellmapsImageDownloader.SAMPLES_FILEKEY: {'guid': '1'}})
        myobj._register_samples_dataset()
        self.assertEqual('1', myobj._samples_datasetid)

    def test_register_unique_dataset_guid_already_set(self):
        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance={CellmapsImageDownloader.UNIQUE_FILEKEY: {'guid': '2'}})
        myobj._register_unique_dataset()
        self.assertEqual('2', myobj._unique_datasetid)

    def test_register_samples_dataset_no_input_data_dict(self):
        prov_dict = {CellmapsImageDownloader.SAMPLES_FILEKEY: {'name': 'x'}}

        imagegen = MagicMock()
        imagegen.write_samples_to_csvfile = MagicMock()

        myobj = CellmapsImageDownloader(outdir='/foo',
                                        imagegen=imagegen,
                                        provenance=prov_dict)
        myobj._add_dataset_to_crate = MagicMock()
        myobj._add_dataset_to_crate.side_effect = ['id1']
        myobj._register_samples_dataset()
        self.assertEqual('id1', myobj._samples_datasetid)
        self.assertEqual(1, imagegen.write_samples_to_csvfile.call_count)
        self.assertEqual('/foo/samplescopy.csv',
                         imagegen.write_samples_to_csvfile.call_args_list[0][1]['csvfile'])
        self.assertEqual(1, myobj._add_dataset_to_crate.call_count)
        self.assertEqual({'name': 'x'},
                         myobj._add_dataset_to_crate.call_args_list[0][1]['data_dict'])
        self.assertEqual('/foo/samplescopy.csv',
                         myobj._add_dataset_to_crate.call_args_list[0][1]['source_file'])
        self.assertEqual(True,
                         myobj._add_dataset_to_crate.call_args_list[0][1]['skip_copy'])

    def test_register_unique_dataset_no_input_data_dict(self):
        prov_dict = {CellmapsImageDownloader.UNIQUE_FILEKEY: {'name': 'x'}}

        imagegen = MagicMock()
        imagegen.write_unique_list_to_csvfile = MagicMock()

        myobj = CellmapsImageDownloader(outdir='/foo',
                                        imagegen=imagegen,
                                        provenance=prov_dict)
        myobj._add_dataset_to_crate = MagicMock()
        myobj._add_dataset_to_crate.side_effect = ['id1']
        myobj._register_unique_dataset()
        self.assertEqual('id1', myobj._unique_datasetid)
        self.assertEqual(1, imagegen.write_unique_list_to_csvfile.call_count)
        self.assertEqual('/foo/uniquecopy.csv',
                         imagegen.write_unique_list_to_csvfile.call_args_list[0][1]['csvfile'])
        self.assertEqual(1, myobj._add_dataset_to_crate.call_count)
        self.assertEqual({'name': 'x'},
                         myobj._add_dataset_to_crate.call_args_list[0][1]['data_dict'])
        self.assertEqual('/foo/uniquecopy.csv',
                         myobj._add_dataset_to_crate.call_args_list[0][1]['source_file'])
        self.assertEqual(True,
                         myobj._add_dataset_to_crate.call_args_list[0][1]['skip_copy'])

    def test_register_samples_dataset(self):
        prov_dict = {CellmapsImageDownloader.SAMPLES_FILEKEY: {'name': 'x'}}

        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_dict,
                                        input_data_dict={CellmapsImageDownloader.SAMPLES_FILEKEY: '/x/samples.csv'})
        myobj._add_dataset_to_crate = MagicMock()
        myobj._add_dataset_to_crate.side_effect = ['id1']
        myobj._register_samples_dataset()
        self.assertEqual('id1', myobj._samples_datasetid)

        self.assertEqual(1, myobj._add_dataset_to_crate.call_count)
        self.assertEqual({'name': 'x'},
                         myobj._add_dataset_to_crate.call_args_list[0][1]['data_dict'])
        self.assertEqual('/x/samples.csv',
                         myobj._add_dataset_to_crate.call_args_list[0][1]['source_file'])
        self.assertEqual(False,
                         myobj._add_dataset_to_crate.call_args_list[0][1]['skip_copy'])

    def test_register_unique_dataset(self):
        prov_dict = {CellmapsImageDownloader.UNIQUE_FILEKEY: {'name': 'x'}}

        myobj = CellmapsImageDownloader(outdir='/foo',
                                        provenance=prov_dict,
                                        input_data_dict={CellmapsImageDownloader.UNIQUE_FILEKEY: '/x/unique.csv'})
        myobj._add_dataset_to_crate = MagicMock()
        myobj._add_dataset_to_crate.side_effect = ['id1']
        myobj._register_unique_dataset()
        self.assertEqual('id1', myobj._unique_datasetid)

        self.assertEqual(1, myobj._add_dataset_to_crate.call_count)
        self.assertEqual({'name': 'x'},
                         myobj._add_dataset_to_crate.call_args_list[0][1]['data_dict'])
        self.assertEqual('/x/unique.csv',
                         myobj._add_dataset_to_crate.call_args_list[0][1]['source_file'])
        self.assertEqual(False,
                         myobj._add_dataset_to_crate.call_args_list[0][1]['skip_copy'])

    def test_run(self):
        """ Tests run()"""
        temp_dir = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(temp_dir, 'run')
            myobj = CellmapsImageDownloader(outdir=run_dir)
            try:
                myobj.run()
                self.fail('Expected CellMapsImageDownloaderError')
            except CellMapsImageDownloaderError as c:
                self.assertTrue('Invalid provenance' in str(c))
            self.assertFalse(os.path.isfile(os.path.join(run_dir, 'output.log')))
            self.assertFalse(os.path.isfile(os.path.join(run_dir, 'error.log')))
        finally:
            shutil.rmtree(temp_dir)

    def test_run_with_skip_logging_false(self):
        """ Tests run()"""
        temp_dir = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(temp_dir, 'run')
            myobj = CellmapsImageDownloader(outdir=run_dir,
                                            skip_logging=False)
            try:
                myobj.run()
                self.fail('Expected CellMapsImageDownloaderError')
            except CellMapsImageDownloaderError as c:
                self.assertTrue('Invalid provenance' in str(c))
            self.assertTrue(os.path.isfile(os.path.join(run_dir, 'output.log')))
            self.assertTrue(os.path.isfile(os.path.join(run_dir, 'error.log')))

        finally:
            shutil.rmtree(temp_dir)

    def test_download_file(self):
        temp_dir = tempfile.mkdtemp()

        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, status_code=200,
                      text='somedata')
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                runner.download_file((mockurl, a_dest_file))
            self.assertTrue(os.path.isfile(a_dest_file))
            with open(a_dest_file, 'r') as f:
                data = f.read()
                self.assertEqual('somedata', data)
        finally:
            shutil.rmtree(temp_dir)

    def test_download_file_failure(self):
        temp_dir = tempfile.mkdtemp()

        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, status_code=500,
                      text='error')
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                rstatus, rtext, rtuple = runner.download_file((mockurl, a_dest_file))
            self.assertEqual(500, rstatus)
            self.assertEqual('error', rtext)
            self.assertEqual((mockurl, a_dest_file), rtuple)
            self.assertFalse(os.path.isfile(a_dest_file))

        finally:
            shutil.rmtree(temp_dir)

    def test_download_raise_httperror(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, exc=requests.exceptions.HTTPError('httperror'))
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                rstatus, rtext, rtuple = runner.download_file((mockurl, a_dest_file))
            self.assertEqual(-1, rstatus)
            self.assertEqual('httperror', rtext)
            self.assertEqual((mockurl, a_dest_file), rtuple)
            self.assertFalse(os.path.isfile(a_dest_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_download_raise_connectionerror(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, exc=requests.exceptions.ConnectionError('conerror'))
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                rstatus, rtext, rtuple = runner.download_file((mockurl, a_dest_file))
            self.assertEqual(-2, rstatus)
            self.assertEqual('conerror', rtext)
            self.assertEqual((mockurl, a_dest_file), rtuple)
            self.assertFalse(os.path.isfile(a_dest_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_download_raise_timeouterror(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, exc=requests.exceptions.Timeout('timeerror'))
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                rstatus, rtext, rtuple = runner.download_file((mockurl, a_dest_file))
            self.assertEqual(-3, rstatus)
            self.assertEqual('timeerror', rtext)
            self.assertEqual((mockurl, a_dest_file), rtuple)
            self.assertFalse(os.path.isfile(a_dest_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_download_raise_requestexceptionerror(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, exc=requests.exceptions.RequestException('error'))
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                rstatus, rtext, rtuple = runner.download_file((mockurl, a_dest_file))
            self.assertEqual(-4, rstatus)
            self.assertEqual('error', rtext)
            self.assertEqual((mockurl, a_dest_file), rtuple)
            self.assertFalse(os.path.isfile(a_dest_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_download_raise_exceptionerror(self):
        temp_dir = tempfile.mkdtemp()
        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, exc=Exception('error'))
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                rstatus, rtext, rtuple = runner.download_file((mockurl, a_dest_file))
            self.assertEqual(-5, rstatus)
            self.assertEqual('error', rtext)
            self.assertEqual((mockurl, a_dest_file), rtuple)
            self.assertFalse(os.path.isfile(a_dest_file))
        finally:
            shutil.rmtree(temp_dir)

    def test_download_file_skip_existing_empty_file_exists(self):
        temp_dir = tempfile.mkdtemp()

        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, status_code=200,
                      text='somedata')
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                open(a_dest_file, 'a').close()

                runner.download_file_skip_existing((mockurl, a_dest_file))
            self.assertTrue(os.path.isfile(a_dest_file))
            with open(a_dest_file, 'r') as f:
                data = f.read()
                self.assertEqual('somedata', data)
        finally:
            shutil.rmtree(temp_dir)

    def test_download_file_skip_existing_file_exists(self):
        temp_dir = tempfile.mkdtemp()

        try:
            mockurl = 'http://fakey.fake.com/ha.txt'

            with requests_mock.mock() as m:
                m.get(mockurl, status_code=200,
                      text='somedata')
                a_dest_file = os.path.join(temp_dir, 'downloadedfile.txt')
                with open(a_dest_file, 'w') as f:
                    f.write('blah')

                self.assertIsNone(runner.download_file_skip_existing((mockurl, a_dest_file)))
            self.assertTrue(os.path.isfile(a_dest_file))
            with open(a_dest_file, 'r') as f:
                data = f.read()
                self.assertEqual('blah', data)
        finally:
            shutil.rmtree(temp_dir)

    def test_create_output_directory(self):
        temp_dir = tempfile.mkdtemp()

        # fail if directory already exists
        try:
            crunner = CellmapsImageDownloader(outdir=temp_dir)
            crunner._create_output_directory()
            self.fail('Expected exception')
        except CellMapsImageDownloaderError as ce:
            self.assertTrue(' already exists' in str(ce))

        try:
            run_dir = os.path.join(temp_dir, 'run')
            crunner = CellmapsImageDownloader(outdir=run_dir)
            crunner._create_output_directory()
            for c in constants.COLORS:
                self.assertTrue(os.path.isdir(os.path.join(run_dir, c)))
        finally:
            shutil.rmtree(temp_dir)

    def test_write_task_start_json(self):
        temp_dir = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(temp_dir, 'run')
            crunner = CellmapsImageDownloader(outdir=run_dir)
            crunner._create_output_directory()
            crunner._write_task_start_json()
            start_file = None
            for entry in os.listdir(run_dir):
                if not entry.endswith('_start.json'):
                    continue
                start_file = os.path.join(run_dir, entry)
            self.assertIsNotNone(start_file)

            with open(start_file, 'r') as f:
                data = json.load(f)

            self.assertEqual(cellmaps_imagedownloader.__version__,
                             data['version'])
            self.assertTrue(data['start_time'] > 0)
            self.assertEqual(run_dir, data['outdir'])
        finally:
            shutil.rmtree(temp_dir)

    def test_get_color_download_map(self):
        temp_dir = tempfile.mkdtemp()
        try:
            run_dir = os.path.abspath(os.path.join(temp_dir, 'run'))
            crunner = CellmapsImageDownloader(outdir=run_dir)
            res = crunner._get_color_download_map()
            self.assertEqual(4, len(res))
            for c in constants.COLORS:
                self.assertTrue(os.path.join(run_dir, c) in res[c])
        finally:
            shutil.rmtree(temp_dir)

    def test_get_download_tuples(self):
        temp_dir = tempfile.mkdtemp()
        try:
            imageurlgen = MagicMock()
            imageurlgen.get_next_image_url = MagicMock()

            def fake_gen(color_map):
                for line in [('url1', '/url1'),
                             ('url2', '/url2')]:
                    yield line

            imageurlgen.get_next_image_url.side_effect = fake_gen

            run_dir = os.path.abspath(os.path.join(temp_dir, 'run'))
            crunner = CellmapsImageDownloader(outdir=run_dir,
                                              imageurlgen=imageurlgen)
            res = crunner._get_download_tuples()
            self.assertEqual(2, len(res))

            self.assertTrue(('url1', '/url1') in res)
            self.assertTrue(('url2', '/url2') in res)
        finally:
            shutil.rmtree(temp_dir)
