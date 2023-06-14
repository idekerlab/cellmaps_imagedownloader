#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import os
import unittest
import tempfile
import shutil
import requests_mock
from unittest.mock import MagicMock
from unittest.mock import Mock
import json
from cellmaps_utils import constants
import cellmaps_imagedownloader
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError
from cellmaps_utils.exceptions import CellMapsProvenanceError
from cellmaps_imagedownloader.runner import CellmapsImageDownloader
from cellmaps_imagedownloader.runner import ImageDownloader
from cellmaps_imagedownloader.gene import ImageGeneNodeAttributeGenerator
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
