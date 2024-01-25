#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Integration Tests for `cellmaps_imagedownloader` package."""

import os
import tempfile
import shutil
import unittest
from cellmaps_imagedownloader import cellmaps_imagedownloadercmd

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


@unittest.skipUnless(os.getenv('CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST') is not None, SKIP_REASON)
class TestIntegrationCellmaps_downloader(unittest.TestCase):
    """Tests for `cellmaps_imagedownloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def get_data_dir(self):
        return os.path.join(os.path.dirname(__file__), 'data')

    def get_test_single_sample(self):
        return os.path.join(self.get_data_dir(), 'test_single_sample.csv')

    def get_test_single_unique(self):
        return os.path.join(self.get_data_dir(), 'test_single_unique.csv')

    def get_test_provenance(self):
        return os.path.join(self.get_data_dir(), 'test_provenance.json')

    def test_single_download(self):
        """Tests parse arguments"""
        temp_dir = tempfile.mkdtemp()
        try:
            run_dir = os.path.join(temp_dir, 'run')
            res = cellmaps_imagedownloadercmd.main(['myprog', run_dir,
                                                    '--poolsize', '1',
                                                    '--samples', self.get_test_single_sample(),
                                                    '--unique', self.get_test_single_unique(),
                                                    '--provenance', self.get_test_provenance(),
                                                    '--skip_logging'])
            self.assertEqual(0, res)
            self.assertTrue(os.path.isfile(os.path.join(run_dir, '1_image_gene_node_attributes.tsv')))
            self.assertTrue(os.path.isfile(os.path.join(run_dir, '2_image_gene_node_attributes.tsv')))
            for c in ['red', 'yellow', 'green', 'blue']:
                self.assertTrue(os.path.isfile(os.path.join(run_dir, c, '1_A1_1_' + c + '.jpg')))
        finally:
            shutil.rmtree(temp_dir)
