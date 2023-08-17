#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `LinkPrefixImageDownloadTupleGenerator` package."""
import unittest
from unittest.mock import MagicMock

from cellmaps_imagedownloader.proteinatlas import LinkPrefixImageDownloadTupleGenerator

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestLinkPrefixImageDownloadTupleGenerator(unittest.TestCase):
    """Tests for `ProteinAtlasReader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_image_prefix_suffix(self):
        gen = LinkPrefixImageDownloadTupleGenerator()
        prefix, suffix = gen._get_image_prefix_suffix('https://ell-vault.stanford.edu/dav/jnhansen/www/B2AI-202308-CancerDrugs/JPGs/HPA061717_TRIM24/B2AI_1_Paclitaxel_H10_R7/z01/B2AI_1_Paclitaxel_H10_R7_z01_blue_red_green.jpg')
        self.assertEqual('https://ell-vault.stanford.edu/dav/jnhansen/www/'
                         'B2AI-202308-CancerDrugs/JPGs/HPA061717_TRIM24/B2'
                         'AI_1_Paclitaxel_H10_R7/z01/B2AI_1_Paclitaxel_H10'
                         '_R7_z01_', prefix)
        self.assertEqual('.jpg', suffix)

    def test_populate_sample_urlmap(self):
        gen = LinkPrefixImageDownloadTupleGenerator(samples_list=[{'antibody': 'HPA01',
                                                                   'filename': 'A_B_C_',
                                                                   'linkprefix': 'http://x_'},
                                                                  {'antibody': 'CAB02',
                                                                   'filename': 'D_E_F_',
                                                                   'linkprefix': 'http://y_'}])
        gen._populate_sample_urlmap()
        self.assertEqual({'1/A_B_C_': 'http://x_blue_red_green.jpg',
                          '2/D_E_F_': 'http://y_blue_red_green.jpg'},
                         gen.get_sample_urlmap())

    def test_get_next_image_url(self):
        gen = LinkPrefixImageDownloadTupleGenerator(samples_list=[{'antibody': 'HPA01',
                                                                   'filename': 'A_B_C_',
                                                                   'if_plate_id': 'A',
                                                                   'position': 'B',
                                                                   'sample': 'C',
                                                                   'linkprefix': 'http://x_'},
                                                                  {'antibody': 'CAB02',
                                                                   'filename': 'D_E_F_',
                                                                   'if_plate_id': 'D',
                                                                   'position': 'E',
                                                                   'sample': 'F',
                                                                   'linkprefix': 'http://y_'}])
        res = [x for x in gen.get_next_image_url(color_download_map={'red': '/red',
                                                                     'blue': '/blue',
                                                                     'green': '/green',
                                                                     'yellow': '/yellow'})]
        self.assertEqual(8, len(res))

        for c in ['red', 'green', 'blue', 'yellow']:
            self.assertTrue(('http://x_' + c + '.jpg', '/' + c + '/A_B_C_' +
                             c + '.jpg') in res)
            self.assertTrue(('http://y_' + c + '.jpg', '/' + c + '/D_E_F_' +
                             c + '.jpg') in res)

