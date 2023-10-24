#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `CM4AIImageCopyTupleGenerator` package."""
import unittest
from unittest.mock import MagicMock

from cellmaps_imagedownloader.proteinatlas import CM4AIImageCopyTupleGenerator

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestCM4AIImageCopyTupleGenerator(unittest.TestCase):
    """Tests for `CM4AIImageCopyTupleGenerator` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_image_prefix_suffix(self):
        reader = CM4AIImageCopyTupleGenerator()
        self.assertEqual(('http://images.proteinatlas.org/4109/1832_C1_2_', '.jpg'),
                         reader._get_image_prefix_suffix('http://images.pro'
                                                         'teinatlas.org/4109'
                                                         '/1832_C1_2_blue_r'
                                                         'ed_green.jpg'))

        self.assertEqual(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_', '.jpg'),
                         reader._get_image_prefix_suffix('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue_red_green.jpg'))

    def test_get_next_image_url(self):
        samples = [{'antibody': 'HPA004109',
                    'filename': '1832_C1_2_',
                    'if_plate_id': '1832',
                    'position': 'C1',
                    'sample': '2',
                    'linkprefix': 'xx/'},
                   {'antibody': 'CAB0004109',
                    'filename': '1843_B2_17_',
                    'if_plate_id': '1843',
                    'position': 'B2',
                    'sample': '17',
                    'linkprefix': 'yy/'}]

        reader = CM4AIImageCopyTupleGenerator(samples_list=samples)
        c_d_map = {'red': '/red', 'blue': '/blue',
                   'green': '/green', 'yellow': '/yellow'}
        res = [a for a in reader.get_next_image_url(color_download_map=c_d_map)]
        self.assertEqual(8, len(res))
        print(res)
        self.assertTrue(('xx/red/1832_C1_2_z01_red.jpg',
                         '/red/1832_C1_2_red.jpg') in res)
        self.assertTrue(('xx/blue/1832_C1_2_z01_blue.jpg',
                         '/blue/1832_C1_2_blue.jpg') in res),
        self.assertTrue(('xx/green/1832_C1_2_z01_green.jpg',
                         '/green/1832_C1_2_green.jpg') in res),
        self.assertTrue(('xx/yellow/1832_C1_2_z01_yellow.jpg',
                         '/yellow/1832_C1_2_yellow.jpg') in res),
        self.assertTrue(('yy/red/1843_B2_17_z01_red.jpg',
                         '/red/1843_B2_17_red.jpg') in res),
        self.assertTrue(('yy/blue/1843_B2_17_z01_blue.jpg',
                         '/blue/1843_B2_17_blue.jpg') in res),
        self.assertTrue(('yy/green/1843_B2_17_z01_green.jpg',
                         '/green/1843_B2_17_green.jpg') in res),
        self.assertTrue(('yy/yellow/1843_B2_17_z01_yellow.jpg',
                         '/yellow/1843_B2_17_yellow.jpg') in res)
        self.assertEqual(2, len(reader.get_sample_urlmap()))
        self.assertTrue('4109/1832_C1_2_' in list(reader.get_sample_urlmap().keys()),
                        list(reader.get_sample_urlmap().keys()))
        self.assertTrue('4109/1843_B2_17_' in list(reader.get_sample_urlmap().keys()),
                        list(reader.get_sample_urlmap().keys()))

        

