#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ImageDownloadTupleGenerator` package."""
import unittest
from unittest.mock import MagicMock

from cellmaps_imagedownloader.proteinatlas import ImageDownloadTupleGenerator

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestImageDownloadTupleGenerator(unittest.TestCase):
    """Tests for `ProteinAtlasReader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_image_prefix_suffix(self):
        mockreader = MagicMock()
        mockreader.get_next_image_id_and_url = MagicMock()

        def fake_generator():
            lines = [('','')]
            yield from lines

        mockreader.get_next_image_id_and_url.side_effect = fake_generator

        reader = ImageDownloadTupleGenerator(reader=mockreader)
        self.assertEqual(('http://images.proteinatlas.org/4109/1832_C1_2_', '.jpg'),
                         reader._get_image_prefix_suffix('http://images.pro'
                                                         'teinatlas.org/4109'
                                                         '/1832_C1_2_blue_r'
                                                         'ed_green.jpg'))

        self.assertEqual(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_', '.jpg'),
                         reader._get_image_prefix_suffix('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue_red_green.jpg'))

    def test_get_next_image_url(self):
        samples = [{'antibody': 'HPA004109',
                    'if_plate_id': '1832',
                    'position': 'C1',
                    'sample': '2'},
                   {'antibody': 'CAB0004109',
                    'if_plate_id': '1843',
                    'position': 'B2',
                    'sample': '17'},
                   {'antibody': 'HPA0001',
                    'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1'}]
        mockreader = MagicMock()
        mockreader.get_next_image_id_and_url = MagicMock()

        def fake_generator():
            lines = [('4109/1832_C1_2_',
                      'http://images.proteinatlas.org/4109/1832_C1_2_blue_red_green.jpg'),
                     ('4109/1843_B2_17_',
                     'http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue_red_green.jpg')]

            yield from lines

        mockreader.get_next_image_id_and_url.side_effect = fake_generator
        reader = ImageDownloadTupleGenerator(samples_list=samples,
                                             reader=mockreader)
        c_d_map = {'red': '/red', 'blue': '/blue',
                   'green': '/green', 'yellow': '/yellow'}
        res = [a for a in reader.get_next_image_url(color_download_map=c_d_map)]
        self.assertEqual(8, len(res))
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_red.jpg',
                         '/red/1832_C1_2_red.jpg') in res)
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_blue.jpg',
                         '/blue/1832_C1_2_blue.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_green.jpg',
                         '/green/1832_C1_2_green.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_yellow.jpg',
                         '/yellow/1832_C1_2_yellow.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_red.jpg',
                         '/red/1843_B2_17_red.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue.jpg',
                         '/blue/1843_B2_17_blue.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_green.jpg',
                         '/green/1843_B2_17_green.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_yellow.jpg',
                         '/yellow/1843_B2_17_yellow.jpg') in res)
        self.assertEqual(2, len(reader.get_sample_urlmap()))
        self.assertTrue('4109/1832_C1_2_' in list(reader.get_sample_urlmap().keys()),
                        list(reader.get_sample_urlmap().keys()))
        self.assertTrue('4109/1843_B2_17_' in list(reader.get_sample_urlmap().keys()),
                        list(reader.get_sample_urlmap().keys()))

    def test_get_next_image_url_with_valid_image_id_list(self):
        samples = [{'antibody': 'HPA004109',
                    'if_plate_id': '1832',
                    'position': 'C1',
                    'sample': '2'},
                   {'antibody': 'CAB0004109',
                    'if_plate_id': '1843',
                    'position': 'B2',
                    'sample': '17'},
                   {'antibody': 'HPA0001',
                    'if_plate_id': '1',
                    'position': 'A1',
                    'sample': '1'}]
        mockreader = MagicMock()
        mockreader.get_next_image_id_and_url = MagicMock()

        def fake_generator():
            lines = [('4109/1832_C1_2_',
                      'http://images.proteinatlas.org/4109/1832_C1_2_blue_red_green.jpg'),
                     ('4109/1843_B2_17_',
                     'http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue_red_green.jpg'),
                     ('1/1_A1_1_',
                      'http://images.proteinatlas.org/1/1_A1_1_blue_red_green.jpg')]

            yield from lines

        mockreader.get_next_image_id_and_url.side_effect = fake_generator
        reader = ImageDownloadTupleGenerator(samples_list=samples,
                                             reader=mockreader,
                                             valid_image_ids=['4109/1832_C1_2_',
                                                              '4109/1843_B2_17_'])
        c_d_map = {'red': '/red', 'blue': '/blue',
                   'green': '/green', 'yellow': '/yellow'}
        res = [a for a in reader.get_next_image_url(color_download_map=c_d_map)]
        self.assertEqual(8, len(res))
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_red.jpg',
                         '/red/1832_C1_2_red.jpg') in res)
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_blue.jpg',
                         '/blue/1832_C1_2_blue.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_green.jpg',
                         '/green/1832_C1_2_green.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1832_C1_2_yellow.jpg',
                         '/yellow/1832_C1_2_yellow.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_red.jpg',
                         '/red/1843_B2_17_red.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue.jpg',
                         '/blue/1843_B2_17_blue.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_green.jpg',
                         '/green/1843_B2_17_green.jpg') in res),
        self.assertTrue(('http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_yellow.jpg',
                         '/yellow/1843_B2_17_yellow.jpg') in res)
        self.assertEqual(2, len(reader.get_sample_urlmap()))
        self.assertTrue('4109/1832_C1_2_' in list(reader.get_sample_urlmap().keys()),
                        list(reader.get_sample_urlmap().keys()))
        self.assertTrue('4109/1843_B2_17_' in list(reader.get_sample_urlmap().keys()),
                        list(reader.get_sample_urlmap().keys()))

        

