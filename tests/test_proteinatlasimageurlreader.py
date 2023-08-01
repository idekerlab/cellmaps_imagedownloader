#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ProteinAtlasImageUrlReader` package."""
import unittest
from unittest.mock import MagicMock

from cellmaps_imagedownloader.proteinatlas import ProteinAtlasImageUrlReader

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestProteinAtlasImageUrlReader(unittest.TestCase):
    """Tests for `ProteinAtlasReader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_url_from_line(self):
        reader = ProteinAtlasImageUrlReader()
        self.assertEqual('xx', reader._get_url_from_line('  <imageUrl>xx</imageUrl>\n'))
        self.assertEqual('', reader._get_url_from_line('<imageUrl></imageUrl>\n'))
        self.assertEqual('bc', reader._get_url_from_line('<a>ab<b>bc<c>\n'))

    def test_get_image_id(self):
        reader = ProteinAtlasImageUrlReader()
        self.assertEqual('39292/1495_A11_1_',
                         reader._get_image_id('http://images.proteinatlas.org/39292/1495_A11_1_blue_red_green.jpg'))

    def test_get_next_image_id_and_url(self):
        mockreader = MagicMock()
        mockreader.readline = MagicMock()

        def fake_generator():
            lines = ['http://images.proteinatlas.org/4109/1843_B2_17_blue\n',
                     '  <imageUrl>http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_selected.jpg</imageUrl>\n',
                     '<imageUrl>http://images.proteinatlas.org/4109/1832_C1_2_blue_red_green.jpg</imageUrl> \n',
                     '<imageUrl>http://images.proteinatlas.org/4109/1843_B2_17_cr5af971a263864_blue_red_green.jpg</imageUrl>']

            yield from lines

        mockreader.readline.side_effect = fake_generator
        reader = ProteinAtlasImageUrlReader(reader=mockreader)
        res = [a for a in reader.get_next_image_id_and_url()]
        self.assertEqual(2, len(res))
        self.assertEqual(('4109/1832_C1_2_',
                          'http://images.proteinatlas.org/4109/1832_C1_2_blue_'
                          'red_green.jpg'), res[0])
        self.assertEqual(('4109/1843_B2_17_',
                          'http://images.proteinatlas.org/4109/1843_B2_17_'
                          'cr5af971a263864_blue_red_green.jpg'), res[1])

