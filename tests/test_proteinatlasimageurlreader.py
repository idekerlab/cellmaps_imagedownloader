#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ProteinAtlasImageUrlReader` package."""
import os
import unittest
import tempfile
import shutil
import gzip
import requests
import requests_mock


from cellmaps_imagedownloader.proteinatlas import ProteinAtlasImageUrlReader
from cellmaps_imagedownloader.proteinatlas import ProteinAtlasReader

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



