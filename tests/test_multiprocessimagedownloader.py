#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import os
import unittest
import tempfile
import shutil
import requests_mock
import json
import cellmaps_imagedownloader
from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError
from cellmaps_imagedownloader.runner import CellmapsImageDownloader
from cellmaps_imagedownloader.runner import ImageDownloader
from cellmaps_imagedownloader import runner


class TestMultiprocessImageDownloader(unittest.TestCase):
    """Tests for `cellmaps_imagedownloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_image_downloader(self):
        dloader = ImageDownloader()
        try:
            dloader.download_images()
            self.fail('Expected Exception')
        except CellMapsImageDownloaderError as ce:
            self.assertEqual('Subclasses should implement this', str(ce))
