#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import os
import tempfile
import shutil

import unittest
from cellmaps_imagedownloader import cellmaps_imagedownloadercmd


class TestCellmapsDownloader(unittest.TestCase):
    """Tests for `cellmaps_imagedownloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_parse_arguments(self):
        """Tests parse arguments"""

        res = cellmaps_imagedownloadercmd._parse_arguments('hi',
                                                           ['foo'])

        self.assertEqual(res.verbose, 0)
        self.assertEqual(res.logconf, None)

        someargs = ['foo', '-vv', '--logconf',
                    'hi']
        res = cellmaps_imagedownloadercmd._parse_arguments('hi',
                                                           someargs)

        self.assertEqual(res.verbose, 2)
        self.assertEqual(res.outdir, 'foo')
        self.assertEqual(res.logconf, 'hi')

    def test_main(self):
        """Tests main function"""

        # try where loading config is successful
        try:
            temp_dir = tempfile.mkdtemp()
            res = cellmaps_imagedownloadercmd.main(['myprog.py', temp_dir])
            self.assertEqual(res, 1)
        finally:
            shutil.rmtree(temp_dir)
