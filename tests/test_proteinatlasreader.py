#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ProteinAtlasReader` package."""
import os
import unittest
import tempfile
import shutil


from cellmaps_imagedownloader.proteinatlas import ProteinAtlasReader

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestProteinAtlasReader(unittest.TestCase):
    """Tests for `ProteinAtlasReader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_readline_with_standard_txt_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            reader = ProteinAtlasReader(temp_dir)
            proteinatlas_file = os.path.join(temp_dir, 'proteinatlas.xml')
            with open(proteinatlas_file, 'w') as f:
                f.write('line1\n')
                f.write('line2\n')
                f.write('line3\n')
            res = [a for a in reader.readline(proteinatlas_file)]
            self.assertEqual([], res)


        finally:
            shutil.rmtree(temp_dir)
