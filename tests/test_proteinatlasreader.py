#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ProteinAtlasReader` package."""
import os
import unittest
import tempfile
import shutil
import gzip
import requests
import requests_mock


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

            proteinatlas_file = os.path.join(temp_dir, 'proteinatlas.xml')
            with open(proteinatlas_file, 'w') as f:
                f.write('line1\n')
                f.write('line2\n')
                f.write('line3\n')
            reader = ProteinAtlasReader(outdir=temp_dir,
                                        proteinatlas=proteinatlas_file)

            res = set([a for a in reader.readline()])
            self.assertEqual(3, len(res))
            self.assertTrue('line1\n' in res)
            self.assertTrue('line2\n' in res)
            self.assertTrue('line3\n' in res)
        finally:
            shutil.rmtree(temp_dir)

    def test_readline_with_standard_gzip_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            proteinatlas_file = os.path.join(temp_dir, 'proteinatlas.xml.gz')
            with gzip.open(proteinatlas_file, 'wt') as f:
                f.write('line1\n')
                f.write('line2\n')
                f.write('line3\n')

            reader = ProteinAtlasReader(outdir=temp_dir,
                                        proteinatlas=proteinatlas_file)

            res = set([a for a in reader.readline()])
            self.assertEqual(3, len(res))
            self.assertTrue('line1\n' in res)
            self.assertTrue('line2\n' in res)
            self.assertTrue('line3\n' in res)
        finally:
            shutil.rmtree(temp_dir)

    def test_readline_with_gzip_url(self):
        temp_dir = tempfile.mkdtemp()
        try:
            proteinatlas_file = os.path.join(temp_dir, 'source.xml.gz')
            with gzip.open(proteinatlas_file, 'wt') as f:
                f.write('line1\n')
                f.write('line2\n')
                f.write('line3\n')

            with requests_mock.Mocker() as m:
                with open(proteinatlas_file, 'rb') as gzfile:
                    p_url = 'https://hpa/proteinatlas.xml.gz'
                    m.get(p_url, body=gzfile)
                    reader = ProteinAtlasReader(outdir=temp_dir,
                                                proteinatlas=p_url)
                    res = set([a for a in reader.readline()])
                    self.assertEqual(3, len(res))
                    self.assertTrue('line1\n' in res)
                    self.assertTrue('line2\n' in res)
                    self.assertTrue('line3\n' in res)
        finally:
            shutil.rmtree(temp_dir)

    def test_readline_with_gzip_url_none_for_proteinatlas(self):
        temp_dir = tempfile.mkdtemp()
        try:
            proteinatlas_file = os.path.join(temp_dir, 'source.xml.gz')
            with gzip.open(proteinatlas_file, 'wt') as f:
                f.write('line1\n')
                f.write('line2\n')
                f.write('line3\n')

            with requests_mock.Mocker() as m:
                with open(proteinatlas_file, 'rb') as gzfile:
                    p_url = ProteinAtlasReader.DEFAULT_PROTEINATLAS_URL
                    m.get(p_url, body=gzfile)
                    reader = ProteinAtlasReader(outdir=temp_dir)
                    res = set([a for a in reader.readline()])
                    self.assertEqual(3, len(res))
                    self.assertTrue('line1\n' in res)
                    self.assertTrue('line2\n' in res)
                    self.assertTrue('line3\n' in res)
        finally:
            shutil.rmtree(temp_dir)

    """
    # @unittest.skipUnless(os.getenv('CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST') is not None, SKIP_REASON)
    def test_real_download_of_url(self):
        temp_dir = tempfile.mkdtemp()
        try:
            reader = ProteinAtlasReader(temp_dir)
            for line in reader.readline('https://www.proteinatlas.org/download/proteinatlas.xml.gz'):
                print(line)
        finally:
            shutil.rmtree(temp_dir)
    """
