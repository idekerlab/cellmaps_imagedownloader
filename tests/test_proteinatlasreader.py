#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ProteinAtlasReader` package."""
import os
import unittest
import tempfile
import shutil
import gzip
from unittest.mock import patch, Mock

import requests
import requests_mock

from cellmaps_imagedownloader.exceptions import CellMapsImageDownloaderError
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

    def test_readline_download_fails_then_succeed(self):
        temp_dir = tempfile.mkdtemp()
        try:
            proteinatlas_file = os.path.join(temp_dir, 'source.xml.gz')
            with gzip.open(proteinatlas_file, 'wt') as f:
                f.write('line1\n')
                f.write('line2\n')
                f.write('line3\n')

            mock_response = Mock()
            mock_response.iter_content = lambda chunk_size: [b'line1\n', b'line2\n', b'line3\n']
            mock_response.headers = {'content-length': '15'}
            mock_exception = requests.exceptions.RequestException(response=Mock())

            p_url = 'https://hpa/proteinatlas.xml.gz'
            with requests_mock.Mocker() as m:
                with open(proteinatlas_file, 'rb') as gzfile:
                    m.get(p_url, [
                        {'exc': mock_exception},
                        {'body': gzfile}
                    ])

                    reader = ProteinAtlasReader(outdir=temp_dir, proteinatlas=p_url)
                    res = set([a for a in reader._readline(reader._proteinatlas, retry_wait=1)])
                    self.assertEqual(3, len(res))
                    self.assertTrue('line1\n' in res)
                    self.assertTrue('line2\n' in res)
                    self.assertTrue('line3\n' in res)
                self.assertEqual(m.call_count, 2)
        finally:
            shutil.rmtree(temp_dir)

    @patch('cellmaps_imagedownloader.proteinatlas.requests.get')
    def test_readline_download_fails_three_times(self, mock_get):
        temp_dir = tempfile.mkdtemp()
        try:
            with requests_mock.Mocker() as m:
                mock_response = Mock()
                mock_response.text = "Mock error response"
                mock_exception = requests.exceptions.RequestException(response=mock_response)
                mock_get.side_effect = [mock_exception] * 3
                p_url = 'https://hpa/proteinatlas.xml.gz'
                reader = ProteinAtlasReader(outdir=temp_dir, proteinatlas=p_url)
                with self.assertRaises(CellMapsImageDownloaderError):
                    for line in reader._readline(reader._proteinatlas, retry_wait=1):
                        yield line
        finally:
            shutil.rmtree(temp_dir)
        self.assertEqual(mock_get.call_count, 3)

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
