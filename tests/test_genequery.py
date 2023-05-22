#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import os
import unittest
import tempfile
import shutil
import json
from unittest.mock import MagicMock
from cellmaps_imagedownloader.gene import GeneQuery

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestGeneQuery(unittest.TestCase):
    """Tests for `cellmaps_downloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_querymany(self):
        mockquery = MagicMock()

        mockquery.querymany = MagicMock(return_value={})
        query = GeneQuery(mygeneinfo=mockquery)
        res = query.querymany(queries=['hi'], scopes='thescope',
                              fields=['field1'],
                              species='human')
        self.assertEqual({}, res)
        mockquery.querymany.assert_called_once_with(['hi'], scopes='thescope',
                                                    fields=['field1'],
                                                    species='human')

    def test_get_symbols_for_genes(self):
        mockquery = MagicMock()

        mockquery.querymany = MagicMock(return_value={})
        query = GeneQuery(mygeneinfo=mockquery)
        res = query.get_symbols_for_genes(genelist=['g1'])
        self.assertEqual({}, res)
        mockquery.querymany.assert_called_once_with(['g1'], scopes='_id',
                                                    fields=['ensembl.gene', 'symbol'],
                                                    species='human')

    @unittest.skipUnless(os.getenv('CELLMAPS_DOWNLOADER_INTEGRATION_TEST') is not None, SKIP_REASON)
    def test_simple_query(self):
        query = GeneQuery()
        res = query.querymany([2, 16], scopes='_id', species='human',
                              fields=['ensembl.gene','symbol'])
        self.assertEqual(2, len(res))
        for entry in res:
            if entry['query'] == '2':
                self.assertEqual({'query': '2',
                                  '_id': '2', '_score': 1.55,
                                  'ensembl': {'gene': 'ENSG00000175899'},
                                  'symbol': 'A2M'}, entry)
            elif entry['query'] == '16':
                self.assertEqual({'query': '16',
                                  '_id': '16', '_score': 1.55,
                                  'ensembl': {'gene': 'ENSG00000090861'},
                                  'symbol': 'AARS1'}, entry)
            else:
                self.fail('Unexpected entry: ' + str(entry))
