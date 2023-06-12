#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `cellmaps_imagedownloader` package."""

import unittest


from cellmaps_imagedownloader.gene import GeneNodeAttributeGenerator

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestGeneNodeAttributeGenerator(unittest.TestCase):
    """Tests for `cellmaps_imagedownloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_get_gene_node_attributes(self):
        gen = GeneNodeAttributeGenerator()
        try:
            gen.get_gene_node_attributes()
            self.fail('Expected Exception')
        except NotImplementedError as ne:
            self.assertEqual('Subclasses should implement', str(ne))
