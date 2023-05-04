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

    def test_add_geneids_to_set(self):
        gen = GeneNodeAttributeGenerator()
        # None checks
        self.assertIsNone(gen.add_geneids_to_set())
        self.assertIsNone(gen.add_geneids_to_set(gene_set=set(), geneid=None))

        # test where id lacks a delimiter
        g_set = set()
        ambiguous_dict = {}
        self.assertEqual(['FOO'],
                         gen.add_geneids_to_set(gene_set=g_set,
                                                ambiguous_gene_dict=ambiguous_dict,
                                                geneid='FOO'))
        self.assertEqual(1, len(g_set))
        self.assertTrue('FOO' in g_set)
        self.assertEqual({}, ambiguous_dict)

        # test where we have two genes
        g_set = set()
        ambiguous_dict = {}
        self.assertEqual(['ISY1', 'ISY1-RAB43'],
                         gen.add_geneids_to_set(gene_set=g_set,
                                                ambiguous_gene_dict=ambiguous_dict,
                                                geneid='ISY1,ISY1-RAB43'))
        self.assertEqual(2, len(g_set))
        self.assertTrue('ISY1' in g_set)
        self.assertTrue('ISY1-RAB43' in g_set)
        self.assertEqual({'ISY1': 'ISY1,ISY1-RAB43',
                          'ISY1-RAB43': 'ISY1,ISY1-RAB43'},
                         ambiguous_dict)






