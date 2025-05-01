#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `CM4AITableConverter` """

import os
import csv
import unittest
import tempfile
import shutil
import json
from unittest.mock import MagicMock

import pandas

from cellmaps_imagedownloader.gene import CM4AITableConverter

SKIP_REASON = 'CELLMAPS_IMAGEDOWNLOADER_INTEGRATION_TEST ' \
              'environment variable not set, cannot run integration ' \
              'tests'


class TestCM4AITableConverter(unittest.TestCase):
    """Tests for `CM4AITableConverter`."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def _write_fake_cm4ai_table_to_file_zeropointsixversion(self, dest_file=None):

        df = pandas.DataFrame(data={'Plate': [3, 3],
                                    'Well': ['C2', 'D1'],
                                    'Treatment': ['Paclitaxel', 'Paclitaxel'],
                                    'HPA_Antibody_ID': ['CAB005042', 'HPA000895'],
                                    'ENSEMBL ID': ['ENSG00000132507',
                                                   'ENSG00000102024;ENSG00000120756;ENSG00000136167'],
                                    'Uniprot': ['P63241', 'P13797;P13796;Q14651'],
                                    'Gene_Name_HPA': ['EIF5A', 'Q14651,PLS3_LCP1_PLS1'],
                                    'Chromosome': [17, 'X;13;3'],
                                    'Baselink': ['https://ell-vault.stanford.edu/dav/jnhansen/www/B2AI/MDA-MB-468/CAB005042_EIF5A/B2AI_3_Paclitaxel_C2_R1/z02/B2AI_3_Paclitaxel_C2_R1_z02',
                                                 'https://ell-vault.stanford.edu/dav/jnhansen/www/B2AI/MDA-MB-468/HPA000895_PLS3_LCP1_PLS1/B2AI_3_Paclitaxel_D1_R1/z00/B2AI_3_Paclitaxel_D1_R1_z00']})

        df.to_csv(dest_file)


    def _write_fake_cm4ai_table_to_file(self, dest_file=None):
        """

        :param dest_file:
        :return:
        """
        with open(dest_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=['Antibody ID',
                                                                   'ENSEMBL ID',
                                                                   'Treatment',
                                                                   'Well',
                                                                   'Region'])
            writer.writeheader()
            writer.writerow({'Antibody ID': 'CAB079904',
                             'ENSEMBL ID': 'ENSG00000187555',
                             'Treatment': 'Paclitaxel',
                             'Well': 'C1',
                             'Region': 'R1'})
            writer.writerow({'Antibody ID': 'CAB079904',
                             'ENSEMBL ID': 'ENSG00000187555',
                             'Treatment': 'Paclitaxel',
                             'Well': 'C1',
                             'Region': 'R2'})
            writer.writerow({'Antibody ID': 'HPA059206',
                             'ENSEMBL ID': 'ENSG00000101266;ENSG00000254598',
                             'Treatment': 'Paclitaxel',
                             'Well': 'A6',
                             'Region': 'R1'})
            writer.writerow({'Antibody ID': 'HPA059206',
                             'ENSEMBL ID': 'ENSG00000101266;ENSG00000254598',
                             'Treatment': 'Paclitaxel',
                             'Well': 'A6',
                             'Region': 'R12'})

    def test_get_samples_from_cm4ai_table_as_dataframe(self):
        temp_dir = tempfile.mkdtemp()
        try:
            cm4ai_table = os.path.join(temp_dir, 'cm4ai.tsv')
            self._write_fake_cm4ai_table_to_file(cm4ai_table)

            converter = CM4AITableConverter(cm4ai=cm4ai_table)
            df = converter._get_samples_from_cm4ai_table_as_dataframe(cm4ai_table)
            self.assertEqual(['filename', 'if_plate_id', 'position',
                              'sample', 'locations',
                              'antibody', 'ensembl_ids', 'gene_names',
                              'linkprefix'],df.columns.tolist())
            self.assertEqual(4, len(df))
        finally:
            shutil.rmtree(temp_dir)

    def test_get_unique_dataframe_from_samples_dataframe(self):
        temp_dir = tempfile.mkdtemp()
        try:
            cm4ai_table = os.path.join(temp_dir, 'cm4ai.tsv')
            self._write_fake_cm4ai_table_to_file(cm4ai_table)

            converter = CM4AITableConverter(cm4ai=cm4ai_table)
            df = converter._get_samples_from_cm4ai_table_as_dataframe(cm4ai_table)

            unique_df = converter._get_unique_dataframe_from_samples_dataframe(samples_df=df)
            self.assertEqual(['antibody', 'ensembl_ids', 'gene_names',
                              'atlas_name', 'locations', 'n_location'],
                             unique_df.columns.tolist())
            self.assertEqual(2, len(unique_df))
        finally:
            shutil.rmtree(temp_dir)


    def test_get_samples_from_cm4ai_tablezeropointsix_as_dataframe(self):
        temp_dir = tempfile.mkdtemp()
        try:
            cm4ai_table = os.path.join(temp_dir, 'cm4ai.csv')
            self._write_fake_cm4ai_table_to_file_zeropointsixversion(cm4ai_table)

            converter = CM4AITableConverter(cm4ai=cm4ai_table)
            df = converter._get_samples_from_cm4ai_table_as_dataframe_from_updatedversion(cm4ai_table)
            self.assertEqual(['filename', 'if_plate_id', 'position',
                              'sample', 'locations',
                              'antibody', 'ensembl_ids', 'gene_names',
                              'linkprefix'],df.columns.tolist())
            self.assertEqual(2, len(df))
        finally:
            shutil.rmtree(temp_dir)

    def test_get_samples_and_unique_lists_unset_in_constructor(self):
        converter = CM4AITableConverter()
        self.assertEqual((None, None), converter.get_samples_and_unique_lists())

    def test_get_samples_and_unique_lists(self):
        temp_dir = tempfile.mkdtemp()
        try:
            cm4ai_table = os.path.join(temp_dir, 'cm4ai.tsv')
            self._write_fake_cm4ai_table_to_file(cm4ai_table)

            converter = CM4AITableConverter(cm4ai=cm4ai_table)
            samples, unique = converter.get_samples_and_unique_lists()
            self.assertEqual(4, len(samples))
            self.assertEqual(2, len(unique))
            # Todo test all fields are set properly
            self.assertEqual('MDA-MB-468', unique[0]['atlas_name'])
        finally:
            shutil.rmtree(temp_dir)






