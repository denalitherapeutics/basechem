import os

import pandas as pd
from django.conf import settings
from django.test import tag

from basechem.common.mmpdb_utils import generate_mmpdb

# from basechem.common.mocks.mock_mmpdb_utils import mock_generate_mmpdb
from basechem.common.tests.base import BasechemNoMockTestCase


@tag("mmpdb")
class MMPDBUtilsTestCase(BasechemNoMockTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.output_filepath = f"{settings.MEDIA_ROOT}/mmpdb_utils_test/generate.csv"
        os.makedirs(os.path.dirname(cls.output_filepath), exist_ok=True)
        cls.smiles = "CC1CC(C)NCN1"

    def test_generate_mmpdb(self):
        """
        Test that `generate_mmpdb` returns the expected output
        """
        self.assertFalse(os.path.exists(self.output_filepath))

        result = generate_mmpdb(self.smiles, self.output_filepath)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(self.output_filepath))

        expected_columns = [
            "start",
            "constant",
            "from_smiles",
            "to_smiles",
            "r",
            "pseudosmiles",
            "final",
            "heavies_diff",
            "#pairs",
            "pair_from_id",
            "pair_from_smiles",
            "pair_to_id",
            "pair_to_smiles",
        ]

        df = pd.read_csv(self.output_filepath, sep="\s+")

        self.assertEqual(expected_columns, df.columns.tolist())
        self.assertEqual((24, len(expected_columns)), df.shape)

        final_smiles = df["final"].tolist()
        self.assertNotIn("CC1CC(C)NCN1", final_smiles)
        self.assertIn("CCCC1CC(C)NCN1", final_smiles)
        self.assertEqual(len(final_smiles), 24)

    def tearDown(self):
        super().tearDown()
        if os.path.exists(self.output_filepath):
            os.remove(self.output_filepath)
