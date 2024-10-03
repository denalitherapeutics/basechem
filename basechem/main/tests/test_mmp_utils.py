import os

from django.test import tag
from rdkit import Chem

from basechem.common.tests.base import BasechemNoMockTestCase
from basechem.main.mmp_utils import generate_mmps, meets_mmp_prop_req
from basechem.main.models.compound_models import Compound, compound_files_path


class MMPUtilsTestCase(BasechemNoMockTestCase):
    def test_meets_mmp_prop_req(self):
        """
        Tests `meets_mmp_prop_req` filters correctly
        """
        not_meets_smiles = "CCCCC1CCC2(C1)CC13CCCC1C1CC(C)CNC14CCCC34C2"
        not_meets_mol = Chem.MolFromSmiles(not_meets_smiles)
        meets_smiles = self.collection_mmp.compounds()[0].smiles
        meets_mol = Chem.MolFromSmiles(meets_smiles)

        self.assertFalse(meets_mmp_prop_req(not_meets_mol))
        self.assertTrue(meets_mmp_prop_req(meets_mol))

    @tag("mmpdb")
    def test_generate_mmps(self):
        """
        Tests that 'generate_mmps' returns a dataframe from mmpdb results
        """
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
            "final_mol",
            "meets_mmp_prop_req",
        ]

        compound = Compound.objects.get(dn_id="DN1000003")
        generate_path = compound_files_path(compound, "generate.csv", local=True)
        mmp_df = generate_mmps(compound.smiles, generate_path)

        self.assertTrue(os.path.exists(generate_path))
        self.assertEqual(mmp_df.shape[0], 10)
        self.assertListEqual(mmp_df.columns.to_list(), expected_columns)
