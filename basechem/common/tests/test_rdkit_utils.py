from rdkit import Chem

from basechem.common.rdkit_utils import RDKitWrappers
from basechem.common.tests.base import BasechemTestCase
from basechem.common.tests.test_constants import BAD_MOLTEXT, MOLTEXT


class TestRDKitUtils(BasechemTestCase):
    def test_is_3d(self):
        """
        Test `RDKitWrappers.is_3d` returns True ifa mol has 3D coordinates and False otherwise
        """
        with self.subTest("3D is correctly identified"):
            mol = self.collection_3d.compound_occurrences.all()[0].mol
            self.assertTrue(RDKitWrappers.is_3d(mol))

        with self.subTest("2D is correctly identified"):
            mol = Chem.MolFromMolBlock(MOLTEXT)
            self.assertFalse(RDKitWrappers.is_3d(mol))

        with self.subTest("2D from smiles is correctly identified"):
            mol = Chem.MolFromSmiles("CCCC")
            self.assertFalse(RDKitWrappers.is_3d(mol))

    def test_clean_mol_object(self):
        """
        Tests `RDKitWrappers.clean_mol_object` cleans mol objects
        """
        romol = Chem.MolFromMolBlock(MOLTEXT)
        mols = RDKitWrappers.clean_mol_object(romol)
        self.assertEqual(len(mols), 1)
        self.assertTrue(type(mols[0][0]), Chem.Mol)
        self.assertTrue(mols[0][1])

        romol = Chem.MolFromMolBlock(BAD_MOLTEXT)
        mols = RDKitWrappers.clean_mol_object(romol)
        self.assertEqual(len(mols), 0)

    def test_ligprep(self):
        """
        Tests `RDKitWrappers.ligprep` generates a set of 3D structures
        """
        for compound in self.collection.compounds():
            with self.subTest(f"Check Gypsum ran for {compound.name}"):
                self.assertFalse(RDKitWrappers.is_3d(compound.mol()))
                ligprep_path = RDKitWrappers.ligprep(
                    compound.mol(), f"/tmp/test_ligprep_{compound.name}.sdf"
                )
                ligprep_mols = [
                    m for m in Chem.SDMolSupplier(ligprep_path, removeHs=False)
                ]
                self.assertTrue(all(RDKitWrappers.is_3d(mol) for mol in ligprep_mols))
                self.assertTrue(
                    all(
                        "Genealogy" in mol.GetPropsAsDict().keys()
                        for mol in ligprep_mols
                    )
                )

    def test_pick_series(self):
        """
        Tests `RDKitWrappers.pick_series`
        """
        series_tuples = [
            (Chem.MolFromSmiles("C1=CC=CC=C1"), "DN0000001"),
            (Chem.MolFromSmiles("C1CCCCC1"), "DN0000002"),
        ]

        with self.subTest("C1CCCCC1"):
            mol = Chem.MolFromSmiles("C1CCCCC1")
            result = RDKitWrappers.pick_series(mol, series_tuples)
            self.assertEqual(result, series_tuples[1][1])

        with self.subTest("C1=CC=CC=C1"):
            mol = Chem.MolFromSmiles("C1=CC=CC=C1")
            result = RDKitWrappers.pick_series(mol, series_tuples)
            self.assertEqual(result, series_tuples[0][1])

        with self.subTest("c1ccccc1"):
            mol = Chem.MolFromSmiles("c1ccccc1")
            result = RDKitWrappers.pick_series(mol, series_tuples)
            self.assertEqual(result, series_tuples[0][1])

    def test_add_hydrogens_to_sdf(self):
        """
        Test that `add_hydrogens_to_sdf` adds explicit hydrogens to all mols in a SDF file
        """
        original_mol = Chem.MolFromSmiles("CCCC")
        filepath = "/tmp/hydrogens.sdf"
        w = Chem.SDWriter(filepath)
        w.write(original_mol)
        w.close()
        self.assertEqual(original_mol.GetNumAtoms(), 4)
        filepath = RDKitWrappers.add_hydrogens_to_sdf(filepath)
        new_mol = [mol for mol in Chem.SDMolSupplier(filepath, removeHs=False)][0]
        self.assertEqual(new_mol.GetNumAtoms(), 14)
