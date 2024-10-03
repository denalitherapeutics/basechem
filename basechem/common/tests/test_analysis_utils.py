import json
import os
import shutil

from django.conf import settings
from django.test import tag
from rdkit import Chem

from basechem.common.analysis_utils import (
    combine_sdf_files,
    convert_sdf_to_mol2,
    generate_rdock_grid,
    generate_rdock_system_prm,
    generate_toklat_scores,
    generate_torsion_alerts,
    mol_to_toklat_annotations,
    mol_to_torsion_alerts,
    run_esp_predict,
)
from basechem.common.constants import TORSION_ALERTS_ENERGY_PROP, TORSION_ALERTS_PROP
from basechem.common.tests.base import BasechemTestCase


class AnalysisUtilsTestCase(BasechemTestCase):
    @tag("local", "esp")
    def test_run_esp_predict(self):
        """
        Check ESP_DNN returns the expected output
        """
        with self.subTest("Test incorrect inputs return none"):
            localpath = run_esp_predict("/fakepath")
            self.assertIsNone(localpath)

        with self.subTest("Test correct inputs return path"):
            input_dir = f"{settings.PROJECT_ROOT}/basechem/main/tests/testdata/esp_dnn"
            expected_output = f"{settings.PROJECT_ROOT}/basechem/main/tests/testdata/esp_dnn/test_onecomp.pdb.pqr"
            # run_esp_predict doesn't return a path, just check for the expected
            run_esp_predict(input_dir)
            self.assertTrue(os.path.exists(expected_output))
            os.remove(expected_output)

    @tag("local", "rdock")
    def test_generate_rdock_system_files(self):
        """
        Check generate_rdock_system_prm generates a file with the correct prefix
        """
        with self.subTest("Test generate_rdock_system_prm generates a file"):
            mol_2_file = self.series1.receptor_file_mol2.path
            ref_lig = self.series1.bound_state_file.path
            prefix = "test"
            outfile_prm_path = generate_rdock_system_prm(mol_2_file, ref_lig, prefix)
            with open(outfile_prm_path) as f:
                text = f.read()

            self.assertIn(f"{prefix}_auto_setup", text)

        with self.subTest("Test generate_rdock_grid generates two files"):
            result = generate_rdock_grid(outfile_prm_path)
            self.assertTrue(result)

            basename = os.path.splitext(outfile_prm_path)[0]
            expected_grid_file = f"{basename}_cav1.grd"
            expected_as_file = f"{basename}.as"

            self.assertTrue(os.path.exists(expected_grid_file))
            self.assertTrue(os.path.exists(expected_as_file))

            with open(expected_grid_file) as f:
                line_1 = f.readline()

            self.assertIn("GRID", line_1)

    def test_combine_sdf_files(self):
        """
        Test that `combine_sdf_files` combines an arbitrary number of SDF files into a single file
        """
        filepath1 = "basechem/main/tests/testdata/test_onecomp.sdf"
        filepath2 = "basechem/main/tests/testdata/test_onecomp_3d.sdf"
        filepath3 = "basechem/main/tests/testdata/test_onecomp_2d.sdf"

        combined_filepath = "combined.sdf"
        empty_filepath = "empty.sdf"
        with open("empty_filepath", "w"):
            pass

        with self.subTest("Combine one file"):
            combine_sdf_files(combined_filepath, [filepath1])
            mols = [mol for mol in Chem.SDMolSupplier(combined_filepath)]
            self.assertEqual(len(mols), 1)
            self.assertEqual(mols[0].GetProp("_Name"), "test_onecomp.sdf")

        with self.subTest("Combine two files"):
            combine_sdf_files(combined_filepath, [filepath1, filepath2])
            mols = [mol for mol in Chem.SDMolSupplier(combined_filepath)]
            self.assertEqual(len(mols), 2)
            self.assertEqual(mols[0].GetProp("_Name"), "test_onecomp.sdf")
            self.assertEqual(mols[1].GetProp("_Name"), "test_onecomp_3d.sdf")

        with self.subTest("Combine two files, one is empty"):
            combine_sdf_files(combined_filepath, [empty_filepath, filepath2])
            mols = [mol for mol in Chem.SDMolSupplier(combined_filepath)]
            self.assertEqual(len(mols), 1)
            self.assertEqual(mols[0].GetProp("_Name"), "test_onecomp_3d.sdf")

        with self.subTest("Combine three files"):
            combine_sdf_files(combined_filepath, [filepath1, filepath2, filepath3])
            mols = [mol for mol in Chem.SDMolSupplier(combined_filepath)]
            self.assertEqual(len(mols), 3)
            self.assertEqual(mols[0].GetProp("_Name"), "test_onecomp.sdf")
            self.assertEqual(mols[1].GetProp("_Name"), "test_onecomp_3d.sdf")
            self.assertEqual(mols[2].GetProp("_Name"), "test_onecomp_2d.sdf")

    def test_mol_to_torsion_alerts(self):
        """
        Test that `mol_to_torsion_alerts` returns a dictionary of torsion alert data
        """
        with self.subTest("No TORSION_ALERTS_PROP"):
            mol = Chem.MolFromSmiles("CC(O)C(N)c1ccc(O)cc1")
            alerts, energy = mol_to_torsion_alerts(mol)
            self.assertEqual(alerts, {})
            self.assertEqual(energy, "")

        with self.subTest("Has TORSION_ALERTS_PROP"):
            mol = Chem.MolFromSmiles("CC(O)C(N)c1ccc(O)cc1")
            alerts = "5,8 9,5,8,11 180.00 0.2606 0.1883 0.3401 CC None/anomeric [OX2:1][CX4:2]!@[CX4:3][N:4] Exact NA NA 6,8 2,6,8,11 -0.00 0.6100 0.4930 0.7541 CC None/Benzyl [a:1][c:2]!@[CX4H1:3][N,O:4] Exact NA NA"
            mol.SetProp(TORSION_ALERTS_PROP, alerts)
            mol.SetProp(TORSION_ALERTS_ENERGY_PROP, "5.48")
            expected_alerts = {"5,8": "0.2606", "6,8": "0.6100"}

            alerts, energy = mol_to_torsion_alerts(mol)

            self.assertEqual(alerts, expected_alerts)
            self.assertEqual(energy, "5.48")

    @tag("local", "mayachem")
    def test_generate_torsion_alerts(self):
        """
        Test `generate_torsion_alerts` returns a file with torsion alert data
        """
        with self.subTest("No rotatable bonds"):
            input_filepath = "basechem/main/tests/testdata/test_onecomp_3d.sdf"
            tmp_input_filepath = (
                f"{settings.MEDIA_ROOT}/torsion_alerts/no_rot_input.sdf"
            )
            tmp_output_filepath = (
                f"{settings.MEDIA_ROOT}/torsion_alerts/no_rot_input_alerts_all.sdf"
            )
            os.makedirs(os.path.dirname(tmp_input_filepath), exist_ok=True)
            shutil.copyfile(input_filepath, tmp_input_filepath)

            generate_torsion_alerts(tmp_input_filepath)
            # Check file is generated, but no torsion alerts are added
            self.assertTrue(os.path.exists(tmp_output_filepath))
            mols = [mol for mol in Chem.SDMolSupplier(tmp_output_filepath)]
            self.assertEqual(len(mols), 1)
            self.assertFalse(mols[0].HasProp(TORSION_ALERTS_PROP))
            self.assertFalse(mols[0].HasProp(TORSION_ALERTS_ENERGY_PROP))
            self.assertEqual(int(mols[0].GetProp("RotBondsCount")), 0)

        with self.subTest("Has rotatable bonds"):
            input_filepath = "basechem/main/tests/testdata/test_torsion_alert_input.sdf"
            tmp_input_filepath = f"{settings.MEDIA_ROOT}/torsion_alerts/rot_input.sdf"
            tmp_output_filepath = (
                f"{settings.MEDIA_ROOT}/torsion_alerts/rot_input_alerts_all.sdf"
            )
            os.makedirs(os.path.dirname(tmp_input_filepath), exist_ok=True)
            shutil.copyfile(input_filepath, tmp_input_filepath)

            generate_torsion_alerts(tmp_input_filepath)
            # Check file is generated and torsion alerts are added
            self.assertTrue(os.path.exists(tmp_output_filepath))
            mols = [mol for mol in Chem.SDMolSupplier(tmp_output_filepath)]
            self.assertEqual(len(mols), 7)
            for mol in mols:
                self.assertTrue(mol.HasProp(TORSION_ALERTS_PROP))
                self.assertTrue(mol.HasProp(TORSION_ALERTS_ENERGY_PROP))
                self.assertTrue(int(mols[0].GetProp("RotBondsCount")) > 0)

    def test_generate_toklat_scores(self):
        """
        Test that `generate_toklat_scores` runs the Toklat scoring model and returns a filepath with results
        """
        receptor_filepath = "/tmp/receptor.pdb"
        shutil.copyfile(
            f"{settings.TOKLAT_DIR}/example_usage/data/2RH1_CAU.pdb", receptor_filepath
        )

        with self.subTest("Test successful run"):
            pose_filepath = "/tmp/poses_success.sdf"
            shutil.copyfile(
                f"{settings.TOKLAT_DIR}/example_usage/data/2RH1_CAU.sdf", pose_filepath
            )
            original_mols = [mol for mol in Chem.SDMolSupplier(pose_filepath)]
            for mol in original_mols:
                for prop in [
                    "toklat_score",
                    "toklat_top_interactions",
                    "toklat_unsatisfied_ligand_atoms",
                    "toklat_unsatisfied_protein_atoms",
                    "toklat_error",
                ]:
                    self.assertFalse(mol.HasProp(prop))

            output_filepath = generate_toklat_scores(pose_filepath, receptor_filepath)
            self.assertEqual(output_filepath, "/tmp/poses_success_scored.sdf")
            for mol in Chem.SDMolSupplier(output_filepath):
                # Check score
                self.assertTrue(mol.HasProp("toklat_score"))
                score = float(mol.GetProp("toklat_score"))
                self.assertTrue(-2 < score < 4)
                # Check annotations
                all_annotations = []
                for prop_name in [
                    "toklat_top_interactions",
                    "toklat_unsatisfied_ligand_atoms",
                    "toklat_unsatisfied_protein_atoms",
                ]:
                    self.assertTrue(mol.HasProp(prop_name))
                    prop_value = json.loads(mol.GetProp(prop_name))
                    self.assertEqual(type(json.loads(mol.GetProp(prop_name))), list)
                    all_annotations.extend(prop_value)
                self.assertTrue(len(all_annotations) > 0)
                # Check errors
                self.assertTrue(mol.HasProp("toklat_error"))
                self.assertEqual(mol.GetProp("toklat_error"), "")

        with self.subTest("test failed run - missing Hs"):
            pose_filepath = "/tmp/poses_fail.sdf"
            shutil.copyfile(
                f"{settings.TOKLAT_DIR}/example_usage/data/2RH1_CAU_with_h_error.sdf",
                pose_filepath,
            )
            original_mols = [mol for mol in Chem.SDMolSupplier(pose_filepath)]
            for mol in original_mols:
                for prop in [
                    "toklat_score",
                    "toklat_top_interactions",
                    "toklat_unsatisfied_ligand_atoms",
                    "toklat_unsatisfied_protein_atoms",
                    "toklat_error",
                ]:
                    self.assertFalse(mol.HasProp(prop))

            output_filepath = generate_toklat_scores(pose_filepath, receptor_filepath)
            self.assertEqual(output_filepath, "/tmp/poses_fail_scored.sdf")
            mols = [mol for mol in Chem.SDMolSupplier(output_filepath)]
            self.assertEqual(len(mols), 40)
            # Check successful mols
            for mol in mols[:38]:
                # Check errors
                self.assertTrue(mol.HasProp("toklat_error"))
                self.assertEqual(mol.GetProp("toklat_error"), "")
                # Check score
                self.assertTrue(mol.HasProp("toklat_score"))
                self.assertTrue(-3 < float(mol.GetProp("toklat_score")) < 4)
                # Check annotations
                all_annotations = []
                for prop_name in [
                    "toklat_top_interactions",
                    "toklat_unsatisfied_ligand_atoms",
                    "toklat_unsatisfied_protein_atoms",
                ]:
                    self.assertTrue(mol.HasProp(prop_name))
                    prop_value = json.loads(mol.GetProp(prop_name))
                    self.assertEqual(type(json.loads(mol.GetProp(prop_name))), list)
                    all_annotations.extend(prop_value)
                self.assertTrue(len(all_annotations) > 0)

            # Check failed mols
            for mol in mols[38:]:
                # Check errors
                self.assertTrue(mol.HasProp("toklat_error"))
                self.assertEqual(
                    mol.GetProp("toklat_error"),
                    "Molecule must have explicit hydrogens.",
                )
                # Check score
                self.assertTrue(mol.HasProp("toklat_score"))
                self.assertEqual(float(mol.GetProp("toklat_score")), 999)
                # Check annotations
                for prop_name in [
                    "toklat_top_interactions",
                    "toklat_unsatisfied_ligand_atoms",
                    "toklat_unsatisfied_protein_atoms",
                ]:
                    self.assertTrue(mol.HasProp(prop_name))
                    prop_value = json.loads(mol.GetProp(prop_name))
                    self.assertEqual(type(prop_value), list)
                    self.assertEqual(prop_value, [])

    def test_mol_to_toklat_annotations(self):
        """
        Test that `mol_to_toklat_annotations` returns a dictionary of toklat annotations
        """
        pose_filepath = f"{settings.TOKLAT_DIR}/example_usage/data/2RH1_CAU.sdf"
        receptor_filepath = f"{settings.TOKLAT_DIR}/example_usage/data/2RH1_CAU.pdb"
        mol = [mol for mol in Chem.SDMolSupplier(pose_filepath, removeHs=False)][0]
        with self.subTest("Toklat properties are missing"):
            annotations = mol_to_toklat_annotations(mol, receptor_filepath)
            self.assertEqual(annotations, {"cylinders": [], "spheres": []})

        with self.subTest("Toklat properties are blank"):
            mol.SetProp("toklat_top_interactions", "[]")
            mol.SetProp("toklat_unsatisfied_ligand_atoms", "[]")
            mol.SetProp("toklat_unsatisfied_protein_atoms", "[]")
            annotations = mol_to_toklat_annotations(mol, receptor_filepath)
            self.assertEqual(annotations, {"cylinders": [], "spheres": []})

        with self.subTest("Toklat properties are present"):
            top_interactions = json.dumps(
                [
                    {
                        "li": 3,
                        "pi": 1326,
                        "pt": "O_hba_anion",
                        "lt": "N_hbd_cation",
                        "summed_interaction": -1.105778061000777,
                    },
                    {
                        "li": 6,
                        "pi": 1325,
                        "pt": "O_hba_anion",
                        "lt": "O_hba_hbd",
                        "summed_interaction": -0.7435372098784564,
                    },
                    {
                        "li": 3,
                        "pi": 6656,
                        "pt": "O_hba",
                        "lt": "N_hbd_cation",
                        "summed_interaction": -0.6751566630051443,
                    },
                    {
                        "li": 6,
                        "pi": 1326,
                        "pt": "O_hba_anion",
                        "lt": "O_hba_hbd",
                        "summed_interaction": -0.2950961484310975,
                    },
                    {
                        "li": 14,
                        "pi": 2757,
                        "pt": "O_hba",
                        "lt": "N_arom_hbd",
                        "summed_interaction": -0.2201932772477118,
                    },
                ]
            )
            unsatisfied_ligand_atoms = json.dumps(
                [
                    {
                        "li": 14,
                        "lt": "N_arom_hbd",
                        "contrib_diff_from_expected": 0.755142818301909,
                    },
                    {
                        "li": 8,
                        "lt": "O_hba",
                        "contrib_diff_from_expected": 0.33203266836914813,
                    },
                ]
            )
            unsatisfied_protein_atoms = json.dumps(
                [
                    {
                        "pi": 6657,
                        "pt": "N_hbd",
                        "contrib_diff_from_expected": 0.7567880974298684,
                    },
                    {
                        "pi": 2759,
                        "pt": "O_hba_hbd",
                        "contrib_diff_from_expected": 0.3143593525900629,
                    },
                ]
            )
            mol.SetProp("toklat_top_interactions", top_interactions)
            mol.SetProp("toklat_unsatisfied_ligand_atoms", unsatisfied_ligand_atoms)
            mol.SetProp("toklat_unsatisfied_protein_atoms", unsatisfied_protein_atoms)

            annotations = mol_to_toklat_annotations(mol, receptor_filepath)
            # Check cylinders
            self.assertEqual(len(annotations["cylinders"]), 5)
            cylinder = annotations["cylinders"][0]
            self.assertEqual(cylinder["start"], {"x": -34.079, "y": 7.375, "z": 7.505})
            self.assertEqual(cylinder["end"], {"x": -35.644, "y": 9.419, "z": 5.727})
            self.assertEqual(cylinder["color"], "blue")
            self.assertEqual(round(cylinder["radius"], 1), 0.2)
            # Check spheres
            self.assertEqual(len(annotations["spheres"]), 4)
            sphere = annotations["spheres"][0]
            self.assertEqual(sphere["center"], {"x": -25.617, "y": 10.569, "z": 5.6})

    def test_convert_sdf_to_mol2(self):
        """
        Test that `convert_sdf_to_mol2` converts an SDF file to a mol2 file
        """
        input_filepath = "basechem/main/tests/testdata/test_onecomp_3d.sdf"
        mol2_filepath = convert_sdf_to_mol2(input_filepath)
        mol_sdf = [m for m in Chem.SDMolSupplier(input_filepath, removeHs=False)]
        mol_mol2 = [Chem.MolFromMol2File(mol2_filepath, sanitize=False, removeHs=False)]
        self.assertTrue(len(mol_sdf) == len(mol_mol2))
        self.assertIsNotNone(mol_mol2)
        for i, _ in enumerate(mol_mol2):
            self.assertEqual(
                Chem.MolToSmiles(mol_sdf[i]), Chem.MolToSmiles(mol_mol2[i])
            )
