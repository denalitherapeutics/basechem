import os
import re
from unittest.mock import patch

from django.conf import settings
from django.core import mail
from django.test import tag
from rdkit import Chem

from basechem.common.mocks.mock_mayachem_utils import (
    mock_has_job_completed,
    mock_process_torsion_job_result_fail_psi4,
    mock_submit_torsion_job_to_slurm,
    mock_submit_torsion_job_to_slurm_fail_slurm_api,
    mock_submit_torsion_job_to_slurm_USE_CACHE,
)
from basechem.common.tests.base import BasechemTestCase
from basechem.main.constants import MMP
from basechem.main.models.compound_models import (
    Compound,
    compound_files_path,
    series_files_path,
)
from basechem.main.tests.factories import (
    CollectionFactory,
    CompoundFactory,
    CompoundOccurrenceFactory,
    ProjectFactory,
    SeriesFactory,
)


class CompoundQuerySetTestCase(BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.not_in = cls.collection_torsion.compounds()[0]  # propanol
        cls.benzene_smiles = "c1ccccc1"
        cls.c1 = CompoundFactory(smiles="Cc1ccccc1")  # is substruct of c2
        cls.c2 = CompoundFactory(smiles="CC(c1ccccc1)=O")  # has substruct c1
        cls.benzene_c = Compound.objects.get(smiles=cls.benzene_smiles)

    def test_has_substruct(self):
        """
        Test has_substruct returns the expected compounds; those that contain
        a benzene ring
        """
        qs_has_ss = Compound.objects.has_substruct(self.benzene_smiles)
        self.assertIn(self.benzene_c, qs_has_ss)
        self.assertIn(self.c1, qs_has_ss)
        self.assertIn(self.c2, qs_has_ss)
        self.assertEqual(len(qs_has_ss), 3)
        self.assertNotIn(self.not_in, qs_has_ss)

    def test_is_substruct(self):
        """
        Test is_substruct returns the expected compounds; those that are substructs
        of the given smiles
        """
        qs_is_ss = Compound.objects.is_substruct(self.c2.smiles)
        self.assertIn(self.benzene_c, qs_is_ss)
        self.assertIn(self.c1, qs_is_ss)
        self.assertIn(self.c2, qs_is_ss)
        self.assertEqual(len(qs_is_ss), 3)
        self.assertNotIn(self.not_in, qs_is_ss)

    def test_is_equal(self):
        """
        Test is_equal returns the expected compounds; those that are an exact match
        """
        qs_is_eq = Compound.objects.is_substruct(self.benzene_smiles)
        self.assertIn(self.benzene_c, qs_is_eq)
        self.assertEqual(len(qs_is_eq), 1)
        self.assertNotIn(self.not_in, qs_is_eq)


class CompoundTestCase(BasechemTestCase):
    def test_name(self):
        """
        Tests either VSM or DN are returned as the name of the Compound
        """
        dn_id = "DN000000001"
        self.compound_a.dn_id = dn_id
        self.compound_a.save()
        self.assertEqual(self.compound_a.name, dn_id)
        self.assertIn("VSM-", self.compound_b.name)

    def test_mol(self):
        """
        Tests an RDKit mol object is returned
        """

        mol = self.compound_a.mol()
        self.assertIsInstance(mol, Chem.Mol)
        self.assertIn("VSM-", mol.GetProp("_Name"))

    def test_mol_w_properties(self):
        """
        Tests properties are attached to the mol object
        """
        for comp in self.compounds:
            with self.subTest("Test each compound to make sure props are generated"):
                self.assertIsInstance(comp, Compound)
                mol = comp.mol_w_properties(inductive=False)
                # Check all props exists in mol
                self.assertTrue(mol.GetProp("acceptors"))
                self.assertTrue(mol.GetProp("aromaticrings"))
                self.assertTrue(mol.GetProp("charge"))
                self.assertTrue(mol.GetProp("chiralcenters"))
                self.assertTrue(mol.GetProp("donors"))
                self.assertTrue(mol.GetProp("fractioncsp3"))
                self.assertTrue(mol.GetProp("heavyatoms"))
                self.assertTrue(mol.GetProp("rings"))
                self.assertTrue(mol.GetProp("rotatablebonds"))
                self.assertTrue(mol.GetProp("clogp"))
                self.assertTrue(mol.GetProp("mw"))
                self.assertTrue(mol.GetProp("solubilityindex"))
                self.assertTrue(mol.GetProp("tpsa"))

    @tag("inductive")
    def test_mol_w_properties_inductive(self):
        """
        Tests properties are attached to the mol object including aLogD
        """
        for comp in self.compounds:
            self.assertIsInstance(comp, Compound)
            mol = comp.mol_w_properties()

            with self.subTest("Test each compound to make sure props are generated"):
                # Check all props exists in mol
                self.assertTrue(mol.GetProp("acceptors"))
                self.assertTrue(mol.GetProp("aromaticrings"))
                self.assertTrue(mol.GetProp("charge"))
                self.assertTrue(mol.GetProp("chiralcenters"))
                self.assertTrue(mol.GetProp("donors"))
                self.assertTrue(mol.GetProp("fractioncsp3"))
                self.assertTrue(mol.GetProp("heavyatoms"))
                self.assertTrue(mol.GetProp("rings"))
                self.assertTrue(mol.GetProp("rotatablebonds"))
                self.assertTrue(mol.GetProp("clogp"))
                self.assertTrue(mol.GetProp("mw"))
                self.assertTrue(mol.GetProp("solubilityindex"))
                self.assertTrue(mol.GetProp("tpsa"))
            with self.subTest("Test each compound to make sure aLogD is predicted"):
                self.assertTrue(mol.GetProp("logd_prediction"))

            with self.subTest(
                "Test LM predictions return both prediction and probabilities"
            ):
                self.assertTrue(mol.GetProp("rlm_prediction"))
                self.assertTrue(mol.GetProp("rlm_probabilities"))
                self.assertTrue(mol.GetProp("hlm_prediction"))
                self.assertTrue(mol.GetProp("hlm_probabilities"))

    def test_get_sdf_file(self):
        """
        Tests an sdf file is created with mol information
        """
        atom_0_coords = (
            "1.8843    1.9425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
        )
        atom_0_adjusted_coords = (
            "1.8846    1.9424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
        )
        with self.subTest(
            "Test that an SDF can be generated when there is no mol specified"
        ):
            filepath, filename = self.compound_a.get_sdf_file()
            self.assertTrue(os.path.exists(filepath))
            with open(filepath, "r") as f:
                r = f.read()
                self.assertNotEqual(r, "")
                self.assertIn("VSM-", r)
                self.assertIn("RDKit", r)
                self.assertIn(atom_0_coords, r)
            self.assertIn(str(self.compound_a.pk), filename)
            self.assertIn(".sdf", filename)

        with self.subTest(
            "Test that the SDF is regenerated if there is a mol specified"
        ):
            molblock = Chem.MolToMolBlock(self.compound_a.mol())
            molblock = molblock.replace(atom_0_coords, atom_0_adjusted_coords)
            mol = Chem.MolFromMolBlock(molblock)
            filepath, filename = self.compound_a.get_sdf_file(mol)
            self.assertTrue(os.path.exists(filepath))
            with open(filepath, "r") as f:
                r = f.read()
                self.assertNotEqual(r, "")
                self.assertIn("VSM-", r)
                self.assertIn("RDKit", r)
                self.assertIn(atom_0_adjusted_coords, r)
            self.assertIn(str(self.compound_a.pk), filename)
            self.assertIn(".sdf", filename)

    def test_generate_svg_file(self):
        """
        Test `generate_svg_file` generates an svg file and saves it to the model
        """
        c = CompoundFactory(smiles="CCCC")
        self.assertEqual(c.svg_file, "")
        c.generate_svg_file()
        filepath = f"{settings.MEDIA_ROOT}/{c.svg_file.name}"
        self.assertTrue(os.path.exists(filepath))
        with open(filepath, "r") as f:
            r = f.read()
            self.assertNotEqual(r, "")
            self.assertIn("<svg", r)

    def test__svg(self):
        """
        Tests mol is returned as an svg
        """
        svg = self.compound_a._svg()
        self.assertIn("<?xml version", svg)
        self.assertIn("<svg", svg)
        self.assertIn("xmlns:rdkit", svg)
        self.assertIn("width='450px'", svg)

        svg = self.compound_a._svg(size_x=150)
        self.assertIn("<?xml version", svg)
        self.assertIn("<svg", svg)
        self.assertIn("xmlns:rdkit", svg)
        self.assertIn("width='150px'", svg)

    def test_inline_svg(self):
        """
        Tests mol is returned as an svg for inline rendering
        """
        svg = self.compound_a.get_inline_svg()
        self.assertNotIn("<?xml version", svg)
        self.assertIn("<svg", svg)
        self.assertIn("xmlns:rdkit", svg)
        self.assertIn("width='450px'", svg)

        svg = self.compound_a.get_inline_svg(size_x=150)
        self.assertNotIn("<?xml version", svg)
        self.assertIn("<svg", svg)
        self.assertIn("xmlns:rdkit", svg)
        self.assertIn("width='150px'", svg)

    @tag("inductive")
    def test_predict_logd(self):
        """
        Tests aLogD is calculated with the InductiveBio model
        """
        ib_data = self.compound_a.predict_logd()
        self.assertTrue(type(ib_data) is dict)
        self.assertEqual(ib_data["ilogd"], "2.30")
        self.assertEqual(ib_data["measured"], "2.22")

    @tag("inductive")
    def test_predict_lms(self):
        """
        Tests result from InductiveBio API is returned
        """
        with self.subTest("Test species as Rat"):
            expected_dict = {
                "prediction": "54",
                "name": self.compound_a.name,
                "measured": "",
                "out_of_domain": "True",
                "latest_data_date": "2022-11-25",
                "model_version": "2.0.0",
                "interp_image": "",
                "probs_image": "",
                "probs_list": [0.049, 0.051, 0.9],
            }
            result_dict = self.compound_a.predict_lms("R", image=False)
            self.assertEqual(result_dict, expected_dict)
        with self.subTest("Test species as Human"):
            expected_dict = {
                "prediction": "20",
                "name": self.compound_a.name,
                "measured": "",
                "out_of_domain": "True",
                "latest_data_date": "2022-11-25",
                "model_version": "2.0.0",
                "interp_image": "",
                "probs_image": "",
                "probs_list": [0.2, 0.4, 0.4],
            }
            result_dict = self.compound_a.predict_lms("H", image=False)
            self.assertEqual(result_dict, expected_dict)

    def test_get_ligprep_path(self):
        """
        Tests a path to the ligprep file is returned
        """
        self.assertEqual(self.compound_a.ligprep_file.name, "")
        self.compound_a.get_sdf_file()
        filepath = self.compound_a.get_ligprep_path()
        self.assertTrue(os.path.exists(filepath))
        suppl = [m for m in Chem.SDMolSupplier(filepath, removeHs=False)]
        for mol in suppl:
            self.assertTrue(mol.GetProp("_Name").startswith("VSM-"))
            self.assertTrue(-5 < float(mol.GetProp("r_mmff_energy")) < 0)

    def test_get_confgen_path(self):
        """
        Tests a path to the confgen file is returned
        """
        self.assertEqual(self.compound_a.confgen_file.name, "")
        self.compound_a.get_ligprep_path()
        filepath = self.compound_a.get_confgen_path()
        self.assertTrue(os.path.exists(filepath))

    def test_pick_series(self):
        """
        Tests a "correct" series is selected for the compound upon original save
        """
        selected = self.compound_a.pick_series()
        self.assertEqual(self.series2, selected)
        selected = self.compound_b.pick_series()
        self.assertEqual(self.series1, selected)

    def test_get_lsalign_path(self):
        """
        Test the filepath to the flexibly aligned file is returned
        """
        self.compound_a.series = self.series2
        superimposed_file = self.compound_a.get_lsalign_path(None)
        expected_fname = f"{self.compound_a.pk}_lsaligned_to_s-{self.series2.pk}.sdf"
        self.assertIn(expected_fname, superimposed_file)
        superimposed_file = self.compound_a.get_lsalign_path(self.series2)
        self.assertIn(expected_fname, superimposed_file)

    def test_get_superimposed_path(self):
        """
        Test the filepath to the superimposed file is returned
        """
        self.compound_a.series = self.series2
        superimposed_file = self.compound_a._get_superimposed_path(None)
        expected_fname = f"{self.compound_a.pk}_superimposed_to_s-{self.series2.pk}.sdf"
        self.assertIn(expected_fname, superimposed_file)
        superimposed_file = self.compound_a._get_superimposed_path(self.series2)
        self.assertIn(expected_fname, superimposed_file)

    def test_pick_reference_file_and_id(self):
        """
        Tests a reference file is returned based on the reference
        """
        with self.subTest("Test with a supplied series as the reference"):
            reffile, ref_id = self.compound_a._pick_reference_file_and_id(self.series2)
            self.assertEqual(reffile, self.series2.bound_state_file.path)
            self.assertEqual(ref_id, f"s-{self.series2.pk}")
        with self.subTest("Test with no reference supplied"):
            self.compound_a.series = self.series2
            reffile, ref_id = self.compound_a._pick_reference_file_and_id(None)
            self.assertEqual(reffile, self.series2.bound_state_file.path)
            self.assertEqual(ref_id, f"s-{self.series2.pk}")

    @tag("local", "dtx")
    def test_update_dn_id(self):
        """
        Test `update_dn_id`
        """
        sdf_filepath, _ = self.compound_a.get_sdf_file()
        self.compound_a.get_confgen_path()

        with open(sdf_filepath, "rb") as sdf_file:
            sdf_mol = Chem.MolFromMolBlock(sdf_file.read())
            self.assertEqual(sdf_mol.GetProp("_Name"), self.compound_a.virtual_id)

        with self.subTest("Keep VSM in file if the assigned DN is the empty string"):
            self.compound_a.update_dn_id("")
            with open(sdf_filepath, "rb") as sdf_file:
                sdf_mol = Chem.MolFromMolBlock(sdf_file.read())
            self.assertEqual(sdf_mol.GetProp("_Name"), self.compound_a.virtual_id)
            self.assertEqual(self.compound_a.dn_id, "")

        with self.subTest(
            "Update the files to use the DN when there was no previous DN"
        ):
            self.compound_a.update_dn_id("DN01")
            self.assertEqual(self.compound_a.dn_id, "DN01")
            # Check SDF file is updated
            sdf_mol = Chem.MolFromMolBlock(self.compound_a.sdf_file.read())
            self.assertEqual(sdf_mol.GetProp("_Name"), "DN01")
            # Check confgen file is updated
            confgen_suppl = Chem.SDMolSupplier()
            confgen_suppl.SetData(self.compound_a.confgen_file.read())
            self.assertTrue(len(confgen_suppl) > 0)
            for confgen_mol in confgen_suppl:
                self.assertEqual(confgen_mol.GetProp("_Name"), "DN01")

        with self.subTest(
            "Update the files to use the DN when there was a different DN previously"
        ):
            self.compound_a.update_dn_id("DN02")
            self.assertEqual(self.compound_a.dn_id, "DN02")
            # Check SDF file is updated
            sdf_mol = Chem.MolFromMolBlock(self.compound_a.sdf_file.read())
            self.assertEqual(sdf_mol.GetProp("_Name"), "DN02")
            # Check confgen file is updated
            confgen_suppl = Chem.SDMolSupplier()
            confgen_suppl.SetData(self.compound_a.confgen_file.read())
            self.assertTrue(len(confgen_suppl) > 0)
            for confgen_mol in confgen_suppl:
                self.assertEqual(confgen_mol.GetProp("_Name"), "DN02")

    @tag("mmpdb")
    def test_mmp_analysis(self):
        """
        Test that `mmp_analysis` returns a list of Compound PKS
        """
        comp = self.collection_mmp.get_cos_for_analysis(MMP)[0].compound
        with self.subTest("No matches"):
            result = comp.mmp_analysis("fake smiles", "fake smiles")
            self.assertEqual(result, [])

        with self.subTest("Has matches"):
            constant = "*C1CCC(CC)CC1"
            variable = "*CF"
            result = comp.mmp_analysis(constant, variable)
            self.assertEqual(len(result), 2)
            self.assertTrue(len(result) < comp.mmps.all().count())
            mmp_comp = Compound.objects.get(pk=result[0])
            self.assertIn(mmp_comp, comp.mmps.all())

        with self.subTest("No constant and variable filters"):
            result = comp.mmp_analysis()
            self.assertEqual(len(result), 7)
            self.assertEqual(len(result), comp.mmps.all().count())
            mmp_comp = Compound.objects.get(pk=result[0])
            self.assertIn(mmp_comp, comp.mmps.all())

    @tag("mmpdb")
    def test_find_mmps(self):
        """
        Tests MMPS are added to this compound correctly
        """
        expected_compounds = Compound.objects.filter(
            dn_id__in=["DN1000004", "DN1000012"]
        )
        compound = self.collection_mmp.compounds()[1]
        expected_path = compound_files_path(compound, "generate.csv", local=True)
        self.assertFalse(compound.mmps.exists())

        compound.find_mmps()
        self.assertTrue(os.path.exists(expected_path))
        self.assertEqual(compound.mmps.all().count(), 14)
        dn_mmps = compound.mmps.exclude(dn_id="")

        self.assertEqual(dn_mmps.count(), 10)
        self.assertIn("similarity", dn_mmps.first().metadata.keys())
        self.assertIn(str(compound.pk), dn_mmps.first().metadata["similarity"].keys())
        # There are additional compounds returned but `expected_compounds` is a subset that we know are "correct" MMPs
        for cpd in expected_compounds:
            self.assertIn(cpd, dn_mmps)

    def test_update_mmp_dtx_avg_assay_data(self):
        """
        Test `update_mmp_dtx_avg_assay_data` updates average assay values for the given Compound
        and all of its MMPs
        """
        dtx_patch = patch("basechem.main.models.compound_models.get_agg_ic50_data")
        dtx_mock = dtx_patch.start()
        dtx_mock.return_value = {
            "DN0000001": {
                "data": {
                    "assay 1": {
                        "analysis 1": {"RESULT_GEOM_MEAN": 10},
                        "analysis 2": {},
                    },
                }
            },
            "DN0000002": {
                "data": {
                    "assay 1": {
                        "analysis 1": {"RESULT_GEOM_MEAN": 11},
                        "analysis 2": {"RESULT_GEOM_MEAN": 12},
                    },
                }
            },
        }
        c = CompoundFactory(dn_id="DN0000001")
        with self.subTest("Compound has no MMPs"):
            self.assertNotIn("assay_results", c.measured_data)
            c.update_mmp_dtx_avg_assay_data()
            self.assertEqual(
                c.measured_data["assay_results"], {"assay 1": {"analysis 1": 10}}
            )

        with self.subTest("Compound has MMPs"):
            mmp = CompoundFactory(dn_id="DN0000002")
            c.measured_data = {}
            c.mmps.add(mmp)
            c.save()
            self.assertNotIn("assay_results", c.measured_data)
            self.assertNotIn("assay_results", mmp.measured_data)
            c.update_mmp_dtx_avg_assay_data()
            self.assertEqual(
                c.measured_data["assay_results"], {"assay 1": {"analysis 1": 10}}
            )
            mmp.refresh_from_db()
            self.assertEqual(
                mmp.measured_data["assay_results"],
                {"assay 1": {"analysis 1": 11, "analysis 2": 12}},
            )
        dtx_patch.stop()

    @tag("local", "dtx", "external")
    def test_update_dtx_avg_assay_data(self):
        """
        Test `update_dtx_avg_assay_data` updates average assay values for the given Compound
        """
        with self.subTest("Compound is not a DN"):
            c = CompoundFactory()
            self.assertNotIn("assay_results", c.measured_data)
            c.update_dtx_avg_assay_data()
            self.assertNotIn("assay_results", c.measured_data)

        with self.subTest("Compound is a DN, no assay data"):
            c.dn_id = "DN9000001"
            c.save()
            self.assertNotIn("assay_results", c.measured_data)
            c.update_dtx_avg_assay_data()
            self.assertEqual(c.measured_data["assay_results"], {})

        with self.subTest("Compound is a DN, has assay data"):
            c.dn_id = "DN0028694"
            c.save()
            self.assertEqual(c.measured_data["assay_results"], {})
            c.update_dtx_avg_assay_data()
            assays = sorted(list(c.measured_data["assay_results"].keys()))
            self.assertEqual(len(assays), 2)
            self.assertRegex(assays[0], r".* NAD-Glo FL WX")
            self.assertRegex(assays[1], r"mouse .* NAD-Glo FL WX")
            self.assertEqual(
                c.measured_data["assay_results"][assays[0]], {"Normalized A": "0.00865"}
            )
            self.assertEqual(
                c.measured_data["assay_results"][assays[1]], {"Normalized A": "0.699"}
            )


class CompoundOccurrenceTestCase(BasechemTestCase):
    def setUp(self):
        super().setUp()
        self.co_a = CompoundOccurrenceFactory(compound=self.compound_a)

    def test_occurrences_created(self):
        """
        Tests that Compound Occurrences have been created after the collection was made
        """
        num_occurrences = self.collection.compound_occurrences.count()
        self.assertEqual(num_occurrences, 3)

    def test_occurrences_owner(self):
        """
        Tests that Compound Occurrences have the correct owner
        """
        comp_occs = self.collection.compound_occurrences.filter(
            owner=self.collection.owner
        )
        for co in comp_occs:
            self.assertEqual(co.owner, self.collection.owner)

    def test_occurrences_compound(self):
        """
        Tests that Compound Occurrences have the expected compounds for the collection
        """
        collection_cos = self.collection.compound_occurrences.all()
        compounds = [co.compound for co in collection_cos]
        self.assertIn(self.compound_a, compounds)
        self.assertIn(self.compound_b, compounds)
        self.assertIn(self.compound_c, compounds)

    def test_mol(self):
        """
        Tests Compound Occurrences can return mol objects
        """
        collection_cos = self.collection.compound_occurrences.all()
        self.assertIsInstance(collection_cos[0].mol, Chem.rdchem.Mol)

    def test_moltext(self):
        """
        Test that `moltext` returns an appropriate string of moltext
        """
        co_2d = self.collection_2d.compound_occurrences.all().first()
        co_3d = self.collection_3d.compound_occurrences.all().first()
        with self.subTest("2D comp"):
            moltext = Chem.MolToMolBlock(co_2d.compound.mol())
            self.assertEqual(self.co_a.moltext(), moltext)

        with self.subTest("3D comp"):
            self.assertEqual(co_3d.moltext(), co_3d.molblock)

        with self.subTest("3D comp with 2D flag"):
            moltext = Chem.MolToMolBlock(co_3d.compound.mol())
            self.assertEqual(co_3d.moltext(twoD=True), moltext)

    def test_get_child_cos(self):
        """
        Test that `get_child_cos` returns the correct CompoundOccurrence objects
        """
        collection = CollectionFactory()
        parent_co = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co)
        with self.subTest("No child COs"):
            qs = parent_co.get_child_cos()
            self.assertEqual(qs.count(), 0)

        with self.subTest("Has child COs"):
            child_co_1 = CompoundOccurrenceFactory(parent_co=parent_co)
            child_co_2 = CompoundOccurrenceFactory(parent_co=parent_co)

            qs = parent_co.get_child_cos()
            expected_pks = sorted([child_co_1.pk, child_co_2.pk])
            self.assertEqual(sorted([co.pk for co in qs]), expected_pks)

        with self.subTest("No child COs in collection"):
            qs = parent_co.get_child_cos(collection)
            self.assertEqual(qs.count(), 0)

        with self.subTest("Has child COs in collection"):
            collection.compound_occurrences.add(child_co_1, child_co_2)
            qs = parent_co.get_child_cos(collection)
            expected_pks = sorted([child_co_1.pk, child_co_2.pk])
            self.assertEqual(sorted([co.pk for co in qs]), expected_pks)

    def test_get_sdf_file(self):
        """
        Test `get_sdf_file` returns a path to an SDF file for both 2d and 3d CompoundOccurrences
        """
        co_2d = self.co_a
        co_3d = CompoundOccurrenceFactory(
            compound=self.co_a.compound, molblock=Chem.MolToMolBlock(self.co_a.mol)
        )
        with self.subTest("2D CO"):
            filepath, filename = co_2d.get_sdf_file()
            # Filepaths are different locally and on github, so only compare the consistent part of the filepath
            filepath = "/".join(filepath.split("/")[-4:])
            filepath_regex = f"^tmp.{{8}}/compounds/VSM-{self.compound_a.pk:011}/{self.compound_a.pk}.sdf$"
            self.assertEqual(filepath, re.match(filepath_regex, filepath)[0])
            self.assertEqual(filename, f"{self.compound_a.pk}.sdf")

        with self.subTest("3D CO"):
            filepath, filename = co_3d.get_sdf_file()
            # Filepaths are different locally and on github, so only compare the consistent part of the filepath
            filepath = "/".join(filepath.split("/")[-5:])
            filepath_regex = f"^tmp.{{8}}/compounds/VSM-{self.compound_a.pk:011}/co_{co_3d.pk}/{co_3d.pk}.sdf$"
            self.assertEqual(filepath, re.match(filepath_regex, filepath)[0])
            self.assertEqual(filename, f"{co_3d.pk}.sdf")
            # Check contents of SDF file:
            mol = Chem.SDMolSupplier(filepath)[0]
            self.assertEqual(co_3d.molblock, Chem.MolToMolBlock(mol))

    def test_get_3d_sdf_file(self):
        """
        Test `get_3d_sdf_file` returns a path to a 3d SDF file for both 2d and 3d CompoundOccurrences
        """
        co_2d = self.co_a
        co_3d = CompoundOccurrenceFactory(
            compound=self.co_a.compound, molblock=Chem.MolToMolBlock(self.co_a.mol)
        )
        comp_pk = self.compound_a.pk
        with self.subTest("2D CO"):
            # Assign a series to compound_a in order to get the default flexibly aligned file
            self.compound_a.series = self.series1
            self.compound_a.save()

            filepath, filename = co_2d.get_3d_sdf_file()
            # Filepaths are different locally and on github, so only compare the consistent part of the filepath
            filepath = "/".join(filepath.split("/")[-4:])
            filepath_regex = f"^tmp.{{8}}/compounds/VSM-{comp_pk:011}/{comp_pk}_lsaligned_to_s-{self.series1.pk}.sdf$"
            self.assertEqual(filepath, re.match(filepath_regex, filepath)[0])
            expected_filename = f"{comp_pk}_lsaligned_to_s-{self.series1.pk}.sdf"
            self.assertEqual(filename, expected_filename)

        with self.subTest("3D CO"):
            filepath, filename = co_3d.get_3d_sdf_file()
            # Filepaths are different locally and on github, so only compare the consistent part of the filepath
            filepath = "/".join(filepath.split("/")[-5:])
            filepath_regex = (
                f"^tmp.{{8}}/compounds/VSM-{comp_pk:011}/co_{co_3d.pk}/{co_3d.pk}.sdf$"
            )
            self.assertEqual(filepath, re.match(filepath_regex, filepath)[0])
            self.assertEqual(filename, f"{co_3d.pk}.sdf")
            # Check contents of SDF file:
            mol = Chem.SDMolSupplier(filepath)[0]
            self.assertEqual(co_3d.molblock, Chem.MolToMolBlock(mol))

    #######################
    ###   Align Tests   ###
    #######################

    def test_superimpose_to_ref(self):
        """
        Test a dict with confs is returned
        """
        self.compound_a.series = self.series2
        co_id = self.co_a.pk
        c_id = self.co_a.compound.pk
        with self.subTest("Test when no reference is given"):
            confs_dict = self.co_a.superimpose_to_ref(None)
            self.assertEqual(len(confs_dict), 2)
            expected_keys = [
                f"c{c_id}-co{co_id}-1-s-{self.series2.pk}",
                f"c{c_id}-co{co_id}-2-s-{self.series2.pk}",
            ]

            self.assertEqual(expected_keys, list(confs_dict.keys()))
            # Check nested dict has two objects
            self.assertEqual(
                len(confs_dict[f"c{c_id}-co{co_id}-1-s-{self.series2.pk}"]), 5
            )
            self.assertEqual(
                list(confs_dict[f"c{c_id}-co{co_id}-1-s-{self.series2.pk}"].keys()),
                [
                    "moltext",
                    "r_mmff_rel_energy",
                    "r_bc_rmsd_to_lsalign",
                    "torsion_alerts",
                    "torsion_alerts_total_energy",
                ],
            )
        with self.subTest("Test when a reference is given"):
            confs_dict = self.co_a.superimpose_to_ref(self.series2)
            self.assertEqual(len(confs_dict), 2)
            expected_keys = [
                f"c{c_id}-co{co_id}-1-s-{self.series2.pk}",
                f"c{c_id}-co{co_id}-2-s-{self.series2.pk}",
            ]

            self.assertEqual(expected_keys, list(confs_dict.keys()))
            # Check nested dict has two objects
            self.assertEqual(
                len(confs_dict[f"c{c_id}-co{co_id}-1-s-{self.series2.pk}"]), 5
            )
            self.assertEqual(
                list(confs_dict[f"c{c_id}-co{co_id}-1-s-{self.series2.pk}"].keys()),
                [
                    "moltext",
                    "r_mmff_rel_energy",
                    "r_bc_rmsd_to_lsalign",
                    "torsion_alerts",
                    "torsion_alerts_total_energy",
                ],
            )

    #######################
    ###   Dock Tests    ###
    #######################

    def test_dock_to_receptor(self):
        """
        Tests a dict of poses is returned
        """
        self.compound_a.series = self.series1
        with self.subTest("Test when no reference is given"):
            confs_dict = self.co_a.dock_to_receptor()
            # Dock results are not exactly the same every run, sometimes the best rmsd
            # and best score poses overlap to result in 3/4 results
            self.assertTrue(3 <= len(confs_dict) <= 8)
            for key in confs_dict.keys():
                self.assertEqual(
                    list(confs_dict[key].keys()),
                    [
                        "moltext",
                        "toklatScore",
                        "rdockScore",
                        "RMSDtoLSAligned",
                        "torsionAlerts",
                        "torsionAlertsTotalEnergy",
                        "toklatAnnotations",
                    ],
                )

        with self.subTest("Test when a reference is given"):
            confs_dict = self.co_a.dock_to_receptor(self.series1)
            # Dock results are not exactly the same every run, sometimes the best rmsd
            # and best score poses overlap to result in 3/4 results
            self.assertTrue(3 <= len(confs_dict) <= 8)
            for key in confs_dict.keys():
                self.assertEqual(
                    list(confs_dict[key].keys()),
                    [
                        "moltext",
                        "toklatScore",
                        "rdockScore",
                        "RMSDtoLSAligned",
                        "torsionAlerts",
                        "torsionAlertsTotalEnergy",
                        "toklatAnnotations",
                    ],
                )

    def test_get_rdock_docking_path(self):
        """
        Tests a path to the rdock output of poses is returned
        """
        self.compound_a.series = self.series1
        expected_path = compound_files_path(
            self.compound_a,
            f"{self.co_a.id}_{self.compound_a.series.dn_id}_rdock_out.sd",
            co=self.co_a,
            local=True,
        )

        localpath = self.co_a._get_rdock_output_path(self.series1)
        self.assertEqual(expected_path, localpath)

    #######################
    ###    ESP Tests    ###
    #######################

    @tag("local", "esp")
    def test_generate_esp_map(self):
        """
        Tests an ESP map is generated
        """
        self.compound_a.series = self.series1
        with self.subTest("Test when no reference is given"):
            confs_dict = self.co_a.generate_esp_map()
            self.assertEqual(len(confs_dict), 3)
            expected_keys = ["pqr", "related_series", "dx"]
            self.assertEqual(expected_keys, list(confs_dict.keys()))

        with self.subTest("Test when a reference is given"):
            confs_dict = self.co_a.generate_esp_map(self.series1)
            self.assertEqual(len(confs_dict), 3)
            expected_keys = ["pqr", "related_series", "dx"]
            self.assertEqual(expected_keys, list(confs_dict.keys()))

    #######################
    ###  Torsion Tests  ###
    #######################

    def test_convert_atoms_to_smarts(self):
        """
        Tests `convert_atoms_to_smarts` for both 2D and 3D compound occurrences
        """
        with self.subTest("2D"):
            co = self.collection_2d.compound_occurrences.first()
            dihedral_smarts = co.convert_atoms_to_smarts("1,2,3,4")
            self.assertEqual(dihedral_smarts, "[#6]-[#6]-[#6]-[#6]")

        with self.subTest("3D"):
            co = self.collection_3d.compound_occurrences.first()
            dihedral_smarts = co.convert_atoms_to_smarts("1,2,3,4")
            self.assertEqual(dihedral_smarts, "[#6]-[#6]-[#6]-[#6]")

    def test_convert_smarts_to_atoms(self):
        """
        Tests `convert_smarts_to_atoms` for both 2D and 3D compound occurrences
        """
        with self.subTest("2D"):
            co = self.collection_2d.compound_occurrences.first()
            atoms = co.convert_smarts_to_atoms("[#6]-[#6]-[#6]-[#6]")
            self.assertEqual(atoms, "0,1,2,3")

        with self.subTest("3D"):
            co = self.collection_3d.compound_occurrences.first()
            atoms = co.convert_smarts_to_atoms("[#6]-[#6]-[#6]-[#6]")
            self.assertEqual(atoms, "0,1,2,3")

    def test_pick_most_relevant_dihedral(self):
        """
        Tests `pick_most_relevant_dihedral`
        """
        pioneer_smarts = "[#6]:[#6]-[#6]-[#6]"
        pioneer_atom_indices = "8,4,3,2"
        TEST_DATA_DIR = "basechem/main/tests/testdata/torsion_dihedral_picking"

        with self.subTest("Exact match, same atom indices"):
            with open(f"{TEST_DATA_DIR}/same_atom_indices.sdf", "rb") as f:
                molblock = f.read()
            co = CompoundOccurrenceFactory(molblock=molblock)
            smarts, indices = co.pick_most_relevant_dihedral(
                pioneer_smarts, pioneer_atom_indices
            )
            self.assertEqual(smarts, pioneer_smarts)
            self.assertEqual(sorted(indices.split(",")), ["2", "3", "4", "8"])

        with self.subTest("Exact match, different atom indices"):
            with open(f"{TEST_DATA_DIR}/diff_atom_indices.sdf", "rb") as f:
                molblock = f.read()
            co = CompoundOccurrenceFactory(molblock=molblock)
            smarts, indices = co.pick_most_relevant_dihedral(
                pioneer_smarts, pioneer_atom_indices
            )
            self.assertEqual(smarts, "[#6]-[#6]-[#6]:[#6]")
            self.assertEqual(sorted(indices.split(",")), ["1", "2", "3", "4"])

        with self.subTest("One element swap"):
            with open(f"{TEST_DATA_DIR}/one_element_swap.sdf", "rb") as f:
                molblock = f.read()
            co = CompoundOccurrenceFactory(molblock=molblock)
            smarts, indices = co.pick_most_relevant_dihedral(
                pioneer_smarts, pioneer_atom_indices
            )
            self.assertEqual(smarts, "[#6]:[#6]-[#6]-[#7]")
            self.assertEqual(sorted(indices.split(",")), ["2", "3", "4", "8"])

        with self.subTest("No match"):
            with open(f"{TEST_DATA_DIR}/no_match.sdf", "rb") as f:
                molblock = f.read()
            co = CompoundOccurrenceFactory(molblock=molblock)
            smarts, indices = co.pick_most_relevant_dihedral(
                pioneer_smarts, pioneer_atom_indices
            )
            self.assertEqual(smarts, None)
            self.assertEqual(indices, None)

    def test_process_torsion_output(self):
        """
        Tests `_process_torsion_output` for an output file from torsion scan
        """
        test_output_file = "basechem/main/tests/testdata/test_torsion_results.sdf"
        co = self.collection_3d.compound_occurrences.first()
        c = co.compound
        initial_dihedral = 22.57
        results = co.process_torsion_output(test_output_file, initial_dihedral)
        # Test a conf where all dihedrals are returned
        self.assertEqual(len(results["torsions"]), 37)
        # Test nested dict is the correct size
        self.assertEqual(len(results["torsions"][f"c{c.pk}-co{co.pk}-0"]), 3)
        self.assertEqual(
            results["torsions"][f"c{c.pk}-co{co.pk}-0"]["rel_energy"], "23.888019"
        )
        self.assertEqual(results["torsions"][f"c{c.pk}-co{co.pk}-0"]["dihedral"], "0")
        # Test energy of closest dihedral to initial is delta_energy
        self.assertEqual(results["delta_energy"], "22.867507")
        self.assertEqual(results["initial_dihedral"], 22.57)

    def test_calculate_dihedral(self):
        """
        Test dihedral is correctly calculated
        """
        co = self.collection_3d.compound_occurrences.first()
        mol = Chem.MolFromMolBlock(co.molblock)
        atoms = "1,2,3,4"
        degree = co._calculate_dihedral(mol, atoms)
        self.assertEqual(degree, "-53.24661788766613")

    def test_find_closest(self):
        """
        Test closest dihedrial to the initial is returned
        """
        co = self.collection_3d.compound_occurrences.first()
        a = co._find_closest(["-90", "-60", "-30", "0", "30", "60", "90"], "22.22")
        self.assertEqual(a, "30")
        b = co._find_closest(["-90", "-60", "-30", "0", "30", "60", "90"], "-57")
        self.assertEqual(b, "-60")

    def test_clean_dihedral(self):
        """
        Test clean dihedral returns False if the torsion dihedral is bad
        """
        co = self.collection_3d.compound_occurrences.first()
        dih_good = "[#6](:[#6])-[#8]-[#6]"
        self.assertTrue(co.clean_dihedrals(dih_good))
        dih_bad = "[#6].[#6]:[#6]:[#6]"
        self.assertFalse(co.clean_dihedrals(dih_bad))

    @tag("local", "mayachem")
    @patch(
        "basechem.common.analysis_utils.submit_torsion_job_to_slurm",
        mock_submit_torsion_job_to_slurm,
    )
    @patch("basechem.common.analysis_utils.has_job_completed", mock_has_job_completed)
    @patch("basechem.common.analysis_utils.SLURM_JOB_STATUS_CHECK_PERIOD", 1)
    def test_generate_torsions(self):
        """
        Test generate torsions returns the expected dictionary
        """
        co = self.collection_torsion.compound_occurrences.all().first()
        dihedral_atoms = "0,1,2,3"
        output = co.generate_torsions(dihedral_atoms, "test job")
        self.assertIn("torsions", output)
        self.assertEqual(len(output["torsions"]), 37)
        # Delete the output file so it doesn't get used in another test
        output_path = f"{settings.MEDIA_ROOT}/compounds/{co.compound.virtual_id}/co_{co.pk}/torsions_{dihedral_atoms.replace(',','-')}.sdf"
        os.remove(output_path)

    @tag("local", "mayachem", "external", "slurm")
    @patch(
        "basechem.common.analysis_utils.submit_torsion_job_to_slurm",
        mock_submit_torsion_job_to_slurm_fail_slurm_api,
    )
    @patch("basechem.common.analysis_utils.has_job_completed", mock_has_job_completed)
    @patch("basechem.common.analysis_utils.SLURM_JOB_STATUS_CHECK_PERIOD", 1)
    def test_generate_torsions_failure_to_start_jobs(self):
        """
        Test generate torsions returns the expected dictionary when some of the dihedrals fail
        to start jobs
        """
        co = self.collection_torsion.compound_occurrences.all().first()
        dihedral_atoms = "0,1,2,3"
        output = co.generate_torsions(dihedral_atoms, "test job")
        self.assertIn("torsions", output)
        # Only dihedrals that are not multiples of 20 succeeded
        self.assertEqual(len(output["torsions"]), 18)
        # Delete the output file so it doesn't get used in another test
        output_path = f"{settings.MEDIA_ROOT}/compounds/{co.compound.virtual_id}/co_{co.pk}/torsions_{dihedral_atoms.replace(',','-')}.sdf"
        os.remove(output_path)

    @tag("local", "mayachem", "external", "slurm")
    @patch(
        "basechem.common.analysis_utils.submit_torsion_job_to_slurm",
        mock_submit_torsion_job_to_slurm_USE_CACHE,
    )
    @patch("basechem.common.analysis_utils.has_job_completed", mock_has_job_completed)
    @patch("basechem.common.analysis_utils.MAX_TORSION_ATTEMPTS", 1)
    @patch("basechem.common.analysis_utils.SLURM_JOB_STATUS_CHECK_PERIOD", 1)
    @patch(
        "basechem.common.analysis_utils.process_torsion_job_result",
        mock_process_torsion_job_result_fail_psi4,
    )
    def test_generate_torsions_job_psi4_failures(self):
        """
        Test generate torsions returns the expected dictionary when some of the dihedrals
        experience Psi4 failures.
        """
        co = self.collection_torsion.compound_occurrences.all().first()
        dihedral_atoms = "0,1,2,3"
        output = co.generate_torsions(dihedral_atoms, "test job")
        self.assertIn("torsions", output)
        # Only dihedrals -180, -90, 0, 90, and 180 succeeded
        self.assertEqual(len(output["torsions"]), 5)
        # Delete the output file so it doesn't get used in another test
        output_path = f"{settings.MEDIA_ROOT}/compounds/{co.compound.virtual_id}/co_{co.pk}/torsions_{dihedral_atoms.replace(',','-')}.sdf"
        os.remove(output_path)


class SeriesTestCase(BasechemTestCase):
    def test_series_files_path(self):
        """
        Tests the path to the series file is the expected path
        """
        path = series_files_path(self.series1, "TESTFILE.sdf")
        self.assertEqual(
            path,
            f"project_{self.series1.project.code}/{self.series1.dn_id}/{'TESTFILE.sdf'}",
        )

    def test_reassign_compounds(self):
        """
        Test that `reassign_compounds` reassigns related Compound objects to new series
        """
        project = ProjectFactory()
        s1 = SeriesFactory(
            smiles="CCCC", active=False, project=project, dn_id="DN9000001"
        )
        s2 = SeriesFactory(
            smiles="FC1=CCC(C2CCNCC2)=C1",
            project=project,
            dn_id="DN9000002",
            name="DN9000002",
        )
        s3 = SeriesFactory(
            smiles="FC1C2C(CCC2)C3=C(C=CN3)C1",
            project=project,
            dn_id="DN9000003",
            name="DN9000003",
        )
        c1 = CompoundFactory(smiles="FC1=CCC(C2CCNCC2)=C1", series=s1, project=project)
        c2 = CompoundFactory(
            smiles="FC1C2C(CCC2)C3=C(C=CN3)C1", series=s1, project=project
        )
        c3 = CompoundFactory(smiles="CC1=CCC(C2CCCCC2)=C1", series=s1)
        self.assertEqual(c1.series, s1)
        self.assertEqual(c2.series, s1)
        self.assertEqual(c3.series, s1)
        mail.outbox = []
        s1.reassign_compounds()
        # Check series assigned correctly
        c1.refresh_from_db()
        self.assertEqual(c1.series, s2)
        c2.refresh_from_db()
        self.assertEqual(c2.series, s3)
        c3.refresh_from_db()
        self.assertEqual(c3.series, s1)  # Could not reassign
        # Check email sent
        self.assertEqual(len(mail.outbox), 1)
        self.assertIn("2/3 Compounds were assigned to new Series:", mail.outbox[0].body)
        self.assertIn("1 assigned to DN9000003", mail.outbox[0].body)
        self.assertIn("1 assigned to DN9000002", mail.outbox[0].body)
        self.assertIn(
            f"The following Compounds (PKs) could not have a new Series assigned:\n{c3.pk}",
            mail.outbox[0].body,
        )
