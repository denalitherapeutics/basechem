import json

from django.forms import ValidationError
from django.test import TestCase
from rdkit import Chem

from basechem.common.tests.base import BasechemFormTestMixin, BasechemTestCase
from basechem.common.tests.test_constants import BAD_MOLTEXT, MMP_MOLTEXT
from basechem.main.constants import ACCEPTORS, ALOGD, CLOGP, TPSA
from basechem.main.forms import CompoundIntakeForm, MMPIntakeField


class CompoundIntakeFormTestCase(BasechemFormTestMixin, BasechemTestCase):
    def test_valid_data(self):
        with self.subTest("Good sdf, moltext, and sketcher"):
            data = {
                "project": self.test_project,
                "upload_file": self.test_sdf_file,
                "moltext": self.moltext,
                "sketcher": self.moltext,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.file_data
            )
            self.assertDictEqual(form.errors, {})
            self.assertTrue(form.is_valid())
            self.assertEqual(len(form.cleaned_data["upload_file"]), 1)
            self.assertEqual(len(form.cleaned_data["moltext"]), 1)
            self.assertEqual(len(form.cleaned_data["sketcher"]), 1)
            self.assertIsInstance(form.cleaned_data["upload_file"][0][0], Chem.Mol)
            self.assertIsInstance(form.cleaned_data["moltext"][0][0], Chem.Mol)
            self.assertIsInstance(form.cleaned_data["sketcher"][0][0], Chem.Mol)

        with self.subTest("Good sdf"):
            data = {
                "project": self.test_project,
                "upload_file": self.test_sdf_file,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.file_data
            )
            self.assertDictEqual(form.errors, {})
            self.assertTrue(form.is_valid())
            self.assertEqual(len(form.cleaned_data["upload_file"]), 1)
            self.assertIsInstance(form.cleaned_data["upload_file"][0][0], Chem.Mol)

        with self.subTest("Good mol text"):
            data = {
                "project": self.test_project,
                "moltext": self.moltext,
            }
            form = CompoundIntakeForm(current_user=self.user, data=data)
            self.assertDictEqual(form.errors, {})
            self.assertTrue(form.is_valid())
            self.assertEqual(len(form.cleaned_data["moltext"]), 1)
            self.assertIsInstance(form.cleaned_data["moltext"][0][0], Chem.Mol)

        with self.subTest("Good sketcher"):
            data = {
                "project": self.test_project,
                "sketcher": self.moltext,
            }
            form = CompoundIntakeForm(current_user=self.user, data=data)
            self.assertDictEqual(form.errors, {})
            self.assertTrue(form.is_valid())
            self.assertEqual(len(form.cleaned_data["sketcher"]), 1)
            self.assertIsInstance(form.cleaned_data["sketcher"][0][0], Chem.Mol)

    def test_invalid_data(self):
        with self.subTest("No Project good compounds"):
            data = {
                "project": "",
                "upload_file": self.test_sdf_file,
                "moltext": self.moltext,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.file_data
            )
            self.assertEqual(form.errors["project"], ["This field is required."])
            self.assertFalse(form.is_valid())
            self.assertEqual(len(form.cleaned_data["upload_file"]), 1)
            self.assertEqual(len(form.cleaned_data["moltext"]), 1)
            self.assertIsInstance(form.cleaned_data["upload_file"][0][0], Chem.Mol)
            self.assertIsInstance(form.cleaned_data["moltext"][0][0], Chem.Mol)

        with self.subTest("No mol or sdf"):
            data = {
                "project": self.test_project,
            }
            form = CompoundIntakeForm(current_user=self.user, data=data)
            self.assertEqual(form.errors["upload_file"], [""])
            self.assertEqual(form.errors["moltext"], [""])
            self.assertFalse(form.is_valid())
            self.assertRaises(KeyError, lambda: form.cleaned_data["upload_file"])
            self.assertRaises(KeyError, lambda: form.cleaned_data["moltext"])

        with self.subTest("Non sdf file"):
            data = {
                "project": self.test_project,
                "upload_file": self.test_bad_file,
                "moltext": self.moltext,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.bad_file_data
            )
            self.assertEqual(form.errors["upload_file"], ["Must be an sdf file."])
            self.assertFalse(form.is_valid())
            self.assertNotIsInstance(form.cleaned_data["upload_file"], list)
            self.assertEqual(len(form.cleaned_data["moltext"]), 1)
            self.assertIsInstance(form.cleaned_data["moltext"][0][0], Chem.Mol)

        with self.subTest("Bad compound data"):
            data = {
                "project": self.test_project,
                "upload_file": self.test_bad_sdf,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.bad_sdf_data
            )
            self.assertEqual(
                form.errors["__all__"],
                ["Your uploaded data did not return any valid compounds"],
            )
            self.assertFalse(form.is_valid())
            self.assertEqual(len(form.cleaned_data["upload_file"]), 0)
            self.assertEqual(len(form.cleaned_data["moltext"]), 0)

        with self.subTest("Bad mol data"):
            data = {
                "project": self.test_project,
                "moltext": self.bad_moltext,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.bad_file_data
            )
            self.assertEqual(
                form.errors["moltext"], ["Are you sure this is valid moltext?"]
            )
            self.assertFalse(form.is_valid())

        with self.subTest("Bad sketcher data"):
            data = {
                "project": self.test_project,
                "sketcher": self.bad_moltext,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.bad_file_data
            )
            self.assertEqual(
                form.errors["sketcher"], ["Are you sure this is a valid molecule?"]
            )
            self.assertFalse(form.is_valid())

        with self.subTest("Bad mol, sdf, and sketcher data"):
            data = {
                "project": self.test_project,
                "upload_file": self.test_bad_sdf,
                "moltext": self.bad_moltext,
                "sketcher": self.bad_moltext,
            }
            form = CompoundIntakeForm(
                current_user=self.user, data=data, files=self.bad_sdf_data
            )
            self.assertEqual(
                form.errors["moltext"], ["Are you sure this is valid moltext?"]
            )
            self.assertEqual(
                form.errors["sketcher"], ["Are you sure this is a valid molecule?"]
            )
            self.assertEqual(
                form.errors["__all__"],
                ["Your uploaded data did not return any valid compounds"],
            )
            self.assertFalse(form.is_valid())
            self.assertIsInstance(form.cleaned_data["upload_file"], list)
            self.assertEqual(len(form.cleaned_data["upload_file"]), 0)


class PropCalcFormTestCase(BasechemFormTestMixin, TestCase):
    def test_valid_form(self):
        data = {
            "project": self.test_project,
            "upload_file": self.test_sdf_file,
            "moltext": self.moltext,
            "counts": [ACCEPTORS],
            "physiochemical": [TPSA, CLOGP, ALOGD],
        }
        form = CompoundIntakeForm(
            current_user=self.user, data=data, files=self.file_data
        )
        self.assertDictEqual(form.errors, {})
        self.assertTrue(form.is_valid())

    def test_invalid_form(self):
        data = {
            "project": self.test_project,
        }
        form = CompoundIntakeForm(current_user=self.user, data=data)
        self.assertIn(
            form.errors["__all__"][0],
            "You must upload an SDF file, paste in moltext, or sketch a compound",
        )
        self.assertFalse(form.is_valid())


class MMPIntakeFieldTestCase(BasechemFormTestMixin, TestCase):
    """
    Test MMPIntakeField successfully accepts widget data from `mmp_ketcher.js` and parses it into
    smiles strings for the constant and variable regions (raising exceptions when the user's input
    will not lead to a successful mmp generate query).
    """

    def setUp(self):
        super().setUp()
        index_map = {i: i for i in range(10)}
        self.data = {
            "moltext": MMP_MOLTEXT,
            "variable": {"atoms": [6, 7], "bonds": [6]},
            "maps": {"atoms": index_map, "bonds": index_map},
        }
        self.field = MMPIntakeField()

    def test_valid(self):
        """
        Test a valid submission does not raise any exceptions
        """
        with self.subTest("CC variable region"):
            result = self.field.clean(json.dumps(self.data))
            expected_keys = ["constant_smiles", "mol", "variable_smiles"]
            self.assertEqual(sorted(list(result.keys())), expected_keys)
            self.assertEqual(Chem.MolToSmiles(result["mol"]), "CCC1CCC(CF)CC1")
            self.assertEqual(result["constant_smiles"], "*C1CCC(CC)CC1")
            self.assertEqual(result["variable_smiles"], "*CF")

        with self.subTest("CF variable region"):
            self.data["variable"]["atoms"] = [8, 9]
            result = self.field.clean(json.dumps(self.data))
            expected_keys = ["constant_smiles", "mol", "variable_smiles"]
            self.assertEqual(sorted(list(result.keys())), expected_keys)
            self.assertEqual(Chem.MolToSmiles(result["mol"]), "CCC1CCC(CF)CC1")
            self.assertEqual(result["constant_smiles"], "*C1CCC(CF)CC1")
            self.assertEqual(result["variable_smiles"], "*CC")

        with self.subTest("Variable region includes ring"):
            self.data["variable"]["atoms"] = [0, 1, 2, 3, 4, 5, 8, 9]
            result = self.field.clean(json.dumps(self.data))
            expected_keys = ["constant_smiles", "mol", "variable_smiles"]
            self.assertEqual(sorted(list(result.keys())), expected_keys)
            self.assertEqual(Chem.MolToSmiles(result["mol"]), "CCC1CCC(CF)CC1")
            self.assertEqual(result["constant_smiles"], "*CF")
            self.assertEqual(result["variable_smiles"], "*C1CCC(CC)CC1")

    def test_valid_requires_map(self):
        """
        Test a valid submission that required use of the index maps (simulate structure being deleted and redrawn)
        """
        self.data["maps"]["atoms"] = {i + 8: i for i in range(10)}
        self.data["variable"]["atoms"] = [14, 15]
        result = self.field.clean(json.dumps(self.data))
        expected_keys = ["constant_smiles", "mol", "variable_smiles"]
        self.assertEqual(sorted(list(result.keys())), expected_keys)
        self.assertEqual(Chem.MolToSmiles(result["mol"]), "CCC1CCC(CF)CC1")
        self.assertEqual(result["constant_smiles"], "*C1CCC(CC)CC1")
        self.assertEqual(result["variable_smiles"], "*CF")

    def test_invalid_no_moltext(self):
        """
        Test that an exception is raised if no molecule is drawn
        """
        self.data["moltext"] = ""
        error = "You must sketch a compound and highlight a variable region"
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))

    def test_invalid_moltext(self):
        """
        Test that an exception is raised if an invalid molecule is drawn
        """
        self.data["moltext"] = BAD_MOLTEXT
        error = "Are you sure this is a valid molecule?"
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))

    def test_invalid_multiple_molecules(self):
        """
        Test that an exception is raised if more than one molecule is drawn
        """
        mol1 = Chem.MolFromSmiles("CCC")
        mol2 = Chem.MolFromSmiles("CCCO")
        self.data["moltext"] = Chem.MolToMolBlock(Chem.CombineMols(mol1, mol2))
        error = "You may only draw one molecule at a time"
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))

    def test_invalid_missing_variable(self):
        """
        Test that an exception is raised if no variable region is highlighted
        """
        self.data["variable"] = {"atoms": [], "bonds": []}
        error = "You must highlight a variable region"
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))

    def test_invalid_variable_is_whole_mol(self):
        """
        Test that an exception is raised if the variable region is the entire molecule
        """
        self.data["variable"]["atoms"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        error = "The variable region cannot be the entire molecule. Please highlight a substructure."
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))

    def test_invalid_disconnected_variable(self):
        """
        Test that an exception is raised if the variable region is not connected
        """
        self.data["variable"]["atoms"] = [6, 7, 8, 9]
        error = "The variable region cannot be disconnected. Please highlight a connected region."
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))

    def test_invalid_split_rings(self):
        """
        Test that an exception is raised if the variable region splits a ring
        """
        self.data["variable"]["atoms"] = [5, 6, 7]
        error = "The variable region cannot include a partial ring"
        with self.assertRaisesMessage(ValidationError, error):
            result = self.field.clean(json.dumps(self.data))
