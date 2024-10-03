from os.path import exists

from django.test import tag
from rdkit import Chem

from basechem.common.propcalc_utils import *
from basechem.common.tests.base import BasechemTestCase
from basechem.main.constants import ALOGD, CLOGP, ROTATABLEBONDS
from basechem.main.models.collection_models import collection_files_path


class PropCalcUtilsTestCase(BasechemTestCase):
    def test_get_dtx_prop_name(self):
        with self.subTest("rot bonds"):
            prop = get_dtx_prop_name(ROTATABLEBONDS)
            self.assertEqual(prop, "RotBonds")
        with self.subTest("alogd"):
            prop = get_dtx_prop_name(ALOGD)
            self.assertEqual(prop, "AlogD")
        with self.subTest("dn_id"):
            prop = get_dtx_prop_name("dn_id")
            self.assertEqual(prop, "Compound ID")
        with self.subTest("normal"):
            prop = get_dtx_prop_name(CLOGP)
            self.assertEqual(prop, CLOGP)

    @tag("local", "inductive", "dtx", "external")
    def test_generate_dtx_propcalc_csv(self):
        """
        Test that `generate_dtx_propcalc_csv` creates a new csv file with the
        correct information. Since this runs propcalc_analysis with inductive props, it must
        be run locally
        """
        filename = "test_generate_dtx_propcalc_file"
        props_filepath = collection_files_path(self.collection, filename, local=True)
        self.assertFalse(exists(props_filepath))
        # Add DN numbers to the collection (so we can test the order they're returned):
        for i, comp in enumerate(self.collection.compounds().order_by("id"), 1):
            comp.update_dn_id(f"DN00{i}")
            # Reset the SDF file because `update_dn_id` tries to use the Dotmatics orientation
            # but there is none because DN00{i} does not actually exist in the DB.
            mol = Chem.MolFromSmiles(comp.smiles)
            mol.SetProp("_Name", comp.name)
            comp.get_sdf_file(mol=mol)
        generate_dtx_propcalc_csv(self.collection, props_filepath)

        self.assertTrue(exists(props_filepath))
        fp = open(props_filepath)
        lines = fp.readlines()

        self.assertEqual(
            lines[0],
            "Compound ID,Acceptors,AromaticRings,Charge,Donors,FractionCSP3,HeavyAtoms,MW,Rings,RotBonds,SolubilityIndex,TPSA,cLogP,AlogD\n",
        )
        self.assertEqual(lines[1], "DN001,0,0,0,0,1.0,6,84,1,0,2.3,0,2.3,2.3\n")
        self.assertEqual(lines[2], "DN002,0,1,0,0,0.0,6,78,1,0,2.7,0,1.7,2.3\n")
        self.assertEqual(lines[3], "DN003,0,0,0,0,0.2,5,66,1,0,1.5,0,1.5,2.3\n")

    @tag("local", "inductive", "dtx")
    def test_generate_dtx_lm_stability_csv(self):
        """
        Test that `generate_dtx_lm_stability_csv` creates a new csv file with the
        lm data. Since this makes calls to the InductiveBio API it must
        be run locally
        """
        date = datetime.datetime.today().strftime("%m/%d/%Y")
        filename = "test_generate_dtx_lm_stability.csv"
        lm_filepath = collection_files_path(self.collection, filename, local=True)
        self.assertFalse(exists(lm_filepath))
        # Add DN numbers to the collection (so we can test the order they're returned):
        for i, comp in enumerate(self.collection.compounds().order_by("id"), 1):
            comp.update_dn_id(f"DN00{i}")
        generate_dtx_lm_stability_csv(self.collection, lm_filepath)

        self.assertTrue(exists(lm_filepath))
        fp = open(lm_filepath)
        lines = fp.readlines()

        self.assertEqual(
            lines[0],
            "name,assay,skip,skip,skip,out_of_domain_flag,skip,pStable,skip,skip,prediction_date,model_version,prediction_hlm,skip,prediction_rlm,skip\n",
        )
        self.assertEqual(
            lines[1],
            f"DN001,Inductive Bio GCNN,,,,out-of-domain,,0.200,,,{date},2.0.0: 2022-11-25,20,,54,\n",
        )
        self.assertEqual(
            lines[2],
            f"DN002,Inductive Bio GCNN,,,,out-of-domain,,0.200,,,{date},2.0.0: 2022-11-25,20,,54,\n",
        )
        self.assertEqual(
            lines[3],
            f"DN003,Inductive Bio GCNN,,,,out-of-domain,,0.200,,,{date},2.0.0: 2022-11-25,20,,54,\n",
        )
