import datetime
import os
import time

from django.conf import settings
from django.core import mail
from django.test import override_settings, tag
from django.urls import reverse
from django_q.tasks import async_task, fetch_group
from rdkit import Chem

from basechem.common.tests.base import BasechemTestCase
from basechem.main.constants import (
    ALIGN,
    ALOGD,
    AUTO,
    CLOGP,
    DOCK,
    ESP,
    HLM,
    MMP,
    MW,
    PROPCALC,
    RLM,
    TORSION,
)
from basechem.main.models.compound_models import Compound
from basechem.main.tasks import get_group_results
from basechem.main.tests.factories import (
    CollectionFactory,
    CompoundFactory,
    CompoundOccurrenceFactory,
    ProjectFactory,
)
from basechem.users.tests.factories import BasechemUserFactory


class CollectionTestCase(BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = BasechemUserFactory()

    def test_get_co_order(self):
        """
        Test `get_co_order` returns self.metadata["co_order"] and adds missing COs to the list
        """
        collection = CollectionFactory()
        with self.subTest("no COs"):
            co_order = collection.get_co_order()
            self.assertEqual(co_order, [])
        parent_co_a = CompoundOccurrenceFactory()
        parent_co_b = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co_a, parent_co_b)

        with self.subTest("No initial co_order"):
            collection.metadata = {}
            collection.save()
            co_order = collection.get_co_order()
            expected_co_order = [parent_co_a.pk, parent_co_b.pk]
            self.assertEqual(co_order, expected_co_order)

        a_child_co1 = CompoundOccurrenceFactory(parent_co=parent_co_a)
        a_child_co2 = CompoundOccurrenceFactory(parent_co=a_child_co1)
        collection.compound_occurrences.add(a_child_co1, a_child_co2)
        with self.subTest("Auto-add child COs"):
            co_order = collection.get_co_order()
            expected_co_order = [
                parent_co_a.pk,
                a_child_co1.pk,
                a_child_co2.pk,
                parent_co_b.pk,
            ]
            self.assertEqual(co_order, expected_co_order)

        with self.subTest("Subset of cos still returns full co_order"):
            co_order = collection.get_co_order(
                cos=collection.compound_occurrences.exclude(parent_co=None)
            )
            expected_co_order = [
                parent_co_a.pk,
                a_child_co1.pk,
                a_child_co2.pk,
                parent_co_b.pk,
            ]
            self.assertEqual(co_order, expected_co_order)

    def test_get_cos_in_order(self):
        """
        Test `get_cos_in_order` returns the expected queryset of COs
        """
        collection = CollectionFactory()
        with self.subTest("no COs"):
            ordered_cos = collection.get_cos_in_order()
            self.assertEqual([co.pk for co in ordered_cos], [])
        cos = []
        for i in range(5):
            co = CompoundOccurrenceFactory()
            cos.append(co)
        collection.compound_occurrences.add(*[co.pk for co in cos])
        collection.save()

        with self.subTest("Original order"):
            cos_in_order = collection.get_cos_in_order()
            expected_co_pks = [co.pk for co in cos]
            self.assertEqual(expected_co_pks, [co.pk for co in cos_in_order])
            self.assertEqual(collection.metadata["co_order"], expected_co_pks)

        with self.subTest("Reverse order"):
            collection.metadata["co_order"] = [co.pk for co in cos][::-1]
            collection.save()
            cos_in_order = collection.get_cos_in_order()
            expected_co_pks = [co.pk for co in cos][::-1]
            self.assertEqual(expected_co_pks, [co.pk for co in cos_in_order])
            self.assertEqual(collection.metadata["co_order"], expected_co_pks)

        with self.subTest("Curated order"):
            expected_cos_in_order = [cos[1], cos[2], cos[0], cos[4], cos[3]]
            expected_co_pks = [co.pk for co in expected_cos_in_order]
            collection.metadata["co_order"] = expected_co_pks
            collection.save()
            cos_in_order = collection.get_cos_in_order()
            self.assertEqual(expected_co_pks, [co.pk for co in cos_in_order])
            self.assertEqual(collection.metadata["co_order"], expected_co_pks)

    def test_update_co_order(self):
        """
        Test `update_co_order` correctly reorders the co_order to accommodate the given order
        """
        collection = CollectionFactory()
        with self.subTest("Add COs for the first time"):
            collection.update_co_order([1, 2, 3])
            collection.refresh_from_db()
            self.assertEqual(collection.metadata["co_order"], [1, 2, 3])

        with self.subTest("Add new CO"):
            collection.update_co_order([4])
            collection.refresh_from_db()
            self.assertEqual(collection.metadata["co_order"], [4, 1, 2, 3])

        with self.subTest("Reorder existing COs"):
            collection.update_co_order([1, 4])
            collection.refresh_from_db()
            self.assertEqual(collection.metadata["co_order"], [1, 4, 2, 3])

    def test_get_cos_for_analysis_propcalc(self):
        """
        Test `get_cos_for_analysis` when the analysis is PROPCALC
        """
        collection = CollectionFactory()
        parent_co_1 = CompoundOccurrenceFactory()
        parent_co_2 = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co_1, parent_co_2)

        with self.subTest("2D COs, no basechem generated COs"):
            qs = collection.get_cos_for_analysis(PROPCALC)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("2D COs w/ basechem generated COs"):
            child_co_1 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            child_co_2 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            collection.compound_occurrences.add(child_co_1, child_co_2)

            qs = collection.get_cos_for_analysis(PROPCALC)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("3D CO w/ basechem generated COs"):
            parent_co_1.molblock = "Not empty"
            parent_co_1.save()

            qs = collection.get_cos_for_analysis(PROPCALC)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

    def test_get_cos_for_analysis_align(self):
        """
        Test `get_cos_for_analysis` when the analysis is ALIGN
        """
        collection = CollectionFactory()
        parent_co_1 = CompoundOccurrenceFactory()
        parent_co_2 = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co_1, parent_co_2)

        with self.subTest("2D COs, no basechem generated COs"):
            qs = collection.get_cos_for_analysis(DOCK)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("2D COs w/ basechem generated COs"):
            child_co_1 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            child_co_2 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=DOCK
            )
            child_co_3 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=DOCK
            )
            collection.compound_occurrences.add(child_co_1, child_co_2, child_co_3)

            qs = collection.get_cos_for_analysis(ALIGN)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("3D CO w/ basechem generated COs"):
            parent_co_1.molblock = "Not empty"
            parent_co_1.save()

            qs = collection.get_cos_for_analysis(ALIGN)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

    def test_get_cos_for_analysis_dock(self):
        """
        Test `get_cos_for_analysis` when the analysis is DOCK
        """
        collection = CollectionFactory()
        parent_co_1 = CompoundOccurrenceFactory()
        parent_co_2 = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co_1, parent_co_2)

        with self.subTest("2D COs, no basechem generated COs"):
            qs = collection.get_cos_for_analysis(DOCK)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("2D COs w/ basechem generated COs"):
            child_co_1 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=DOCK
            )
            child_co_2 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            child_co_3 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            collection.compound_occurrences.add(child_co_1, child_co_2, child_co_3)

            qs = collection.get_cos_for_analysis(DOCK)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("3D CO w/ basechem generated COs"):
            parent_co_1.molblock = "Not empty"
            parent_co_1.save()

            qs = collection.get_cos_for_analysis(DOCK)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

    def test_get_cos_for_analysis_esp(self):
        """
        Test `get_cos_for_analysis` when the analysis is ESP
        """
        collection = CollectionFactory()
        parent_co_1 = CompoundOccurrenceFactory()
        parent_co_2 = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co_1, parent_co_2)

        with self.subTest("2D COs, no basechem generated COs"):
            qs = collection.get_cos_for_analysis(ESP)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("2D COs w/ basechem generated COs"):
            child_co_1 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            child_co_2 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ESP
            )
            collection.compound_occurrences.add(child_co_1, child_co_2)

            qs = collection.get_cos_for_analysis(ESP)
            expected_pks = sorted([parent_co_2.pk, child_co_1.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("3D CO w/ basechem generated COs"):
            parent_co_1.molblock = "Not empty"
            parent_co_1.save()

            qs = collection.get_cos_for_analysis(ESP)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk, child_co_1.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

    def test_get_cos_for_analysis_torsion(self):
        """
        Test `get_cos_for_analysis` when the analysis is TORSION
        """
        collection = CollectionFactory()
        parent_co_1 = CompoundOccurrenceFactory()
        parent_co_2 = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(parent_co_1, parent_co_2)

        with self.subTest("2D COs, no basechem generated COs"):
            qs = collection.get_cos_for_analysis(TORSION)
            expected_pks = sorted([parent_co_1.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("2D COs w/ basechem generated COs"):
            child_co_1 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=TORSION
            )
            child_co_2 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            child_co_3 = CompoundOccurrenceFactory(
                parent_co=parent_co_1, saved_from=ALIGN
            )
            collection.compound_occurrences.add(child_co_1, child_co_2, child_co_3)

            qs = collection.get_cos_for_analysis(TORSION)
            expected_pks = sorted([child_co_2.pk, child_co_3.pk, parent_co_2.pk])
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("3D CO w/ basechem generated COs"):
            parent_co_1.molblock = "Not empty"
            parent_co_1.save()

            qs = collection.get_cos_for_analysis(TORSION)
            expected_pks = sorted(
                [parent_co_1.pk, child_co_2.pk, child_co_3.pk, parent_co_2.pk]
            )
            qs_pks = sorted([co.pk for co in qs])
            self.assertEqual(qs_pks, expected_pks)

    def test_get_cos_for_analysis_mmp(self):
        """
        Test `get_cos_for_analysis` when the analysis is MMP
        """
        collection = CollectionFactory()
        co1 = CompoundOccurrenceFactory()
        co2 = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(co1, co2)

        with self.subTest("No MMP metadata"):
            qs_pks = [co.pk for co in collection.get_cos_for_analysis(MMP)]
            self.assertEqual(len(qs_pks), 1)
            # Order is not guaranteed, a CompoundOccurrence in the mmp_analysis metadata is chosen at random
            self.assertIn(qs_pks[0], [co1.pk, co2.pk])

        with self.subTest("Has MMP metadata"):
            collection.metadata["mmp_analysis"] = {str(co2.compound.pk): {}}
            expected_pks = [co2.pk]
            qs_pks = [co.pk for co in collection.get_cos_for_analysis(MMP)]
            self.assertEqual(qs_pks, expected_pks)

        with self.subTest("Multiple COs with MMP metadata"):
            collection.metadata["mmp_analysis"][str(co1.compound.pk)] = {}
            qs_pks = [co.pk for co in collection.get_cos_for_analysis(MMP)]
            self.assertEqual(len(qs_pks), 1)
            # Order is not guaranteed, a CompoundOccurrence in the mmp_analysis metadata is chosen at random
            self.assertIn(qs_pks[0], [co1.pk, co2.pk])

    def test_get_co_pks_used_in_group(self):
        """
        Test `get_co_pks_used_in_group` returns only the pks of the CompoundOccurrences that
        were included in a task group
        """
        collection = CollectionFactory()
        co_1 = CompoundOccurrenceFactory()
        co_2 = CompoundOccurrenceFactory()
        co_3 = CompoundOccurrenceFactory(parent_co=co_1, saved_from=ALIGN)
        collection.compound_occurrences.add(co_1, co_2, co_3)

        dihedral_atoms = "1,2,3,4"
        group_name = collection.get_torsion_group_name(co_1.pk, dihedral_atoms)
        # Run fast sleep tasks with the same names as torsion tasks to cut down on testing time
        for co in [co_1, co_2, co_3]:
            task_name = co.get_torsion_task_name("SMARTS", dihedral_atoms)
            async_task(time.sleep, 0.01, task_name=task_name, group=group_name)

        co_4 = CompoundOccurrenceFactory(parent_co=co_1, saved_from=ALIGN)

        pks = sorted(collection.get_co_pks_used_in_group(group_name))
        expected_pks = sorted([co_1.pk, co_2.pk, co_3.pk])
        self.assertEqual(pks, expected_pks)

    def test_handle_sdf_upload(self):
        """
        Tests sdf files are parsed into mols
        """
        new_collection = CollectionFactory(owner=self.owner, project=self.project1)
        with self.subTest("Test uploads that are correctly defined as 2D"):
            mols = new_collection.handle_sdf_upload(self.one_uploaded_file)
            self.assertEqual(len(mols), 1)
            self.assertTrue(type(mols[0]), Chem.Mol)
            self.assertTrue(mols[0][1])  # TwoD boolean

            mols = new_collection.handle_sdf_upload(self.three_uploaded_file)
            self.assertEqual(len(mols), 3)
            self.assertTrue(type(mols[0]), Chem.Mol)
            self.assertTrue(mols[0][1])  # TwoD boolean

        with self.subTest("Test uploads that are correctly defined as 3D"):
            mols = new_collection.handle_sdf_upload(self.one_3d_uploaded_file)
            self.assertEqual(len(mols), 1)
            self.assertTrue(type(mols[0]), Chem.Mol)
            self.assertFalse(mols[0][1])  # TwoD boolean

        with self.subTest("Test uploads that are correctly defined as 2D"):
            mols = new_collection.handle_sdf_upload(self.one_2d_uploaded_file)
            self.assertEqual(len(mols), 1)
            self.assertTrue(type(mols[0]), Chem.Mol)
            self.assertTrue(mols[0][1])  # TwoD boolean

    def test_handle_romols(self):
        """
        Tests sdf uploaded files are handled to create compound occurrences
        """
        new_collection = CollectionFactory(owner=self.owner, project=self.project1)
        with self.subTest(
            "Handle file with one compound in 2D that is labeled as 2D in the sdf"
        ):
            one_romol = new_collection.handle_sdf_upload(self.one_uploaded_file)
            original_cos = new_collection.compound_occurrences.all()
            self.assertEqual(len(original_cos), 0)
            new_collection.handle_romols(one_romol, test=True)
            self.assertEqual(len(new_collection.compound_occurrences.all()), 1)
            self.assertEqual(new_collection.compound_occurrences.first().molblock, "")

        with self.subTest("Handle file with multiple compounds"):
            original_cos = new_collection.compound_occurrences.all()
            self.assertEqual(len(original_cos), 1)
            three_romols = new_collection.handle_sdf_upload(self.three_uploaded_file)
            new_collection.handle_romols(three_romols, test=True)
            # This test returns three CO's because one of the compound occurrences is the same
            self.assertEqual(len(new_collection.compound_occurrences.all()), 3)
            self.assertEqual(new_collection.compound_occurrences.first().molblock, "")

        with self.subTest("Handle file with 3D conformations"):
            self.assertEqual(len(new_collection.compound_occurrences.all()), 3)
            one_romol = new_collection.handle_sdf_upload(self.one_3d_uploaded_file)
            new_collection.handle_romols(one_romol, test=True)
            self.assertEqual(len(new_collection.compound_occurrences.all()), 4)
            self.assertEqual(new_collection.compound_occurrences.first().molblock, "")
            self.assertNotEqual(new_collection.compound_occurrences.last().molblock, "")
            self.assertEqual(len(new_collection.compounds()), 3)

        with self.subTest("Handle file with 2D conformation of the same 3D compound"):
            one_romol = self.collection.handle_sdf_upload(self.one_2d_uploaded_file)
            new_collection.handle_romols(one_romol, test=True)
            self.assertEqual(len(new_collection.compound_occurrences.all()), 4)
            self.assertEqual(len(new_collection.compounds()), 3)

    def test_get_sdf_file(self):
        """
        Tests `get_sdf_file`
        """
        with self.subTest("SDF doesn't exist"):
            filepath, filename = self.collection_one_cpd.get_sdf_file()
            self.assertTrue(os.path.exists(filepath))
            with open(filepath, "r") as f:
                r = f.read()
                self.assertNotEqual(r, "")
                self.assertIn("RDKit", r)
            self.assertIn(str(self.collection_one_cpd.pk), filename)
            self.assertIn(".sdf", filename)

        with self.subTest("SDF already exists"):
            self.assertNotEqual(self.collection_one_cpd.sdf_file, "")
            filepath, filename = self.collection_one_cpd.get_sdf_file()
            self.assertTrue(os.path.exists(filepath))
            with open(filepath, "r") as f:
                r = f.read()
                self.assertNotEqual(r, "")
                self.assertIn("RDKit", r)
            self.assertIn(str(self.collection_one_cpd.pk), filename)
            self.assertIn(".sdf", filename)

    def test_most_relevant_series(self):
        """
        Tests `most_relevant_series`
        """
        new_collection = CollectionFactory(owner=self.owner, project=self.project1)
        three_mols = new_collection.handle_sdf_upload(self.three_uploaded_file)
        new_collection.handle_romols(three_mols, test=True)
        new_collection.save()
        # Assign series 1 to two compounds and series 2 to one
        for i, comp in enumerate(new_collection.compounds()):
            comp.series = self.series1
            if i == 2:
                comp.series = self.series2
            comp.save()

        s1_reference = f"s-{self.series1.pk}"
        s2_reference = f"s-{self.series2.pk}"
        with self.subTest("valid reference"):
            valid_series = [s2_reference, s1_reference]
            result = new_collection.most_relevant_series(s1_reference, valid_series)
            # series 1 is the given reference and is valid, so pick it
            self.assertEqual(result, s1_reference)

        with self.subTest("invalid reference"):
            valid_series = [s2_reference]
            result = new_collection.most_relevant_series(s1_reference, valid_series)
            # series 1 is the given reference but is not valid, so pick the first valid series
            self.assertEqual(result, s2_reference)

        with self.subTest("valid default reference"):
            valid_series = [s2_reference, s1_reference]
            result = new_collection.most_relevant_series("default", valid_series)
            # series 1 is valid and is the assigned series for the majority of compounds
            self.assertEqual(result, s1_reference)

        with self.subTest("invalid default reference"):
            valid_series = [s2_reference]
            result = new_collection.most_relevant_series("default", valid_series)
            # series 1 is the assigned series for the majority of compounds, but is not valid
            self.assertEqual(result, s2_reference)

    @tag("local", "dtx")
    def test_create_compounds(self):
        """
        Tests compounds are created from the collection
        """
        new_collection = CollectionFactory(owner=self.owner, project=self.project1)
        one_romol = new_collection.handle_sdf_upload(self.one_uploaded_file)
        one_comp = new_collection._create_compound(one_romol[0][0], test=False)
        self.assertIsInstance(one_comp, Compound)

        three_romols = new_collection.handle_sdf_upload(self.three_uploaded_file)
        three_comp = []
        for mol, _ in three_romols:
            comp = new_collection._create_compound(mol, test=False)
            three_comp.append(comp)

        self.assertEqual(len(three_comp), 3)
        for c in three_comp:
            self.assertIsInstance(c, Compound)

        with self.subTest(
            "Test the second time a compound is seen, the project gets updated if it is AUTO"
        ):
            for c in three_comp:
                c.project = ProjectFactory(code=AUTO)
                c.series = None
                c.save()

            collection_second_upload = CollectionFactory(
                owner=self.user, project=new_collection.project
            )
            three_romols = collection_second_upload.handle_sdf_upload(
                self.three_uploaded_file
            )
            three_comp = []
            for mol, _ in three_romols:
                c = collection_second_upload._create_compound(mol, test=False)
                three_comp.append(c)
            for c in three_comp:
                self.assertIsInstance(c, Compound)
                self.assertIsNotNone(c.series)
                self.assertEqual(c.project, collection_second_upload.project)

    def test_generate_file(self):
        """
        Test result file for download is generated with content
        """
        with self.subTest("Check no file is generated"):
            # This tests the edge case of homepage being the current_view
            # that will never occur in real use,
            filepath, _ = self.collection.generate_file("homepage", test=True)
            self.assertFalse(os.path.exists(filepath))

        with self.subTest("Check romols are returned in the file"):
            filepath, filename = self.collection.generate_file(PROPCALC, test=True)
            self.assertTrue(os.path.exists(filepath))
            with open(filepath, "r") as f:
                r = f.read()
                self.assertNotEqual(r, "")

    def test_run_analysis(self):
        """
        Test run analysis executes the correct analysis
        """
        with self.subTest("Test propcalc analysis"):
            co_ids = [
                co.id for co in self.collection_one_cpd.compound_occurrences.all()
            ]
            self.collection_one_cpd.metadata["props_to_show"] = [MW, CLOGP]
            group_name = self.collection_one_cpd.get_propcalc_group_name()
            self.collection_one_cpd.run_analysis(PROPCALC)
            failed, completed, property_results = get_group_results(group_name)
            self.assertFalse(failed)
            self.assertTrue(completed)
            self.assertEqual(len(property_results), 1)
            self.assertListEqual(list(property_results.keys()), co_ids)
            for co_id in co_ids:
                v = property_results[co_id]
                self.assertIn("mw", v)

        with self.subTest("Test align analysis"):
            co_ids = []
            for co in self.collection_one_cpd.compound_occurrences.all():
                co_ids.append(co.id)
                co.compound.series = self.series1
                co.compound.save()

            ref_string = f"s-{self.series1.pk}"
            group_name = self.collection_one_cpd.get_align_group_name(ref_string)

            self.collection_one_cpd.run_analysis(ALIGN, ref_string=ref_string)
            failed, completed, result = get_group_results(group_name)
            self.assertFalse(failed)
            self.assertTrue(completed)
            self.assertEqual(len(result), 3)
            self.assertIn(ref_string, list(result["references"].keys()))
            self.assertIn(ref_string, list(result["receptors"].keys()))
            self.assertIn(f"co-{co_ids[0]}", list(result["compounds"].keys()))

        with self.subTest("Test dock analysis"):
            co_ids = []
            for co in self.collection_one_cpd.compound_occurrences.all():
                co_ids.append(co.id)
                co.compound.series = self.series1
                co.compound.save()

            ref_string = f"s-{self.series1.pk}"
            group_name = self.collection_one_cpd.get_dock_group_name(ref_string)

            self.collection_one_cpd.run_analysis(DOCK, ref_string=ref_string)
            failed, completed, result = get_group_results(group_name)
            self.assertFalse(failed)
            self.assertTrue(completed)
            self.assertEqual(len(result), 3)
            self.assertIn(ref_string, list(result["references"].keys()))
            self.assertIn(f"co-{co_ids[0]}", list(result["compounds"].keys()))
            self.assertIn(ref_string, list(result["receptors"].keys()))

    @tag("mmpdb")
    def test_find_mmps(self):
        """
        Tests MMPS are added to this collection correctly
        """
        for c in self.collection_mmp.compounds():
            self.assertFalse(c.mmps.exists())
        self.collection_mmp.find_mmps()
        compounds = self.collection_mmp.compounds()
        self.assertEqual(compounds[0].mmps.all().count(), 15)
        self.assertEqual(compounds[0].mmps.exclude(dn_id="").count(), 10)
        self.assertEqual(compounds[1].mmps.all().count(), 14)
        self.assertEqual(compounds[1].mmps.exclude(dn_id="").count(), 10)
        self.assertEqual(compounds[2].mmps.all().count(), 7)
        self.assertEqual(compounds[2].mmps.exclude(dn_id="").count(), 1)

    @tag("mmpdb")
    def test_update_mmp_dtx_avg_assay_data(self):
        """
        Test `update_mmp_dtx_avg_assay_data` pulls updated assay data for all Compounds and
        their MMPs in the given collection
        """
        mock = self.collection_get_agg_ic50_data_mock
        mock.return_value = {
            "DN0000001": {
                "data": {"assay 1": {"analysis 1": {"RESULT_GEOM_MEAN": 10}}}
            },
            "DN0000002": {
                "data": {"assay 1": {"analysis 1": {"RESULT_GEOM_MEAN": 11}}}
            },
        }
        # Create collection with two COs, one with an mmp and one without
        mmp = CompoundFactory(dn_id="DN0000002")
        c_w_mmps = CompoundFactory(dn_id="DN0000001")
        c_w_mmps.mmps.set([mmp])
        c_wo_mmps = CompoundFactory(dn_id="DN0000003")
        co_w_mmps = CompoundOccurrenceFactory(compound=c_w_mmps)
        co_wo_mmps = CompoundOccurrenceFactory(compound=c_wo_mmps)
        collection = CollectionFactory()
        collection.compound_occurrences.set([co_w_mmps, co_wo_mmps])

        self.assertEqual(c_wo_mmps.measured_data, {})
        self.assertEqual(c_w_mmps.measured_data, {})
        self.assertEqual(mmp.measured_data, {})
        collection.update_mmp_dtx_avg_assay_data()
        c_wo_mmps.refresh_from_db()
        self.assertEqual(c_wo_mmps.measured_data, {"assay_results": {}})
        c_w_mmps.refresh_from_db()
        expected = {"assay_results": {"assay 1": {"analysis 1": 10}}}
        self.assertEqual(c_w_mmps.measured_data, expected)
        mmp.refresh_from_db()
        expected = {"assay_results": {"assay 1": {"analysis 1": 11}}}
        self.assertEqual(mmp.measured_data, expected)

    @override_settings(INDUCTIVE_BIO_ENABLED=True)
    def test_inductive_in_props_ENABLED(self):
        """
        Test that `inductive_in_props` correctly identifies if any of the selected properties require
        InductiveBio
        """
        with self.subTest("No props metadata"):
            self.collection.metadata = {}
            self.assertFalse(self.collection.inductive_in_props())
        with self.subTest("No inductive props"):
            self.collection.metadata["props_to_show"] = [MW, CLOGP]
            self.assertFalse(self.collection.inductive_in_props())
        with self.subTest("Inductive props"):
            self.collection.metadata["props_to_show"] = [MW, CLOGP, ALOGD]
            self.assertTrue(self.collection.inductive_in_props())

    @override_settings(INDUCTIVE_BIO_ENABLED=False)
    def test_inductive_in_props_NOT_ENABLED(self):
        """
        Test that `inductive_in_props` correctly returns False regardless of the selected properties
        if InductiveBio is not enabled.
        """
        with self.subTest("No props metadata"):
            self.collection.metadata = {}
            self.assertFalse(self.collection.inductive_in_props())
        with self.subTest("No inductive props"):
            self.collection.metadata["props_to_show"] = [MW, CLOGP]
            self.assertFalse(self.collection.inductive_in_props())
        with self.subTest("Inductive props"):
            self.collection.metadata["props_to_show"] = [MW, CLOGP, ALOGD]
            self.assertFalse(self.collection.inductive_in_props())

    @tag("inductive", "external")
    def test_propcalc_analysis(self):
        """
        Test propcalc analysis returns expected property dict
        """
        with self.subTest("Run aLogD prediction"):
            cos = self.collection.compound_occurrences.all()
            co_ids = [co.id for co in cos]
            self.collection.metadata["props_to_show"] = [MW, CLOGP, ALOGD]
            self.collection.save()
            self.collection.propcalc_analysis()

            group_name = self.collection.get_propcalc_group_name()
            _, _, prop_dict = get_group_results(group_name)

            self.assertEqual(len(prop_dict), 3)
            self.assertListEqual(sorted(list(prop_dict.keys())), sorted(co_ids))

            for co in cos:
                v = prop_dict[co.pk]
                self.assertIn("mw", v)
                self.assertIn("logd_prediction", v)
                self.assertNotEqual(co.compound.predict_logd(), None)

        with self.subTest("Run LM predictions"):
            cos = self.collection_one_cpd.compound_occurrences.all()
            co_ids = [co.id for co in cos]
            self.collection_one_cpd.metadata["props_to_show"] = [MW, CLOGP, RLM, HLM]
            self.collection_one_cpd.save()
            self.collection_one_cpd.propcalc_analysis()

            group_name = self.collection_one_cpd.get_propcalc_group_name()
            _, _, prop_dict = get_group_results(group_name)

            self.assertEqual(len(prop_dict), 1)
            self.assertListEqual(sorted(list(prop_dict.keys())), sorted(co_ids))

            v = prop_dict[co_ids[0]]
            # Check RLM prediction has both prediction and probabilities image
            self.assertTrue(float(v["rlm_prediction"]) > 39)
            self.assertTrue(len(v["rlm_probabilities"]) > 50)
            self.assertTrue(bool(v["rlm_ood"]))
            # Check HLM prediction has both prediction and probabilities image
            self.assertTrue(float(v["hlm_prediction"]) > 15)
            self.assertTrue(len(v["hlm_probabilities"]) > 50)
            self.assertTrue(bool(v["hlm_ood"]))

    def test_align_analysis(self):
        """
        Tests the align analysis
        """
        co_ids = []
        c_ids = []
        for co in self.collection_one_cpd.compound_occurrences.all():
            co_ids.append(co.id)
            c_ids.append(co.compound.id)
            co.compound.series = self.series1
            co.compound.save()

        ref_string = f"s-{self.series1.pk}"
        self.collection_one_cpd.align_analysis(ref_string)

        group_name = self.collection_one_cpd.get_align_group_name(ref_string)
        _, _, result = get_group_results(group_name)

        # Check Series exists in result dict
        self.assertEqual(len(result), 3)
        series_key = f"s-{self.series1.pk}"
        self.assertIn(series_key, list(result["references"].keys()))
        self.assertIn(self.series1.dn_id, result["references"][series_key])
        self.assertIn(series_key, list(result["receptors"].keys()))

        # Check Comp has conformers and is in result dict
        confs_dict = result["compounds"][f"co-{co_ids[0]}"]

        self.assertEqual(len(confs_dict), 2)
        expected_key1 = f"c{c_ids[0]}-co{co_ids[0]}-1-{ref_string}"
        self.assertIn(expected_key1, list(confs_dict.keys()))
        self.assertIn("VSM-", confs_dict[expected_key1]["moltext"])
        self.assertIn("0.000", confs_dict[expected_key1]["r_bc_rmsd_to_lsalign"])
        self.assertIn("1.696", confs_dict[expected_key1]["r_mmff_rel_energy"])

        expected_key2 = f"c{c_ids[0]}-co{co_ids[0]}-2-{ref_string}"
        self.assertIn(expected_key2, list(confs_dict.keys()))
        self.assertIn("VSM-", confs_dict[expected_key2]["moltext"])
        self.assertIn("1.032", confs_dict[expected_key2]["r_bc_rmsd_to_lsalign"])
        self.assertIn("0.000", confs_dict[expected_key2]["r_mmff_rel_energy"])

    def test_send_task_notifications(self):
        """
        Test that `send_task_notifications` sends an email with the appropriate info
        """
        with self.subTest("Send no email for short tasks"):
            self.assertEqual(len(mail.outbox), 0)
            self.collection_one_cpd.metadata["props_to_show"] = [MW, CLOGP]
            self.collection_one_cpd.save()
            self.assertEqual(len(mail.outbox), 0)

        with self.subTest("Check emails for long tasks"):
            self.assertEqual(len(mail.outbox), 0)
            ref_string = f"s-{self.series1.pk}"
            # Run an alignment task. This may or may not send emails depending on how long
            # it takes to run
            self.collection.align_analysis(ref_string)
            original_outbox_length = len(mail.outbox)
            # "Hack" the tasks to send emails
            group_name = self.collection.get_align_group_name(ref_string)
            tasks = fetch_group(group_name, failures=True)
            failed_recipients = [email for _, email in settings.ADMINS]
            for i, t in enumerate(tasks, 1):
                # Make sure each task takes at least 60 seconds
                t.started = datetime.datetime.now()
                t.stopped = t.started + datetime.timedelta(0, 65)
                t.success = False
                t.save()  # This calls send_task_notifications as part of the hook

                # Sends two emails, one is an error email sent to django admins because the
                # task failed, and one is a success/failure notification to the requestor
                # because all tasks are complete
                self.assertEqual(len(mail.outbox), original_outbox_length + 2 * i)
                # All failed tasks send an admin failure notification
                self.assertEqual(mail.outbox[-2].subject, "[Django] Admin Failure")
                self.assertIn(str(round(t.time_taken())), mail.outbox[-2].body)
                self.assertIn(str(self.collection.owner), mail.outbox[-2].body)
                self.assertIn(t.name, mail.outbox[-2].body)
                self.assertEqual(sorted(mail.outbox[-2].to), sorted(failed_recipients))
                if i == 4:
                    # Last task sends a failed notification since all tasks failed
                    self.assertEqual(mail.outbox[-1].subject, "Basechem Task Failed")
                    failed_recipients.append(self.collection.owner.email)
                    self.assertEqual(
                        sorted(mail.outbox[-1].to), sorted(failed_recipients)
                    )
                else:
                    # First 3 tasks send a complete notification since at least one task succeeded
                    self.assertEqual(mail.outbox[-1].subject, "Basechem Task Complete")
                    self.assertEqual(mail.outbox[-1].to, [self.collection.owner.email])

    def test_get_url(self):
        """
        Test that `get_url` gets the full url
        """
        full_url = self.collection.get_url(PROPCALC)
        url_suffix = reverse(PROPCALC, kwargs={"collection_id": self.collection.id})
        self.assertEqual(full_url, settings.BASE_URL + url_suffix)
        self.assertIn("http", settings.BASE_URL)

    def test_get_group_name(self):
        """
        Tests appropriate group names are returned
        """
        with self.subTest("propcalc"):
            expected_name = f"propcalc_{self.collection_one_cpd.id}"
            self.assertEqual(
                self.collection_one_cpd.get_propcalc_group_name(), expected_name
            )
        with self.subTest("align"):
            expected_name = f"align_{self.collection_one_cpd.id}_default"
            self.assertEqual(
                self.collection_one_cpd.get_align_group_name("default"), expected_name
            )
            expected_name = f"align_{self.collection_one_cpd.id}_s-1"
            self.assertEqual(
                self.collection_one_cpd.get_align_group_name("s-1"), expected_name
            )
        with self.subTest("dock"):
            expected_name = f"dock_{self.collection_one_cpd.id}_default"
            self.assertEqual(
                self.collection_one_cpd.get_dock_group_name("default"), expected_name
            )
            expected_name = f"dock_{self.collection_one_cpd.id}_s-1"
            self.assertEqual(
                self.collection_one_cpd.get_dock_group_name("s-1"), expected_name
            )

        with self.subTest("esp"):
            expected_name = f"esp_{self.collection_one_cpd.id}_default"
            self.assertEqual(
                self.collection_one_cpd.get_esp_group_name("default"), expected_name
            )
            expected_name = f"esp_{self.collection_one_cpd.id}_s-1"
            self.assertEqual(
                self.collection_one_cpd.get_esp_group_name("s-1"), expected_name
            )

        with self.subTest("torsion"):
            comp = self.collection_torsion.compounds()[0]
            expected_name = f"torsion_{self.collection_torsion.pk}_{comp.pk}_1-2-3-4"
            self.assertEqual(
                self.collection_torsion.get_torsion_group_name(comp.pk, "1,2,3,4"),
                expected_name,
            )

        with self.subTest("mmp"):
            comp = self.collection_mmp.compounds()[2]
            expected_name = f"mmp_{self.collection_mmp.pk}_{comp.pk}"
            self.assertEqual(
                self.collection_mmp.get_mmp_group_name(comp.pk), expected_name
            )

    def test_dock_analysis(self):
        """
        Tests the dock analysis
        """
        co_ids = []
        for co in self.collection_one_cpd.compound_occurrences.all():
            co_ids.append(co.id)
            co.compound.series = self.series1
            co.compound.save()

        ref_string = f"s-{self.series1.pk}"
        self.collection_one_cpd.dock_analysis(ref_string)

        group_name = self.collection_one_cpd.get_dock_group_name(ref_string)
        _, _, result = get_group_results(group_name)

        # Check result dict
        self.assertEqual(len(result), 3)

        self.assertEqual(
            sorted(["receptors", "references", "compounds"]),
            sorted(list(result.keys())),
        )
        self.assertIn(f"s-{self.series1.pk}", list(result["references"].keys()))
        for co_id in co_ids:
            # Dock results are not exactly the same every run, sometimes the best rmsd
            # and best score poses overlap to result in 3/4 results
            self.assertTrue(3 <= len(result["compounds"][f"co-{co_id}"]) <= 8)

    def test_collect_dock_references(self):
        """
        Tests collecting references for the dock analysis
        """
        for c in self.collection_one_cpd.compounds():
            c.series = self.series1
            c.save()

        result = self.collection_one_cpd._collect_dock_references()

        self.assertEqual(len(result), 2)
        self.assertIn(f"s-{self.series1.pk}", list(result["references"].keys()))
        self.assertIn(f"s-{self.series1.pk}", list(result["receptors"].keys()))
        # There should only be one receptor b/c one compound/series pair
        self.assertEqual(len(result["receptors"]), 1)

    @tag("local", "esp")
    def test_esp_analysis(self):
        """
        Tests the ESP analysis
        """
        co_ids = []
        for co in self.collection_one_cpd.compound_occurrences.all():
            co_ids.append(co.id)
            co.compound.series = self.series1
            co.compound.save()

        ref_string = f"s-{self.series1.pk}"
        self.collection_one_cpd.esp_analysis(ref_string)

        group_name = self.collection_one_cpd.get_esp_group_name(ref_string)
        _, _, result = get_group_results(group_name)

        # Check result dict
        self.assertEqual(len(result), 3)
        self.assertEqual(
            sorted(["receptors", "references", "compounds"]),
            sorted(list(result.keys())),
        )
        self.assertIn(f"s-{self.series1.pk}", list(result["references"].keys()))
        self.assertIn(f"s-{self.series1.pk}", list(result["receptors"].keys()))
        for co_id in co_ids:
            self.assertEqual(len(result["compounds"][f"co-{co_id}"]), 3)

    @tag("local", "esp")
    def test_collect_esp_references(self):
        """
        Tests collecting references for the ESP analysis
        """
        co_ids = []
        for co in self.collection_one_cpd.compound_occurrences.all():
            co_ids.append(co.id)
            co.compound.series = self.series1
            co.compound.save()

        result = self.collection_one_cpd._collect_esp_references()
        self.assertEqual(len(result), 2)
        self.assertIn(f"s-{self.series1.pk}", list(result["references"].keys()))
        self.assertIn(f"s-{self.series1.pk}", list(result["receptors"].keys()))

    @tag("mmpdb")
    def test_mmp_analysis(self):
        """
        Test that `Collection.mmp_analysis` runs `Compound.mmp_analysis` for all
        Compounds in the given Collection
        """
        co = self.collection_mmp.get_cos_for_analysis(MMP)[0]
        self.assertEqual(co.compound.mmps.all().count(), 0)
        self.collection_mmp.mmp_analysis()

        co.compound.refresh_from_db()
        self.assertNotEqual(co.compound.mmps.all().count(), 0)

        group_name = self.collection_mmp.get_mmp_group_name(co.compound.pk)
        _, _, result = get_group_results(group_name)
        self.assertEqual(list(result.keys()), [f"co-{co.pk}"])
        self.assertEqual(len(result[f"co-{co.pk}"]), 2)
