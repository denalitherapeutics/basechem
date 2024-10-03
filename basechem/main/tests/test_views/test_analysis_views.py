from unittest.mock import patch

from django.test import tag
from django.urls import reverse

from basechem.common.tests.base import BasechemTestCase, BasechemViewTestMixin
from basechem.main.constants import (
    ACCEPTORS,
    ALIGN,
    CLOGP,
    DOCK,
    DTX_MMP,
    ESP,
    HLM,
    MMP,
    MW,
    PROPCALC,
    RLM,
    TORSION,
)
from basechem.main.forms import ChooseReferenceForm
from basechem.main.tasks import get_analysis_results, get_group_results
from basechem.main.tests.factories import SeriesFactory


class PropCalcTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.url_one_cpd = reverse(
            PROPCALC, kwargs={"collection_id": cls.collection_one_cpd.id}
        )
        cls.url_no_cpds = reverse(
            PROPCALC, kwargs={"collection_id": cls.collection_no_cpds.id}
        )
        cls.collection_one_cpd.metadata = {"props_to_show": [ACCEPTORS, MW, CLOGP]}
        cls.collection_one_cpd.save()

        cls.collection_no_cpds.metadata = {"props_to_show": [MW, CLOGP]}
        cls.collection_no_cpds.save()

    def test_unauthenticated(self):
        """
        User must be logged in order to see the propcalc view
        """
        self.client.logout()
        response = self.client.get(self.url_one_cpd)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_one_cpd}")
        self.assertEqual(self.client.session["next"], self.url_one_cpd)
        response = self.client.post(self.url_one_cpd)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_one_cpd}")

    def test_propcalc_view(self):
        """
        Test logged in user can see the propcalc view
        """
        response = self.client.get(self.url_one_cpd, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/propcalc/propcalc.html")

        response = self.client.get(self.url_no_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/propcalc/propcalc.html")

    def test_get_context_data(self):
        """
        Tests context data includes the properties when requested
        """
        with self.subTest("Test context data for collection with no compounds"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(response.context["current_view"], PROPCALC)
            self.assertIn("collection", response.context)
            self.assertIn("props", response.context)
            self.assertEqual(response.context["props"], [MW, CLOGP])
            self.assertEqual(response.context["property_results"], {})

        with self.subTest("Test context data for collection with compounds"):
            response = self.client.get(self.url_one_cpd)
            self.assertEqual(response.context["current_view"], PROPCALC)
            self.assertIn("collection", response.context)
            self.assertIn("props", response.context)
            self.assertEqual(response.context["props"], [ACCEPTORS, MW, CLOGP])
            group_name = self.collection_one_cpd.get_propcalc_group_name()
            failed, completed, result = get_group_results(group_name)
            self.assertFalse(failed)
            self.assertTrue(completed)
            self.assertEqual(len(result), 1)

    @tag("inductive", "external")
    def test_get_context_data_rlm(self):
        """
        Tests context data includes the rlm info when requested
        """
        self.collection_one_cpd.metadata["props_to_show"] = ["MW", "RLM"]
        self.collection_one_cpd.save()

        response = self.client.get(self.url_one_cpd)
        self.assertEqual(response.context["current_view"], PROPCALC)
        self.assertIn("collection", response.context)
        self.assertIn("props", response.context)
        self.assertListEqual(
            response.context["props"], [MW, "RLM Prediction", "RLM Probabilities"]
        )
        group_name = self.collection_one_cpd.get_propcalc_group_name()
        failed, completed, result = get_group_results(group_name)
        self.assertFalse(failed)
        self.assertTrue(completed)
        self.assertEqual(len(result), 1)
        key = self.collection_one_cpd.compound_occurrences.first().pk
        self.assertTrue(float(result[key]["rlm_prediction"]) > 39)
        self.assertTrue(len(result[key]["rlm_probabilities"]) > 50)
        self.assertEqual("True", result[key]["rlm_ood"])

    @tag("inductive", "external")
    def test_get_context_data_lm(self):
        """
        Tests context data includes rlm and hlm info when requested
        """
        self.collection_one_cpd.metadata["props_to_show"] = [MW, RLM, HLM]
        self.collection_one_cpd.save()

        response = self.client.get(self.url_one_cpd)
        self.assertEqual(response.context["current_view"], PROPCALC)
        self.assertIn("collection", response.context)
        self.assertIn("props", response.context)
        self.assertListEqual(
            response.context["props"],
            [
                MW,
                "RLM Prediction",
                "RLM Probabilities",
                "HLM Prediction",
                "HLM Probabilities",
            ],
        )
        group_name = self.collection_one_cpd.get_propcalc_group_name()
        failed, completed, result = get_group_results(group_name)
        self.assertFalse(failed)
        self.assertTrue(completed)
        self.assertEqual(len(result), 1)
        key = self.collection_one_cpd.compound_occurrences.first().pk
        self.assertTrue(result[key]["rlm_prediction"] > 15)
        self.assertTrue(len(result[key]["hlm_probabilities"]) > 50)
        self.assertEqual("True", result[key]["hlm_ood"])

    def test_get_queryset(self):
        """
        Tests queryset is correct
        """
        with self.subTest("Test empty queryset"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(len(response.context["compounds"]), 0)

        with self.subTest("Test non-empty queryset"):
            response = self.client.get(self.url_one_cpd)
            self.assertEqual(len(response.context["compounds"]), 1)


class LigandAlignTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.url_three_cpds = reverse(ALIGN, kwargs={"collection_id": cls.collection.id})
        cls.url_no_cpds = reverse(
            ALIGN, kwargs={"collection_id": cls.collection_no_cpds.id}
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to see the ligand align view
        """
        self.client.logout()
        response = self.client.get(self.url_three_cpds)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_three_cpds}")
        self.assertEqual(self.client.session["next"], self.url_three_cpds)
        response = self.client.post(self.url_three_cpds)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_three_cpds}")

    def test_align_view(self):
        """
        Test logged in user can see the ligand align view
        """
        response = self.client.get(self.url_three_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/align/align.html")

        response = self.client.get(self.url_no_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/align/align.html")

    def test_context_data_before_submit(self):
        """
        Tests context data before the user selects compounds for alignment
        """
        with self.subTest("No compounds"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(response.context["current_view"], ALIGN)
            self.assertEqual(type(response.context["form"]), ChooseReferenceForm)
            self.assertEqual(response.context["form"].instance, self.collection_no_cpds)

        with self.subTest("Three compounds"):
            response = self.client.get(self.url_three_cpds)
            self.assertEqual(response.context["current_view"], ALIGN)
            self.assertEqual(type(response.context["form"]), ChooseReferenceForm)
            self.assertEqual(response.context["form"].instance, self.collection)

    def test_context_data_after_submit(self):
        """
        Tests context data after the user selects compounds for alignment
        """
        c_ids = [c.id for c in self.collection.compounds()]
        # Assign series to compounds
        for c in self.collection.compounds():
            c.series = self.series1
            c.save()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "compounds": c_ids,
            "reference": "default",
        }
        json_response = self.client.post(self.url_three_cpds, data=form_data)
        redirect_url = json_response.json()["redirect_url"]
        response = self.client.get(redirect_url)
        self.assertEqual(response.context["current_view"], ALIGN)

        # Check that tasks ran
        group_name = self.collection.get_align_group_name(form_data["reference"])

        failed, _, result = get_group_results(group_name)

        self.assertFalse(failed)
        # Top level "references", "receptors", and "compounds" keys
        self.assertEqual(len(result), 3)
        self.assertEqual(len(result["receptors"]), 1)
        self.assertEqual(len(result["references"]), 2)
        self.assertEqual(len(result["compounds"]), 3)

        expected_url = reverse(ALIGN, args=[self.collection.id, group_name])
        self.assertEqual(redirect_url, expected_url)

    def test_get_analysis_results(self):
        """
        Tests context data after the user selects compounds for alignment
        """
        c_ids = [c.id for c in self.collection.compounds()]
        # Assign series to compounds
        for c in self.collection.compounds():
            c.series = self.series1
            c.save()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "compounds": c_ids,
            "reference": "default",
        }
        json_response = self.client.post(self.url_three_cpds, data=form_data)
        redirect_url = json_response.json()["redirect_url"]
        response = self.client.get(redirect_url)
        self.assertEqual(response.context["current_view"], ALIGN)

        # Check that tasks ran
        group_name = self.collection.get_align_group_name(form_data["reference"])
        analysis_kwargs = {"ref_string": form_data["reference"]}

        failed, result = get_analysis_results(
            self.collection, ALIGN, group_name, analysis_kwargs
        )

        self.assertFalse(failed)
        # Top level "references", "receptors", and "compounds" keys
        self.assertEqual(len(result), 3)
        self.assertEqual(len(result["receptors"]), 1)
        self.assertEqual(len(result["references"]), 2)
        self.assertEqual(len(result["compounds"]), 3)

        expected_url = reverse(ALIGN, args=[self.collection.id, group_name])
        self.assertEqual(expected_url, redirect_url)


class DockTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        for comp in cls.collection.compounds():
            comp.series = cls.series1
            comp.save()

        cls.url_three_cpds = reverse(DOCK, kwargs={"collection_id": cls.collection.id})
        cls.url_no_cpds = reverse(
            DOCK, kwargs={"collection_id": cls.collection_no_cpds.id}
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to see the dock view
        """
        self.client.logout()
        response = self.client.get(self.url_three_cpds)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_three_cpds}")
        self.assertEqual(self.client.session["next"], self.url_three_cpds)
        response = self.client.post(self.url_three_cpds)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_three_cpds}")

    def test_dock_view(self):
        """
        Test logged in user can see the dock view
        """
        response = self.client.get(self.url_three_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/dock/dock.html")

        response = self.client.get(self.url_no_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/dock/dock.html")

    def test_context_data_before_submit(self):
        """
        Tests context data before the user submits compounds for docking
        """
        with self.subTest("No compounds"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(response.context["current_view"], DOCK)
            self.assertEqual(type(response.context["form"]), ChooseReferenceForm)
            self.assertEqual(response.context["form"].instance, self.collection_no_cpds)

        with self.subTest("Three compounds"):
            response = self.client.get(self.url_three_cpds)
            self.assertEqual(response.context["current_view"], DOCK)
            self.assertEqual(type(response.context["form"]), ChooseReferenceForm)
            self.assertEqual(response.context["form"].instance, self.collection)

    def test_context_data_after_submit(self):
        """
        Tests context data after the user submits compounds for docking
        """
        # Assign series to compounds
        for c in self.collection.compounds():
            c.series = self.series1
            c.save()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "reference": "default",
        }
        json_response = self.client.post(self.url_three_cpds, data=form_data)
        redirect_url = json_response.json()["redirect_url"]
        response = self.client.get(redirect_url)
        self.assertEqual(response.context["current_view"], DOCK)

        # Check that tasks ran
        group_name = self.collection.get_dock_group_name(form_data["reference"])
        failed, _, result = get_group_results(group_name)

        self.assertFalse(failed)
        # Top level "receptors", "references", and "compounds" keys
        self.assertEqual(len(result), 3)
        self.assertIn("receptors", result)
        self.assertEqual(len(result["references"]), 1)
        self.assertEqual(len(result["compounds"]), 3)

        expected_url = reverse(DOCK, args=[self.collection.id, group_name])
        self.assertEqual(redirect_url, expected_url)

    def test_failed_post(self):
        """
        Test that a user-readable error is returned when a user tries to "default" dock
        a collection where not every CompoundOccurrence has a series w/ a mol2 file
        """
        # Assign bad series to compounds
        for i, c in enumerate(self.collection.compounds()):
            c.series = SeriesFactory(name=f"series {i}")
            c.save()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "reference": "default",
        }
        json_response = self.client.post(self.url_three_cpds, data=form_data)
        response_data = json_response.json()
        expected_errors = {
            "__all__": [
                "At least one of the assigned Series don't have the necessary files to run Dock. Please select a reference from this dropdown. Offending Series: series 0, series 1, series 2"
            ]
        }
        self.assertEqual(response_data["errors"], expected_errors)


class EspTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        for comp in cls.collection_one_cpd.compounds():
            comp.series = cls.series1
            comp.save()

        cls.url_one_cpd = reverse(
            ESP, kwargs={"collection_id": cls.collection_one_cpd.id}
        )
        cls.url_no_cpds = reverse(
            ESP, kwargs={"collection_id": cls.collection_no_cpds.id}
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to see the dock view
        """
        self.client.logout()
        response = self.client.get(self.url_one_cpd)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_one_cpd}")
        self.assertEqual(self.client.session["next"], self.url_one_cpd)
        response = self.client.post(self.url_one_cpd)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url_one_cpd}")

    @tag("local", "esp")
    def test_esp_view(self):
        """
        Test logged in user can see the ESP view
        """
        response = self.client.get(self.url_one_cpd, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/esp/esp.html")

        response = self.client.get(self.url_no_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/esp/esp.html")

    @tag("local", "esp")
    def test_context_data_before_submit(self):
        """
        Tests context data before the user submits compounds for docking
        """
        with self.subTest("No compounds"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(response.context["current_view"], ESP)
            self.assertEqual(type(response.context["form"]), ChooseReferenceForm)
            self.assertEqual(response.context["form"].instance, self.collection_no_cpds)

        with self.subTest("Three compounds"):
            response = self.client.get(self.url_one_cpd)
            self.assertEqual(response.context["current_view"], ESP)
            self.assertEqual(type(response.context["form"]), ChooseReferenceForm)
            self.assertEqual(response.context["form"].instance, self.collection_one_cpd)

    @tag("local", "esp")
    def test_context_data_after_submit(self):
        """
        Tests context data after the user submits compounds for ESP map generation
        """
        # Assign series to compounds
        for c in self.collection_one_cpd.compounds():
            c.series = self.series1
            c.save()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "reference": "default",
        }
        json_response = self.client.post(self.url_one_cpd, data=form_data)
        redirect_url = json_response.json()["redirect_url"]
        response = self.client.get(redirect_url)
        self.assertEqual(response.context["current_view"], ESP)

        # Check that tasks ran
        group_name = self.collection_one_cpd.get_esp_group_name(form_data["reference"])

        failed, _, result = get_group_results(group_name)

        self.assertFalse(failed)
        # Top level "references", "receptors", and "compounds" keys
        self.assertEqual(len(result), 3)
        self.assertEqual(len(result["receptors"]), 1)
        self.assertEqual(len(result["references"]), 2)
        self.assertEqual(len(result["compounds"]), 1)

        expected_url = reverse(ESP, args=[self.collection_one_cpd.id, group_name])
        self.assertEqual(expected_url, redirect_url)


class TorsionTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.url_three_cpds = reverse(
            TORSION, kwargs={"collection_id": cls.collection.id}
        )
        cls.url_torsion_cpd = reverse(
            TORSION, kwargs={"collection_id": cls.collection_torsion.id}
        )
        cls.url_no_cpds = reverse(
            TORSION, kwargs={"collection_id": cls.collection_no_cpds.id}
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to see the torsion view
        """
        self.client.logout()
        response = self.client.get(self.url_torsion_cpd)
        self.assertRedirects(
            response, f"{reverse('login')}?next={self.url_torsion_cpd}"
        )
        self.assertEqual(self.client.session["next"], self.url_torsion_cpd)
        response = self.client.post(self.url_torsion_cpd)
        self.assertRedirects(
            response, f"{reverse('login')}?next={self.url_torsion_cpd}"
        )

    def test_torsion_view(self):
        """
        Test logged in user can see the torsion view
        """
        response = self.client.get(self.url_torsion_cpd, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/torsion/torsion.html")

        response = self.client.get(self.url_no_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/torsion/torsion.html")

    def test_context_data_before_submit(self):
        """
        Tests context data before the user submits compounds for scanning
        """
        with self.subTest("No compounds"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(response.context["current_view"], TORSION)
            self.assertEqual(response.context["form"].instance, self.collection_no_cpds)

        with self.subTest("Three compounds"):
            response = self.client.get(self.url_torsion_cpd)
            self.assertEqual(response.context["current_view"], TORSION)
            self.assertEqual(response.context["form"].instance, self.collection_torsion)

    @patch("basechem.main.models.compound_models.run_mc_torsion_scan")
    def test_context_data_after_submit(self, mock_run_mc_torsion_scan):
        """
        Tests context data after the user submits compounds for Torsion scan generation
        """
        mock_run_mc_torsion_scan.return_value = (
            "basechem/main/tests/testdata/test_torsion_results.sdf"
        )
        # Assign series to compounds
        for c in self.collection_torsion.compounds():
            c.series = self.series1
            c.save()
        co = self.collection_torsion.compound_occurrences.all().first()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "pioneer": f"co-{co.pk}",
            "selected_atoms": "0,1,2,3",
        }
        json_response = self.client.post(self.url_torsion_cpd, data=form_data)
        redirect_url = json_response.json()["redirect_url"]
        response = self.client.get(redirect_url)
        self.assertEqual(response.context["current_view"], TORSION)

        # Check that tasks ran
        dihedral_smarts = co.convert_atoms_to_smarts("0,1,2,3")
        group_name = self.collection_torsion.get_torsion_group_name(co.pk, "0,1,2,3")

        failed, _, result = get_group_results(group_name)

        # Check returned data
        self.assertFalse(failed)
        self.assertEqual(len(result), 1)

        expected_url = reverse(TORSION, args=[self.collection_torsion.id, group_name])
        self.assertEqual(expected_url, redirect_url)

    def test_failed_post(self):
        """
        Test that a user-readable error is returned when a user tries to run torsion on a
        nonsensical dihedral (disconnected atoms)
        """
        co = self.collection.compound_occurrences.all().first()
        form_data = {
            "csrfmiddlewaretoken": "token",
            "pioneer": f"co-{co.pk}",
            "selected_atoms": "0,1,2,4",
        }
        json_response = self.client.post(self.url_three_cpds, data=form_data)
        response_data = json_response.json()
        expected_errors = {
            "__all__": [
                "The 4 atoms you selected did not result in a valid dihedral to scan. Please try again."
            ]
        }
        self.assertEqual(response_data["errors"], expected_errors)


@tag("mmpdb")
class DTXMMPTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.url_collection_mmp = reverse(
            DTX_MMP, kwargs={"collection_id": cls.collection_mmp.pk}
        )
        cls.url_collection_3_comps = reverse(
            DTX_MMP, kwargs={"collection_id": cls.collection.pk}
        )
        cls.url_no_cpds = reverse(
            DTX_MMP, kwargs={"collection_id": cls.collection_no_cpds.pk}
        )

    def test_unauthenticated(self):
        """
        User must be logged in to see the DTX MMPs view
        """
        self.client.logout()
        response = self.client.get(self.url_collection_mmp)
        self.assertRedirects(
            response, f"{reverse('login')}?next={self.url_collection_mmp}"
        )
        self.assertEqual(self.client.session["next"], self.url_collection_mmp)
        response = self.client.post(self.url_collection_mmp)
        self.assertRedirects(
            response, f"{reverse('login')}?next={self.url_collection_mmp}"
        )

    def test_dtx_mmp_view(self):
        """
        Test logged in user can see the DTX MMP view
        """
        response = self.client.get(self.url_collection_mmp, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/mmps/dtx_mmps.html")

        response = self.client.get(self.url_no_cpds, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/mmps/dtx_mmps.html")

    def test_context_data(self):
        """
        Tests the context data contains the correct Compounds
        """
        with self.subTest("Collection has no Compounds"):
            response = self.client.get(self.url_no_cpds)
            self.assertEqual(len(response.context["compounds"]), 0)
            self.assertIn(
                "There are no compounds in this collection, contact MnI to report this if you think it is an error",
                response.content.decode(),
            )
        with self.subTest("Collection has Compounds, no mmps"):
            response = self.client.get(self.url_collection_3_comps)
            self.assertEqual(
                len(response.context["compounds"]), self.collection.compounds().count()
            )
            self.assertEqual(
                response.content.decode().count("No Dotmatics MMPs found with MMPDB"), 3
            )
        with self.subTest("Collection has Compounds, both Compounds have mmps"):
            self.collection_mmp.find_mmps()
            response = self.client.get(self.url_collection_mmp)
            self.assertEqual(
                len(response.context["compounds"]),
                self.collection_mmp.compounds().count(),
            )
            # Check mmps are sorted by similarity
            assayed_comp, mmps = response.context["compounds"][1]
            similarity = [
                mmp.metadata["similarity"][str(assayed_comp.pk)] for mmp in mmps
            ]
            # If mmps overlap, the similarity scores will be sorted incorrectly for those mmps that don't overlap
            self.assertEqual(similarity[1:], sorted(similarity[1:], reverse=True))
            self.assertEqual(
                response.content.decode().count("No Dotmatics MMPs found with MMPDB"), 0
            )


@tag("mmpdb")
class MMPSearchTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.collection = cls.collection_mmp
        cls.co = cls.collection.get_cos_for_analysis(MMP)[0]
        cls.url = reverse(MMP, kwargs={"collection_id": cls.collection.pk})

    def test_unauthenticated(self):
        """
        User must be logged in to use this view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")
        self.assertEqual(self.client.session["next"], self.url)

    def test_context_data(self):
        """
        Test logged in user can see the view
        """
        response = self.client.get(self.url, follow=True)
        group_name = self.collection.get_mmp_group_name(self.co.compound.pk)
        url_w_group = reverse(
            MMP, kwargs={"collection_id": self.collection.pk, "group_name": group_name}
        )
        self.assertRedirects(response, url_w_group)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/mmps/mmps.html")

        # Check task result
        failed, _, result = get_group_results(group_name)
        self.assertFalse(failed)
        self.assertEqual(len(result), 1)
        self.assertEqual(len(result[f"co-{self.co.pk}"]), 2)
