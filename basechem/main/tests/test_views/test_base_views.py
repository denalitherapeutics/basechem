import os
from unittest.mock import patch

from django.core import signing
from django.test import tag
from django.urls import reverse
from django_q.models import Task
from rdkit import Chem

from basechem.common.tests.base import BasechemTestCase, BasechemViewTestMixin
from basechem.common.tests.test_constants import BAD_MOLTEXT, MOLTEXT
from basechem.main.constants import (
    ALIGN,
    ALOGD,
    CLOGP,
    DOCK,
    DTX_MMP,
    ESP,
    MMP,
    MW,
    PROPCALC,
    TORSION,
    TPSA,
)
from basechem.main.forms import CompoundIntakeForm, PropCalcForm
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, CompoundOccurrence
from basechem.main.tasks import get_analysis_results
from basechem.main.tests.factories import CompoundFactory, CompoundOccurrenceFactory
from basechem.main.views.analysis_views import MMPSearchView, PropCalcView
from basechem.main.views.base_views import HomePageView, SubmitCompoundsView


class HomePageTestCase(BasechemViewTestMixin, BasechemTestCase):

    url = reverse("homepage")

    def test_unauthenticated(self):
        """
        User must be logged in order to use the homepage view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")
        self.assertEqual(self.client.session["next"], self.url)
        response = self.client.post(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_authenticated(self):
        """
        Testing the homepage view redirects to the correct page
        """
        response = self.client.get(self.url)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/homepage.html")


class ProjectTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.url = reverse("project", kwargs={"project_code": cls.project1.code})
        cls.cyclohexane = "C1CCCCC1"
        cls.e_cyclohexane = signing.dumps(cls.cyclohexane)
        # Create a Compound that contains cyclohexane but belongs to a different project
        cls.c_wrong_proj = CompoundFactory(smiles="CC1CCCCC1", series=cls.series1)
        # Create a Compound in project1 that contains cyclohexane but is not an exact match
        cls.c_not_exact = CompoundFactory(
            smiles="CC1CCCCC1", project=cls.project1, series=cls.series1
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to use the Project view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")
        self.assertEqual(self.client.session["next"], self.url)
        response = self.client.post(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_authenticated(self):
        """
        Testing the project view redirects to the correct page
        """
        response = self.client.get(self.url)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/project_home.html")
        self.assertIn("compounds", response.context)
        self.assertIn("project", response.context)

    def test_unfiltered(self):
        """
        Tests context data for project is correct when there is no structure filter
        """
        response = self.client.get(self.url)
        expected_cids = [
            c.id for c in Compound.objects.filter(project=self.project1).order_by("id")
        ]
        c_ids = [c.id for c in response.context["compounds"].order_by("id")]
        self.assertEqual(expected_cids, c_ids)
        self.assertEqual(self.project1, response.context["project"])

    def test_filtered_sss(self):
        """
        Tests context data for project is correct when the results are filtered by a substructure search
        """
        with self.subTest("No matches"):
            smiles = "C12CCCCC1CCCN2"
            encrypted_smiles = signing.dumps(smiles)
            url = reverse("project", args=[self.project1.code, "sss", encrypted_smiles])
            response = self.client.get(url)
            self.assertEqual(response.context["smiles"], smiles)
            self.assertEqual(response.context["compounds"].count(), 0)
            self.assertEqual(self.project1, response.context["project"])

        with self.subTest("Some matches"):
            url = reverse(
                "project", args=[self.project1.code, "sss", self.e_cyclohexane]
            )
            response = self.client.get(url)
            expected_cids = sorted(
                [c.id for c in Compound.objects.has_substruct(self.cyclohexane)]
            )
            expected_cids.remove(self.c_wrong_proj.pk)
            c_ids = [c.id for c in response.context["compounds"].order_by("id")]
            self.assertEqual(expected_cids, c_ids)
            self.assertEqual(len(expected_cids), 8)
            self.assertEqual(self.project1, response.context["project"])
            self.assertEqual(response.context["smiles"], self.cyclohexane)

    def test_filtered_exact(self):
        """
        Tests context data for project is correct when the results are filtered by an exact structure search
        """
        with self.subTest("No matches"):
            smiles = "C12CCCCC1CCCN2"
            encrypted_smiles = signing.dumps(smiles)
            url = reverse(
                "project", args=[self.project1.code, "exact", encrypted_smiles]
            )
            response = self.client.get(url)
            self.assertEqual(response.context["compounds"].count(), 0)
            self.assertEqual(self.project1, response.context["project"])
            self.assertEqual(response.context["smiles"], smiles)

        with self.subTest("Some matches"):
            url = reverse(
                "project", args=[self.project1.code, "exact", self.e_cyclohexane]
            )
            response = self.client.get(url)
            expected_cids = sorted(
                [c.id for c in Compound.objects.is_equal(self.cyclohexane)]
            )
            c_ids = [c.id for c in response.context["compounds"].order_by("id")]
            self.assertEqual(expected_cids, c_ids)
            self.assertEqual(len(expected_cids), 1)
            self.assertEqual(self.project1, response.context["project"])
            self.assertEqual(response.context["smiles"], self.cyclohexane)


class SubstructureSearchTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.url = reverse("sss", args=[cls.project1.code])

    def test_unauthenticated(self):
        """
        User must be logged in order to use this view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_fail_invalid_moltext(self):
        """
        Test that the form returns errors when moltext is invalid
        """
        data = {"search_type": "sss", "sketcher": "hello"}
        response = self.client.post(self.url, data)
        self.assertEqual(response.status_code, 200)
        expected_errors = {"sketcher": ["Are you sure this is a valid molecule?"]}
        self.assertEqual(response.json()["errors"], expected_errors)

    def test_success_sss(self):
        """
        Test that the form posts successfully with a substructure structure search
        """
        data = {"search_type": "sss", "sketcher": MOLTEXT}
        response = self.client.post(self.url, data)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["errors"], {})
        redirect_response = self.client.get(response.json()["redirect_url"])
        self.assertEqual(redirect_response.status_code, 200)
        self.assertEqual(redirect_response.context["project"], self.project1)
        self.assertEqual(redirect_response.context["smiles"], "C1CCCCC1")
        self.assertEqual(redirect_response.context["search_type"], "sss")

    def test_success_exact(self):
        """
        Test that the form posts successfully with an exact structure search
        """
        data = {"search_type": "exact", "sketcher": MOLTEXT}
        response = self.client.post(self.url, data)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["errors"], {})
        redirect_response = self.client.get(response.json()["redirect_url"])
        self.assertEqual(redirect_response.status_code, 200)
        self.assertEqual(redirect_response.context["project"], self.project1)
        self.assertEqual(redirect_response.context["smiles"], "C1CCCCC1")
        self.assertEqual(redirect_response.context["search_type"], "exact")


@tag("local", "schrodigner", "dtx")
class SubmitCompoundsTestCase(BasechemViewTestMixin, BasechemTestCase):

    # using homepage as the default next view
    url = reverse("submit", kwargs={"nextview": "homepage"})

    def test_unauthenticated(self):
        """
        User must be logged in order to use the submitcompounds view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")
        self.assertEqual(self.client.session["next"], self.url)
        response = self.client.post(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_submit_compounds_view(self):
        """
        Test logged in user can see the submitcompounds view
        """
        response = self.client.get(self.url, follow=True)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "main/submitcompounds.html")

    def test_generic_upload(self):
        """
        Tests upload to the generic form
        """
        with self.subTest("Invalid Upload"):
            data = {"project": ""}
            response = self.client.post(self.url, data=data, follow=True)
            self.assertEqual(len(response.context["form"].fields), 6)
            self.assertIsNotNone(response.context["form"].errors)
            self.assertEqual(
                ["This field is required."], response.context["form"].errors["project"]
            )
            # Check field error exists but has no message
            self.assertEqual([""], response.context["form"].errors["upload_file"])
            self.assertEqual([""], response.context["form"].errors["moltext"])

            data = {"project": self.test_project}
            response = self.client.post(self.url, data=data, follow=True)
            self.assertEqual(len(response.context["form"].fields), 6)
            self.assertIsNotNone(response.context["form"].errors)
            # Check field error exists but has no message
            self.assertEqual([""], response.context["form"].errors["upload_file"])
            self.assertEqual([""], response.context["form"].errors["moltext"])

            response = self.upload_bad_sdf(self.url)
            self.assertFalse(response.context["form"].is_valid())
            self.assertIsNotNone(response.context["form"].errors)
            self.assertEqual(
                ["Your uploaded data did not return any valid compounds"],
                response.context["form"].errors["__all__"],
            )

            response = self.upload_bad_mol(self.url)
            self.assertFalse(response.context["form"].is_valid())
            self.assertIsNotNone(response.context["form"].errors)
            self.assertEqual(
                ["Are you sure this is valid moltext?"],
                response.context["form"].errors["moltext"],
            )
        with self.subTest("Valid Upload redirects to homepage"):
            response = self.upload_sdf(self.url)
            self.assertEqual(response.status_code, 200)
            self.assertIsInstance(response.context["view"], HomePageView)

            response = self.upload_mol(self.url)
            self.assertEqual(response.status_code, 200)
            self.assertIsInstance(response.context["view"], HomePageView)

    def test_use_existing_collection(self):
        """
        Test submit compounds with an existing collection
        """
        num_old_collections = len(Collection.objects.all())
        collection = Collection.objects.last()
        self.client.login(username=collection.owner.username, password="testpassword")
        url = reverse("submit", kwargs={"nextview": PROPCALC})
        data = {"collection": collection.pk, "physiochemical": [MW, TPSA, CLOGP]}
        response = self.client.post(url, data=data, follow=True)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context.get("collection"), collection)
        self.assertIn("property_results", response.context)
        new_collections = Collection.objects.all()
        self.assertEqual(num_old_collections, len(new_collections))

    def test_propcalc_invalid_upload(self):
        """
        Tests upload to the propcalc form, expects errors
        """
        url = reverse("submit", kwargs={"nextview": PROPCALC})
        data = {
            "project": "",
        }
        response = self.client.post(url, data=data, follow=True)
        self.assertEqual(len(response.context["form"].fields), 8)
        self.assertIsNotNone(response.context["form"].errors)
        self.assertEqual(
            ["This field is required."], response.context["form"].errors["project"]
        )
        # Check field error exists but has no message
        self.assertEqual([""], response.context["form"].errors["upload_file"])
        self.assertEqual([""], response.context["form"].errors["moltext"])

        data = {"project": self.test_project}
        response = self.client.post(url, data=data, follow=True)
        self.assertEqual(len(response.context["form"].fields), 8)
        self.assertIsNotNone(response.context["form"].errors)
        # Check field error exists but has no message
        self.assertEqual([""], response.context["form"].errors["upload_file"])
        self.assertEqual([""], response.context["form"].errors["moltext"])

        response = self.upload_bad_sdf(self.url)
        self.assertFalse(response.context["form"].is_valid())
        self.assertIsNotNone(response.context["form"].errors)
        self.assertEqual(
            ["Your uploaded data did not return any valid compounds"],
            response.context["form"].errors["__all__"],
        )

        response = self.upload_bad_mol(self.url)
        self.assertFalse(response.context["form"].is_valid())
        self.assertIsNotNone(response.context["form"].errors)
        self.assertEqual(
            ["Your uploaded data did not return any valid compounds"],
            response.context["form"].errors["__all__"],
        )

    @tag("local", "inductive", "external")
    def test_propcalc_valid_upload_inductive(self):
        """
        Tests upload to the propcalc form with Inductive Props requested
        """
        url = reverse("submit", kwargs={"nextview": PROPCALC})
        data = {"physiochemical": [MW, TPSA, CLOGP, ALOGD]}
        response = self.upload_sdf(url, data)

        self.assertEqual(response.status_code, 200)
        self.assertIsInstance(response.context["view"], PropCalcView)
        self.assertIn("collection", response.context)
        self.assertIn("property_results", response.context)

    def test_propcalc_valid_upload(self):
        """
        Tests upload to the propcalc form with only RDKit props
        """
        url = reverse("submit", kwargs={"nextview": PROPCALC})
        data = {"physiochemical": [MW, TPSA, CLOGP]}
        response = self.upload_mol(url, data)

        self.assertEqual(response.status_code, 200)
        self.assertIsInstance(response.context["view"], PropCalcView)
        self.assertIn("collection", response.context)
        self.assertIn("property_results", response.context)

    def test_mmp_valid_upload(self):
        """
        Test a valid upload to the MMPSubmitForm
        """
        response = self.upload_mmp_mol(variable={"atoms": [6, 7]})
        self.assertEqual(response.status_code, 200)
        self.assertIsInstance(response.context["view"], MMPSearchView)
        self.assertIn("collection", response.context)
        collection = response.context["collection"]
        compound_pk = collection.compounds().first().pk
        mmp_metadata_dict = {
            str(compound_pk): {
                "constant_smiles": "*C1CCC(CC)CC1",
                "variable_smiles": "*CF",
            }
        }
        self.assertEqual(collection.metadata["mmp_analysis"], mmp_metadata_dict)

    def test_mmp_invalid_upload(self):
        """
        Test an invalid upload to the MMPSubmitForm
        """
        url = reverse("submit", kwargs={"nextview": MMP})
        response = self.upload_mmp_mol(variable={"atoms": []})
        self.assertIsNotNone(response.context["form"].errors)
        expected_errors = ["You must highlight a variable region"]
        self.assertEqual(response.context["form"].errors["__all__"], expected_errors)

    def test_get(self):
        """
        Tests different nextviews will have different forms attached to the view
        """
        with self.subTest("Test compound intake form is attached"):
            url = reverse("submit", kwargs={"nextview": "homepage"})
            response = self.client.get(url)
            self.assertIsInstance(response.context["view"], SubmitCompoundsView)
            self.assertIsInstance(response.context["form"], CompoundIntakeForm)

        with self.subTest("Test prop calc form is attached"):
            url = reverse("submit", kwargs={"nextview": PROPCALC})
            response = self.client.get(url)
            self.assertIsInstance(response.context["view"], SubmitCompoundsView)
            self.assertIsInstance(response.context["form"], PropCalcForm)

    def test_get_context_data(self):
        """
        Tests different nextviews will have different context data
        """
        with self.subTest("Test homepage context data"):
            url = reverse("submit", kwargs={"nextview": "homepage"})
            response = self.client.get(url)
            self.assertIn("form", response.context)
            self.assertNotIn("nextview", response.context)
            self.assertNotIn("title", response.context)
            self.assertNotIn("directions", response.context)

        with self.subTest("Test propcalc context data"):
            url = reverse("submit", kwargs={"nextview": PROPCALC})
            response = self.client.get(url)
            self.assertIn("form", response.context)
            self.assertIn("nextview", response.context)
            self.assertIn("title", response.context)
            self.assertIn("directions", response.context)

    def test_form_valid(self):
        """
        Tests valid forms redirect correctly
        """
        with self.subTest("Test form valid directs to homepage"):
            # In practice homepage should never be a nextview but for test purposes we
            # are checking it redirects correctly and that as an else case, it behaves
            # as expected, with no collections being created.
            old_collections = Collection.objects.all()
            url = reverse("submit", kwargs={"nextview": "homepage"})
            response = self.upload_sdf(url)
            self.assertEqual(response.status_code, 200)
            self.assertIsInstance(response.context["view"], HomePageView)
            new_collections = Collection.objects.all()
            self.assertEqual(len(old_collections), len(new_collections))

        with self.subTest("Test form valid directs to propcalc"):
            num_old_collections = len(Collection.objects.all())
            url = reverse("submit", kwargs={"nextview": PROPCALC})
            response = self.upload_sdf(url)
            self.assertEqual(response.status_code, 200)
            self.assertIsInstance(response.context["view"], PropCalcView)
            self.assertIn("collection", response.context)
            self.assertIn("property_results", response.context)
            new_collections = Collection.objects.all()
            self.assertEqual(num_old_collections + 1, len(new_collections))


class DownloadResultFileTestCase(BasechemViewTestMixin, BasechemTestCase):
    def test_get(self):
        """
        Tests get returns response object with filename attachment
        """
        with self.subTest("propcalc"):
            url = reverse(
                "downloadresultfile", args=[PROPCALC, self.collection_one_cpd.id]
            )
            response = self.client.get(url)
            filepath, filename = self.collection_one_cpd.generate_file(
                PROPCALC, test=True
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        with self.subTest("ligand align download all"):
            # Run align analysis (download relies on completed task)
            align_url = reverse(ALIGN, kwargs={"collection_id": self.collection.id})
            c_ids = [c.id for c in self.collection.compounds()]
            form_data = {
                "csrfmiddlewaretoken": "token",
                "compounds": c_ids,
                "reference": f"s-{self.series1.id}",
            }
            json_response = self.client.post(align_url, data=form_data, follow=True)
            response = self.client.get(json_response.json()["redirect_url"])
            group_name = response.context.get("group_name")
            # Test download align file
            url = reverse(
                "downloadresultfile", args=[ALIGN, self.collection.id, group_name]
            )
            response = self.client.get(url)
            filepath, filename = self.collection.generate_file(
                ALIGN, test=True, group_name=group_name
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        with self.subTest("ligand align download selected"):
            # Use same analysis/group_name as previous test
            url = reverse(
                "downloadresultfile",
                args=[ALIGN, self.collection.id, group_name, f"s-{self.series1.id}"],
            )
            response = self.client.get(url)
            filepath, filename = self.collection.generate_file(
                ALIGN,
                test=True,
                group_name=group_name,
                selected_ids=f"s-{self.series1.id}",
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        with self.subTest("dtx mmp download all"):
            # Add measured data to a Compound to confirm that it is included in the downloaded file
            comp = self.collection_mmp.compounds().order_by("dn_id")[1]
            comp.measured_data["assay_results"] = {"assay 1": {"analysis 1": 1}}
            comp.save()

            url = reverse("downloadresultfile", args=[DTX_MMP, self.collection_mmp.id])
            response = self.client.get(url)
            filepath, filename = self.collection_mmp.generate_file(DTX_MMP, test=True)
            mols = [mol for mol in Chem.SDMolSupplier(filepath)]
            # Check just compounds, no mmps are downloaded
            self.assertEqual(len(mols), 3)
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])
            # Check assay data is included
            self.assertEqual(mols[1].GetProp("assay 1: analysis 1"), "1")

            # Assign mmps and regenerate file
            self.collection_mmp.find_mmps()
            mmp = comp.mmps.exclude(dn_id="").order_by("dn_id")[0]
            mmp.measured_data["assay_results"] = {"assay 1": {"analysis 1": 2}}
            mmp.save()
            response = self.client.get(url)
            filepath, filename = self.collection_mmp.generate_file(DTX_MMP, test=True)

            mols = [mol for mol in Chem.SDMolSupplier(filepath)]
            self.assertEqual(len(mols), 12)
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])
            # TODO: Check assay data is included
            # self.assertEqual(mols[2].GetProp("assay 1: analysis 1"), "2")

    @tag("local", "esp")
    def test_get_local(self):
        """
        [For analyses that can only be tested locally because they rely on the efs mount]
        Tests get returns response object with filename attachment
        """
        with self.subTest("dock download all"):
            # Run dock analysis (download relies on completed task)
            dock_url = reverse(DOCK, kwargs={"collection_id": self.collection.id})
            c_ids = [c.id for c in self.collection.compounds()]
            form_data = {
                "csrfmiddlewaretoken": "token",
                "compounds": c_ids,
                "reference": f"s-{self.series1.id}",
            }
            json_response = self.client.post(dock_url, data=form_data, follow=True)
            response = self.client.get(json_response.json()["redirect_url"])
            group_name = response.context.get("group_name")
            # Test download dock file
            url = reverse(
                "downloadresultfile", args=[DOCK, self.collection.id, group_name]
            )
            response = self.client.get(url)
            filepath, filename = self.collection.generate_file(
                DOCK, test=True, group_name=group_name
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        with self.subTest("dock download selected"):
            # Use same analysis/group_name as previous test
            url = reverse(
                "downloadresultfile",
                args=[DOCK, self.collection.id, group_name, f"s-{self.series1.id}"],
            )
            response = self.client.get(url)
            filepath, filename = self.collection.generate_file(
                DOCK,
                test=True,
                group_name=group_name,
                selected_ids=f"s-{self.series1.id}",
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        with self.subTest("esp download all"):
            # Run esp predict (download relies on completed task)
            esp_url = reverse(ESP, kwargs={"collection_id": self.collection.id})
            c_ids = [c.id for c in self.collection.compounds()]
            form_data = {
                "csrfmiddlewaretoken": "token",
                "compounds": c_ids,
                "reference": f"s-{self.series1.id}",
            }
            json_response = self.client.post(esp_url, data=form_data, follow=True)
            response = self.client.get(json_response.json()["redirect_url"])
            group_name = response.context.get("group_name")
            # Test download dock file
            url = reverse(
                "downloadresultfile", args=[ESP, self.collection.id, group_name]
            )
            response = self.client.get(url)
            filepath, filename = self.collection.generate_file(
                ESP, test=True, group_name=group_name
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        with self.subTest("esp download selected"):
            # Use same analysis/group_name as previous test
            url = reverse(
                "downloadresultfile",
                args=[ESP, self.collection.id, group_name, f"s-{self.series1.id}"],
            )
            response = self.client.get(url)
            filepath, filename = self.collection.generate_file(
                ESP,
                test=True,
                group_name=group_name,
                selected_ids=f"s-{self.series1.id}",
            )
            self.assertTrue(os.path.exists(filepath))
            self.assertIn(filename, response["Content-Disposition"])

        # Patch the function that runs the torsion scan so we get the expected output file
        # without having to run a job on the Slurm cluster.
        with patch(
            "basechem.main.models.compound_models.run_mc_torsion_scan"
        ) as mock_run_mc_torsion_scan:
            mock_run_mc_torsion_scan.return_value = (
                "basechem/main/tests/testdata/test_torsion_results.sdf"
            )

            with self.subTest("torsion download all"):
                # Run torsion analysis (download relies on completed task)
                torsion_url = reverse(
                    TORSION, kwargs={"collection_id": self.collection_torsion.id}
                )
                co = self.collection_torsion.compound_occurrences.first()
                c = co.compound
                form_data = {
                    "csrfmiddlewaretoken": "token",
                    "pioneer": f"co-{co.pk}",
                    "selected_atoms": "0,1,2,3",
                }
                json_response = self.client.post(torsion_url, data=form_data)
                response = self.client.get(json_response.json()["redirect_url"])
                group_name = response.context.get("group_name")

                # Test download torsion file
                url = reverse(
                    "downloadresultfile",
                    args=[TORSION, self.collection_torsion.id, group_name],
                )
                response = self.client.get(url)
                filepath, filename = self.collection_torsion.generate_file(
                    TORSION, test=True, group_name=group_name
                )
                mols = [mol for mol in Chem.SDMolSupplier(filepath)]
                self.assertEqual(len(mols), 37)

                self.assertTrue(os.path.exists(filepath))
                self.assertIn(filename, response["Content-Disposition"])

            with self.subTest("torsion download selected viewer"):
                # Use same analysis/group_name as previous test
                url = reverse(
                    "downloadresultfile",
                    args=[
                        TORSION,
                        self.collection_torsion.id,
                        group_name,
                        f"c{c.pk}-co{co.pk}-120_c{c.pk}-co{co.pk}-180",
                    ],
                )
                response = self.client.get(url)
                filepath, filename = self.collection_torsion.generate_file(
                    TORSION,
                    test=True,
                    group_name=group_name,
                    selected_ids=f"c{c.pk}-co{co.pk}-120_c{c.pk}-co{co.pk}-180",
                )
                mols = [mol for mol in Chem.SDMolSupplier(filepath)]
                self.assertEqual(len(mols), 2)
                self.assertTrue(os.path.exists(filepath))
                self.assertIn(filename, response["Content-Disposition"])


class HikeTestCase(BasechemViewTestMixin, BasechemTestCase):
    def setUp(self):
        super().setUp()
        self.url = reverse("hike", kwargs={"collection_id": self.collection.id})
        self.client.login(
            username=self.collection.owner.username, password="testpassword"
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to use the hike view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")
        self.assertEqual(self.client.session["next"], self.url)
        response = self.client.post(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_hike_to_propcalc(self):
        form_data = {"analysis": PROPCALC}
        response = self.client.post(self.url, form_data)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(f"/{PROPCALC}/submit/{self.collection.id}", response.url)

    def test_hike_to_align(self):
        form_data = {"analysis": ALIGN}
        response = self.client.post(self.url, form_data)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(f"/{ALIGN}/{self.collection.id}", response.url)

    def test_hike_to_dock(self):
        # You can only hike when the structure is available
        form_data = {"analysis": DOCK}
        response = self.client.post(self.url, form_data)
        self.assertEqual(response.status_code, 302)
        self.assertEqual(f"/{DOCK}/{self.collection.id}", response.url)


@tag("local")
class SaveGemsTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.collection = cls.collection_one_cpd
        cls.co = cls.collection.compound_occurrences.all().first()
        cls.c = cls.co.compound
        # Run align analysis to have conformers to choose from
        co_ids = []
        for co in cls.collection.compound_occurrences.all():
            co_ids.append(co.id)
            co.compound.series = cls.series1
            co.compound.save()
        ref_string = f"s-{cls.series1.pk}"
        analysis_kwargs = {"ref_string": ref_string}
        cls.collection.align_analysis(ref_string)
        cls.group_name = cls.collection.get_align_group_name(ref_string)
        cls.task_name = cls.co.get_align_task_name(ref_string)
        # Duplicate the conformers so that there are >3 to choose from
        task = Task.objects.get(name=cls.task_name)
        conf_prefix = f"c{cls.c.pk}-co{cls.co.pk}"
        cls.conf_id1 = f"{conf_prefix}-1-s-{cls.series1.pk}"
        cls.conf_id2 = f"{conf_prefix}-2-s-{cls.series1.pk}"
        cls.conf_id3 = f"{conf_prefix}-3-s-{cls.series1.pk}"
        cls.conf_id4 = f"{conf_prefix}-4-s-{cls.series1.pk}"
        task.result[cls.conf_id3] = task.result[cls.conf_id1]
        task.result[cls.conf_id4] = task.result[cls.conf_id2]
        task.save()

        _, cls.align_results = get_analysis_results(
            cls.collection, ALIGN, cls.group_name, analysis_kwargs
        )
        cls.url = reverse(
            "save_gems",
            kwargs={
                "collection_id": cls.collection.pk,
                "compound_occurrence_id": cls.co.pk,
                "group_name": cls.group_name,
                "task_name": cls.task_name,
            },
        )

    def setUp(self):
        super().setUp()
        self.client.login(
            username=self.collection.owner.username, password="testpassword"
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to use the hike view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")
        self.assertEqual(self.client.session["next"], self.url)
        response = self.client.post(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_failure_too_many_new(self):
        """
        Test saving more than 3 new conformers fails
        """
        expected_last_co = CompoundOccurrence.objects.all().order_by("pk").last()
        data = {
            "gems": [self.conf_id1, self.conf_id2, self.conf_id3, self.conf_id4],
            "other_gems": [],
        }
        response = self.client.post(self.url, data=data, follow=True)
        self.assertEqual(response.status_code, 200)
        content = response.json()
        # Check errors are as expected
        expected_errors = {"__all__": ["You cannot save more than 3 gems at a time."]}
        self.assertEqual(content["errors"], expected_errors)

        # Check no new COs were added to the db
        last_co = CompoundOccurrence.objects.all().order_by("pk").last()
        self.assertEqual(expected_last_co, last_co)
        # Check the collection's COs are unchanged
        current_co_pks = [co.pk for co in self.collection.compound_occurrences.all()]
        self.assertEqual(current_co_pks, [self.co.pk])

    def test_failure_too_many_old_and_new(self):
        """
        Test saving more than 3 new and old gems fails
        """
        dock_gem_id = f"c{self.c.pk}-co{self.co.pk}-5"
        # Add existing DOCK CompoundOccurrence
        dock_co = CompoundOccurrenceFactory(
            parent_co=self.co,
            compound=self.co.compound,
            saved_from=DOCK,
            gem_id=dock_gem_id,
            molblock=self.align_results["compounds"][f"co-{self.co.pk}"][self.conf_id1][
                "moltext"
            ],
        )
        self.collection.compound_occurrences.add(dock_co)
        expected_last_co = dock_co
        data = {
            "gems": [self.conf_id1, self.conf_id2, self.conf_id3],
            "other_gems": [dock_gem_id],
        }
        response = self.client.post(self.url, data=data, follow=True)
        self.assertEqual(response.status_code, 200)
        content = response.json()
        # Check errors are as expected
        expected_errors = {"__all__": ["You cannot save more than 3 gems at a time."]}
        self.assertEqual(content["errors"], expected_errors)

        # Check no new COs were added to the db
        last_co = CompoundOccurrence.objects.all().order_by("pk").last()
        self.assertEqual(expected_last_co, last_co)

        # Check the collection's COs are unchanged
        current_co_pks = [co.pk for co in self.collection.compound_occurrences.all()]
        self.assertEqual(sorted(current_co_pks), sorted([self.co.pk, dock_co.pk]))

    def test_success_saved(self):
        """
        Test successfully saving gems
        """
        with self.subTest("None"):
            expected_last_co = CompoundOccurrence.objects.all().order_by("pk").last()
            data = {"gems": [], "other_gems": []}
            response = self.client.post(self.url, data=data, follow=True)
            self.assertEqual(response.status_code, 200)
            content = response.json()
            # Check no errors
            self.assertEqual(content["errors"], {})

            # Check no new COs were added to the db
            last_co = CompoundOccurrence.objects.all().order_by("pk").last()
            self.assertEqual(expected_last_co, last_co)
            # Check the collection's COs are unchanged
            current_co_pks = [
                co.pk for co in self.collection.compound_occurrences.all()
            ]
            self.assertEqual(current_co_pks, [self.co.pk])

        with self.subTest("Save 3"):
            prev_num_cos = CompoundOccurrence.objects.all().count()
            data = {
                "gems": [self.conf_id1, self.conf_id2, self.conf_id3],
                "other_gems": [],
            }
            response = self.client.post(self.url, data=data, follow=True)
            self.assertEqual(response.status_code, 200)
            content = response.json()
            # Check no errors
            self.assertEqual(content["errors"], {})

            # Check new COs have been created and saved to the collection
            cos = list(self.collection.compound_occurrences.all().order_by("pk"))
            self.assertEqual(cos[0], self.co)
            for i in range(1, 4):
                self.assertEqual(
                    cos[i].gem_id,
                    f"c{self.c.pk}-co{self.co.pk}-{i}-s-{self.series1.pk}",
                )
                self.assertEqual(cos[i].saved_from, ALIGN)
                self.assertEqual(cos[i].parent_co, self.co)
                self.assertEqual(cos[i].compound, self.co.compound)

            # Check only 3 COs made:
            num_cos = CompoundOccurrence.objects.all().count()
            self.assertEqual(prev_num_cos + 3, num_cos)

        with self.subTest("Change saved"):
            original_last_co = CompoundOccurrence.objects.all().order_by("pk").last()
            data = {
                "gems": [self.conf_id1, self.conf_id4],
                "other_gems": [],
            }
            response = self.client.post(self.url, data=data, follow=True)
            self.assertEqual(response.status_code, 200)
            content = response.json()
            # Check no errors
            self.assertEqual(content["errors"], {})

            # Check 1 new CO was created and saved to the collection
            collection_cos = self.collection.compound_occurrences.all()
            new_co = CompoundOccurrence.objects.all().last()
            self.assertEqual(original_last_co.pk + 1, new_co.pk)
            self.assertEqual(collection_cos.count(), 3)
            self.assertIn(new_co, collection_cos)

            # Check conf_id1 remains saved to the collection
            conf_1 = CompoundOccurrence.objects.get(gem_id=self.conf_id1)
            self.assertIn(conf_1, collection_cos)

            # Check 2 unsaved COs are not in the collection anymore
            conf_2 = CompoundOccurrence.objects.get(gem_id=self.conf_id2)
            self.assertNotIn(conf_2, collection_cos)
            conf_3 = CompoundOccurrence.objects.get(gem_id=self.conf_id3)
            self.assertNotIn(conf_3, collection_cos)


@tag("local", "dtx")
class AddCompoundTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.collection = cls.collection_2d
        new_mol_smiles = "C1=CCC=C1"
        cls.moltext = Chem.MolToMolBlock(Chem.MolFromSmiles("C1=CCC=C1"))
        cls.data = {"sketcher": cls.moltext}
        cls.url = reverse(
            "add_comp",
            kwargs={"current_view": ALIGN, "collection_id": cls.collection.pk},
        )

    def test_unauthenticated(self):
        """
        User must be logged in order to use the hike view
        """
        self.client.logout()
        response = self.client.post(self.url, data=self.data)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_failure_comp_already_exists(self):
        """
        Test that the form comes back with errors if the drawn compound is already in the collection
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        data = {"sketcher": self.collection.compound_occurrences.first().moltext()}
        response = self.client.post(self.url, data=data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 1)
        self.assertEqual(
            content["errors"],
            {"__all__": ["This molecule is already in the collection"]},
        )
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        self.assertEqual(Task.objects.all().count(), num_tasks)

    def test_failure_invalid_moltext(self):
        """
        Test that the form comes back with errors if the drawn compound is invalid moltext
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        data = {"sketcher": BAD_MOLTEXT}
        url = reverse(
            "add_comp",
            kwargs={"current_view": ALIGN, "collection_id": self.collection.pk},
        )
        response = self.client.post(self.url, data=data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 1)
        self.assertEqual(
            content["errors"], {"sketcher": ["Are you sure this is a valid molecule?"]}
        )
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        self.assertEqual(Task.objects.all().count(), num_tasks)

    def test_success_propcalc(self):
        """
        Test successful form posts when adding a compound in propcalc
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        group_name = self.collection.get_propcalc_group_name()
        url = reverse(
            "add_comp",
            kwargs={
                "current_view": PROPCALC,
                "collection_id": self.collection.pk,
                "group_name": group_name,
            },
        )
        response = self.client.post(url, data=self.data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 3)
        self.assertEqual(content["errors"], {})
        self.assertEqual(self.collection.compound_occurrences.count(), 2)
        self.assertEqual(Task.objects.all().count(), num_tasks + 1)

    def test_success_align_before_submit(self):
        """
        Test successful form posts when adding a compound in align before starting the analysis
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        url = reverse(
            "add_comp",
            kwargs={"current_view": ALIGN, "collection_id": self.collection.pk},
        )
        response = self.client.post(url, data=self.data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 3)
        self.assertEqual(content["errors"], {})
        self.assertEqual(self.collection.compound_occurrences.count(), 2)
        self.assertEqual(Task.objects.all().count(), num_tasks)

    def test_success_align_after_submit(self):
        """
        Test successful form posts when adding a compound in align after starting the analysis
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        group_name = self.collection.get_align_group_name("default")
        url = reverse(
            "add_comp",
            kwargs={
                "current_view": ALIGN,
                "collection_id": self.collection.pk,
                "group_name": group_name,
            },
        )
        response = self.client.post(url, data=self.data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 3)
        self.assertEqual(content["errors"], {})
        self.assertEqual(self.collection.compound_occurrences.count(), 2)
        self.assertEqual(Task.objects.all().count(), num_tasks + 1)

    def test_success_dock(self):
        """
        Test a successful form post when adding a compound in dock
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        group_name = self.collection.get_dock_group_name("default")
        url = reverse(
            "add_comp",
            kwargs={
                "current_view": DOCK,
                "collection_id": self.collection.pk,
                "group_name": group_name,
            },
        )
        response = self.client.post(url, data=self.data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 3)
        self.assertEqual(content["errors"], {})
        self.assertEqual(self.collection.compound_occurrences.count(), 2)
        self.assertEqual(Task.objects.all().count(), num_tasks + 1)

    def test_success_esp(self):
        """
        Test a successful form post when adding a compound in esp
        """
        num_tasks = Task.objects.all().count()
        self.assertEqual(self.collection.compound_occurrences.count(), 1)
        group_name = self.collection.get_esp_group_name("default")
        url = reverse(
            "add_comp",
            kwargs={
                "current_view": ESP,
                "collection_id": self.collection.pk,
                "group_name": group_name,
            },
        )
        response = self.client.post(url, data=self.data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 3)
        self.assertEqual(content["errors"], {})
        self.assertEqual(self.collection.compound_occurrences.count(), 2)
        self.assertEqual(Task.objects.all().count(), num_tasks + 1)

    @patch("basechem.main.models.compound_models.run_mc_torsion_scan")
    def test_success_torsion(self, mock_run_mc_torsion_scan):
        """
        Test a successful form post when adding a compound in torsion
        """
        # Create another CO to be the pioneer so we can add the cached torsion compound (CCCO)
        # to the collection
        collection = self.collection_torsion
        self.torsion_co = collection.compound_occurrences.first()
        collection.compound_occurrences.clear()
        compound = Compound.objects.create(smiles="CCCCO", project=collection.project)
        self.pioneer = CompoundOccurrence.objects.create(
            compound=compound, owner=collection.owner
        )
        collection.compound_occurrences.add(self.pioneer)
        mock_run_mc_torsion_scan.return_value = (
            "basechem/main/tests/testdata/test_torsion_results.sdf"
        )

        data = {"sketcher": self.torsion_co.moltext()}
        group_name = collection.get_torsion_group_name(self.pioneer.pk, "1,2,3,4")
        num_tasks = Task.objects.all().count()
        self.assertEqual(collection.compound_occurrences.count(), 1)
        url = reverse(
            "add_comp",
            kwargs={
                "current_view": TORSION,
                "collection_id": collection.pk,
                "group_name": group_name,
            },
        )
        response = self.client.post(url, data=data, follow=True)
        content = response.json()

        self.assertEqual(len(content), 3)
        self.assertEqual(content["errors"], {})
        self.assertEqual(collection.compound_occurrences.count(), 2)
        self.assertEqual(Task.objects.all().count(), num_tasks + 1)
