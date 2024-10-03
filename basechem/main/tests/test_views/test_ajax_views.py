import json
from time import sleep
from unittest import mock
from unittest.mock import patch

from django.test import tag
from django.urls import reverse
from django_q.models import Task
from django_q.tasks import async_task

from basechem.common.tests.base import BasechemTestCase, BasechemViewTestMixin
from basechem.main.constants import (
    ALIGN,
    ALOGD,
    CLOGP,
    COMPLETE,
    DOCK,
    DROPPED,
    ERROR,
    ESP,
    IN_PROGRESS,
    MMP,
    MW,
)
from basechem.main.models.compound_models import CompoundOccurrence
from basechem.main.tests.factories import (
    CollectionFactory,
    CompoundFactory,
    CompoundOccurrenceFactory,
    OrmQFactory,
)
from basechem.main.views.ajax_views import AjaxCollectTaskView


class CheckTaskGroupTestCase(BasechemViewTestMixin, BasechemTestCase):
    def setUp(self):
        super().setUp()
        self.url = reverse("ajax_check_task_group")

    def test_get_group(self):
        """
        Tests that the `get` method on the `CheckTaskView` returns the appropriate
        response when the task has not completed yet when `task_type` is "group".
        """
        group_name = "TEST_1"
        async_task(sleep(1), sync=False, group=group_name)
        response = self.client.get(
            self.url, data={"identifier": group_name, "task_type": "group"}
        )
        self.assertEqual(response.content, b"{}")

    def test_get_task(self):
        """
        Tests that the `get` method on the `CheckTaskView` returns the appropriate
        response when the task has not completed yet when `task_type` is "task".
        """
        group_name = "TEST_1"
        task_id = async_task(sleep(1), sync=False, group=group_name)
        response = self.client.get(
            self.url, data={"identifier": task_id, "task_type": "task"}
        )
        self.assertEqual(response.content, b"{}")


class SaveViewerGemsTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.url = reverse("ajax_save_viewer_gems", args=(cls.collection.pk,))
        cls.original_cos = cls.collection.compound_occurrences.all()
        # Run align analysis
        for co in cls.original_cos:
            co.compound.series = cls.series1
            co.compound.save()
        cls.collection.align_analysis("default")
        cls.group_name = cls.collection.get_align_group_name("default")
        cls.co_a = cls.original_cos[0]
        cls.co_b = cls.original_cos[1]
        cls.co_b_conf1 = f"c{cls.co_b.compound.pk}-co{cls.co_b.pk}-1-s-{cls.series1.pk}"
        # Duplicate the conformers for co_a so there are more than 3 to choose from
        task = Task.objects.get(name=f"align_{cls.co_a.pk}_default")
        conf_prefix = f"c{cls.co_a.compound.pk}-co{cls.co_a.pk}"
        cls.co_a_conf1 = f"{conf_prefix}-1-s-{cls.series1.pk}"
        cls.co_a_conf2 = f"{conf_prefix}-2-s-{cls.series1.pk}"
        cls.co_a_conf3 = f"{conf_prefix}-3-s-{cls.series1.pk}"
        cls.co_a_conf4 = f"{conf_prefix}-4-s-{cls.series1.pk}"
        task.result[cls.co_a_conf3] = task.result[cls.co_a_conf1]
        task.result[cls.co_a_conf4] = task.result[cls.co_a_conf2]
        task.save()
        cls.task = task

    def setUp(self):
        super().setUp()
        self.post_data = {"groupName": self.group_name, "viewerConfIds": []}

    def test_failure_too_many_confs(self):
        """
        Test that SaveViewerGems returns errors if a user attempts to save too many conformers
        on a single CompoundOccurrence
        """
        self.post_data["viewerConfIds"] = [
            self.co_a_conf1,
            self.co_a_conf2,
            self.co_a_conf3,
            self.co_a_conf4,
        ]
        response = self.client.generic(
            "POST", self.url, data=json.dumps(self.post_data)
        )
        expected_errors = [
            f"Attempting to save 4 conformers on {self.co_a.compound.name}.",
            "Maximum is 3.",
        ]
        self.assertEqual(response.json()["errors"], expected_errors)
        self.assertEqual(
            list(
                self.collection.compound_occurrences.all().values_list("pk", flat=True)
            ),
            list(self.original_cos.values_list("pk", flat=True)),
        )

    def test_success_blank(self):
        """
        Test that SaveViewerGems successfully does nothing without errors when the viewer is blank
        """
        with self.subTest("Blank save w/ no existing gems"):
            response = self.client.generic(
                "POST", self.url, data=json.dumps(self.post_data)
            )
            expected_errors = []
            self.assertEqual(response.json()["errors"], expected_errors)
            self.assertEqual(
                list(
                    self.collection.compound_occurrences.all().values_list(
                        "pk", flat=True
                    )
                ),
                list(self.original_cos.values_list("pk", flat=True)),
            )
        with self.subTest("Blank save doesn't override existing gems"):
            # Save gem
            gem = CompoundOccurrenceFactory(
                compound=self.co_a.compound,
                parent_co=self.co_a,
                gem_id=self.co_a_conf1,
                molblock=self.task.result[self.co_a_conf1],
            )
            self.collection.compound_occurrences.add(gem)
            self.original_cos = self.collection.compound_occurrences.all()
            # Test blank viewer save
            self.post_data["viewerConfIds"] = []
            response = self.client.generic(
                "POST", self.url, data=json.dumps(self.post_data)
            )
            self.assertEqual(response.json()["errors"], [])
            self.assertEqual(
                list(
                    self.collection.compound_occurrences.all().values_list(
                        "pk", flat=True
                    )
                ),
                list(self.original_cos.values_list("pk", flat=True)),
            )

    def test_success_one_parent_co(self):
        """
        Test that SaveViewerGems can be used to save conformers for a single CompoundOccurrence
        """
        with self.subTest("Save confs for the first time"):
            self.post_data["viewerConfIds"] = [self.co_a_conf1]
            response = self.client.generic(
                "POST", self.url, data=json.dumps(self.post_data)
            )
            self.assertEqual(response.json()["errors"], [])
            self.assertEqual(
                self.collection.compound_occurrences.all().count(),
                self.original_cos.count() + 1,
            )
            new_co = CompoundOccurrence.objects.all().last()
            self.assertEqual(new_co.parent_co, self.co_a)
            self.assertEqual(new_co.gem_id, self.co_a_conf1)

        with self.subTest("Save more confs"):
            self.post_data["viewerConfIds"] = [self.co_a_conf2, self.co_a_conf3]
            response = self.client.generic(
                "POST", self.url, data=json.dumps(self.post_data)
            )
            self.assertEqual(response.json()["errors"], [])
            self.assertEqual(
                self.collection.compound_occurrences.all().count(),
                self.original_cos.count() + 3,
            )
            conf2_co = CompoundOccurrence.objects.get(gem_id=self.co_a_conf2)
            self.assertEqual(conf2_co.parent_co, self.co_a)
            self.assertEqual(conf2_co.gem_id, self.co_a_conf2)
            conf3_co = CompoundOccurrence.objects.get(gem_id=self.co_a_conf3)
            self.assertEqual(conf3_co.parent_co, self.co_a)
            self.assertEqual(conf3_co.gem_id, self.co_a_conf3)

        with self.subTest("Re-save same confs"):
            self.post_data["viewerConfIds"] = [
                self.co_a_conf1,
                self.co_a_conf2,
                self.co_a_conf3,
            ]
            response = self.client.generic(
                "POST", self.url, data=json.dumps(self.post_data)
            )
            self.assertEqual(response.json()["errors"], [])
            self.assertEqual(
                self.collection.compound_occurrences.all().count(),
                self.original_cos.count() + 3,
            )

    def test_success_multiple_parent_cos(self):
        """
        Test that SaveViewerGemsView can be used to save conformers for multiple CompoundOccurrences
        """
        self.post_data["viewerConfIds"] = [self.co_a_conf1, self.co_b_conf1]
        response = self.client.generic(
            "POST", self.url, data=json.dumps(self.post_data)
        )
        self.assertEqual(response.json()["errors"], [])
        self.assertEqual(
            self.collection.compound_occurrences.all().count(),
            self.original_cos.count() + 2,
        )
        new_co_a = CompoundOccurrence.objects.filter(compound=self.co_a.compound).last()
        self.assertEqual(new_co_a.parent_co, self.co_a)
        self.assertEqual(new_co_a.gem_id, self.co_a_conf1)
        new_co_b = CompoundOccurrence.objects.filter(compound=self.co_b.compound).last()
        self.assertEqual(new_co_b.parent_co, self.co_b)
        self.assertEqual(new_co_b.gem_id, self.co_b_conf1)


class UpdateCollectionTestCase(BasechemViewTestMixin, BasechemTestCase):
    def setUp(self):
        super().setUp()
        self.collection = CollectionFactory()
        self.url = reverse("ajax_update_collection", args=(self.collection.pk,))

    def test_update_co_order(self):
        """
        Test that this view can be used to update the co_order in a Collection's metadata
        """
        with self.subTest("First co_order"):
            self.assertNotIn("co_order", self.collection.metadata)
            post_data = {"update_field": "co_order", "co_order": [1, 2, 3]}
            response = self.client.generic("POST", self.url, data=json.dumps(post_data))
            self.assertEqual(response.status_code, 200)
            self.collection.refresh_from_db()
            self.assertEqual(self.collection.metadata["co_order"], [1, 2, 3])

        with self.subTest("Add all new COs to co_order"):
            post_data = {"update_field": "co_order", "co_order": [4]}
            self.client.generic("POST", self.url, data=json.dumps(post_data))
            self.collection.refresh_from_db()
            self.assertEqual(self.collection.metadata["co_order"], [4, 1, 2, 3])

        with self.subTest("Reorder existing COs"):
            post_data = {"update_field": "co_order", "co_order": [1, 4]}
            self.client.generic("POST", self.url, data=json.dumps(post_data))
            self.collection.refresh_from_db()
            self.assertEqual(self.collection.metadata["co_order"], [1, 4, 2, 3])


class AjaxCollectTaskTestCase(BasechemViewTestMixin, BasechemTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.url = reverse("ajax_collect_task")
        cls.collection = cls.collection_one_cpd
        cls.collection.metadata["props_to_show"] = [MW, CLOGP, ALOGD]
        cls.collection.save()
        cls.co = cls.collection.compound_occurrences.all().first()
        cls.c = cls.co.compound
        cls.c.series = cls.series1
        cls.c.save()

        # Run analyses
        cls.collection.propcalc_analysis()
        cls.collection.align_analysis("default")
        cls.collection.dock_analysis("default")
        cls.collection.esp_analysis("default")

    def test_update_task_status(self):
        """
        Test that `update_task_status` updates `self.task_status` appropriately
        """
        view = AjaxCollectTaskView()
        view.task_name = "a_task"
        view.group_name = "a_group"
        view.task = None

        with self.subTest("Task does not exist"):
            view.update_task_status()
            self.assertEqual(view.task_status, DROPPED)

        with self.subTest("Task is on the queue"):
            OrmQFactory()
            with patch("django_q.models.OrmQ.task") as mock_task:
                mock_task.__get__ = mock.Mock(return_value={"name": view.task_name})
                view.update_task_status()
            self.assertEqual(view.task_status, IN_PROGRESS)

        with self.subTest("Task failed"):
            task_id = async_task(
                sleep, 1, task_name=view.task_name, group=view.group_name
            )
            view.task = Task.objects.get(id=task_id)
            view.task.success = False
            view.update_task_status()
            self.assertEqual(view.task_status, ERROR)

        with self.subTest("Task succeeded"):
            view.task.success = True
            view.task.result = "something"
            view.update_task_status()
            self.assertEqual(view.task_status, COMPLETE)

    def test_get_task_DNE(self):
        """
        Test that a get request returns a task result of DROPPED if the requested task does not exist.
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": "a_group",
                "task_name": "this task doesn't exist",
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), {"taskStatus": DROPPED})

    @tag("inductive", "external")
    def test_get_propcalc(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is a propcalc task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_propcalc_group_name(),
                "task_name": self.co.get_propcalc_task_name(self.collection),
                "co_pk": self.co.pk,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()),
            ["taskStatus", "taskResult", "tableRow", "gridItem"],
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        expected_props = [
            "logd_prediction",
            "logd_measured",
            "latest_logd_data_date",
            "mw",
            "clogp",
        ]
        returned_props = list(response.json()["taskResult"].keys())
        self.assertTrue(set(expected_props).issubset(set(returned_props)))
        self.assertIn("</tr>", response.json()["tableRow"])
        self.assertIn("</div>", response.json()["gridItem"])

    def test_get_alignrefs(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is an alignrefs task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_align_group_name("default"),
                "task_name": self.collection.get_refs_task_name(ALIGN),
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()), ["taskStatus", "taskResult", "refOptions"]
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(
            list(response.json()["taskResult"].keys()), ["references", "receptors"]
        )
        self.assertIn("</option>", response.json()["refOptions"])

    def test_get_align(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is an align task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_align_group_name("default"),
                "task_name": self.co.get_align_task_name("default"),
                "co_pk": self.co.pk,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()),
            ["taskStatus", "taskResult", "saveGemsModal", "confOptions"],
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(
            list(response.json()["taskResult"].keys()),
            [
                f"c{self.c.pk}-co{self.co.pk}-1-s-{self.series1.pk}",
                f"c{self.c.pk}-co{self.co.pk}-2-s-{self.series1.pk}",
            ],
        )
        self.assertIn("</option>", response.json()["confOptions"])
        self.assertIn("</form>", response.json()["saveGemsModal"])

    def test_get_dockrefs(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is a dockrefs task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_dock_group_name("default"),
                "task_name": self.collection.get_refs_task_name(DOCK),
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()),
            ["taskStatus", "taskResult", "refOptions", "recOptions"],
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(
            list(response.json()["taskResult"].keys()), ["references", "receptors"]
        )
        self.assertIn("</option>", response.json()["refOptions"])
        self.assertIn("</option>", response.json()["recOptions"])

    def test_get_dock(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is a dock task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_dock_group_name("default"),
                "task_name": self.co.get_dock_task_name("default"),
                "co_pk": self.co.pk,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()),
            ["taskStatus", "taskResult", "saveGemsModal", "confOptions"],
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertTrue(len(response.json()["taskResult"]) > 1)
        for pose_id in response.json()["taskResult"].keys():
            self.assertIn(f"c{self.c.pk}-co{self.co.pk}", pose_id)
        self.assertIn("</option>", response.json()["confOptions"])
        self.assertIn("</form>", response.json()["saveGemsModal"])

    def test_get_esprefs(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is an esprefs task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_esp_group_name("default"),
                "task_name": self.collection.get_refs_task_name(ESP),
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()), ["taskStatus", "taskResult", "refOptions"]
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(
            list(response.json()["taskResult"].keys()), ["references", "receptors"]
        )
        self.assertIn("</option>", response.json()["refOptions"])

    def test_get_esp(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is an esp task
        """
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": self.collection.get_esp_group_name("default"),
                "task_name": self.co.get_esp_task_name("default"),
                "co_pk": self.co.pk,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(list(response.json().keys()), ["taskStatus", "taskResult"])
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(
            list(response.json()["taskResult"].keys()), ["pqr", "related_series", "dx"]
        )

    @patch("basechem.main.models.compound_models.run_mc_torsion_scan")
    def test_get_torsion(self, mock_run_mc_torsion_scan):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is a torsion task
        """
        mock_run_mc_torsion_scan.return_value = (
            "basechem/main/tests/testdata/test_torsion_results.sdf"
        )
        self.collection = self.collection_torsion
        self.co = self.collection.compound_occurrences.all().first()
        self.c = self.co.compound
        dihedral_atoms = "1,2,3,4"
        dihedral_smarts = self.co.convert_atoms_to_smarts(dihedral_atoms)
        group_name = self.collection.get_torsion_group_name(self.co.pk, dihedral_atoms)
        self.collection.torsion_analysis(self.co.pk, dihedral_atoms)

        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection.pk,
                "group_name": group_name,
                "task_name": self.co.get_torsion_task_name(
                    dihedral_smarts, dihedral_atoms
                ),
                "co_pk": self.co.pk,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()), ["taskStatus", "taskResult", "saveGemsModal"]
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(
            list(response.json()["taskResult"].keys()),
            ["torsions", "initial_dihedral", "delta_energy"],
        )
        self.assertIn("</form>", response.json()["saveGemsModal"])

    @tag("mmpdb")
    def test_get_mmp(self):
        """
        Test that a get request returns the expected JsonResponse when the provided task_name is an mmp task
        """
        # Create Compound objects for 2 mmps we expect to find in the mmp_analysis task
        self.mmp1 = CompoundFactory(smiles="CCCCC1CCC(CC)CC1", dn_id="DN9900001")
        self.mmp2 = CompoundFactory(smiles="CCC1CCC(CCl)CC1")
        self.collection_mmp.mmp_analysis()

        self.co = self.collection_mmp.get_cos_for_analysis(MMP)[0]
        response = self.client.get(
            self.url,
            data={
                "collection_pk": self.collection_mmp.pk,
                "group_name": self.collection_mmp.get_mmp_group_name(
                    self.co.compound.pk
                ),
                "task_name": self.co.get_mmp_task_name(),
                "co_pk": self.co.pk,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            list(response.json().keys()),
            ["taskStatus", "taskResult", "dnGrid", "ideaGrid"],
        )
        self.assertEqual(response.json()["taskStatus"], COMPLETE)
        self.assertEqual(type(response.json()["taskResult"]), list)
        self.assertEqual(
            set(response.json()["taskResult"]), {self.mmp1.pk, self.mmp2.pk}
        )
        self.assertIn("</div>", response.json()["dnGrid"])
        self.assertIn(self.mmp1.dn_id, response.json()["dnGrid"])
        self.assertIn("</div>", response.json()["ideaGrid"])
        self.assertIn(self.mmp2.virtual_id, response.json()["ideaGrid"])
