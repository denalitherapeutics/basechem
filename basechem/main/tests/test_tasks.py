import datetime
import os
from io import BytesIO
from time import sleep
from zipfile import ZipFile

import pytz
from django.core import mail
from django.test import tag
from django_q.models import Task
from django_q.tasks import async_task
from rdkit import Chem

from basechem.common.dtx_utils import get_logd_agg_data
from basechem.common.tests.base import BasechemTestCase
from basechem.main.models.collection_models import Collection
from basechem.main.tasks import (
    hide_outdated_results,
    hide_task_group,
    monitor_toklat_scoring,
    update_logd_model_data,
)
from basechem.main.tests.factories import CollectionFactory, CompoundOccurrenceFactory


@tag("local", "dtx")
class TasksTestCase(BasechemTestCase):
    def test_hide_outdated_results(self):
        """
        Test `hide_outdated_results` accurately checks if the required CompoundOccurrences
        are in the currently cached results
        """
        collection = CollectionFactory()
        co = CompoundOccurrenceFactory()
        collection.compound_occurrences.add(co)
        group_name = "torsion_test_group"
        with self.subTest("outdated results"):
            results = {f"co-{co.pk + 1}": {}}
            hidden = hide_outdated_results(collection, group_name, results)
            self.assertTrue(hidden)

        with self.subTest("up-to-date results"):
            results = {f"co-{co.pk + 1}": {}, f"co-{co.pk}": {}}
            hidden = hide_outdated_results(collection, group_name, results)
            self.assertFalse(hidden)

    def test_hide_task_group(self):
        """
        Test `hide_task_group` removes the group name from the relevant tasks
        """
        group_name = "torsion_test_group"
        with self.subTest("Group does not exist"):
            # Assert no error raised
            hide_task_group(group_name)

        with self.subTest("Group does exist"):
            t_id1 = async_task(sleep, 0.01, task_name="task_1", group=group_name)
            t_id2 = async_task(sleep, 0.01, task_name="task_2", group=group_name)

            self.assertEqual(Task.objects.filter(group=group_name).count(), 2)

            hide_task_group(group_name)
            # Check group gone, but tasks still exist
            self.assertEqual(Task.objects.filter(group=group_name).count(), 0)
            self.assertEqual(Task.objects.filter(id__in=[t_id1, t_id2]).count(), 2)

    def test_update_logd_model_data(self):
        """
        Test LogD data is pulled from DTX but not PUT to Inductive (since TEST)
        when there are new mols (variable based on DTX refreshes)
        """
        last_week = datetime.datetime.today() - datetime.timedelta(days=8)
        last_week_fmt = last_week.strftime("%Y%m%d")
        # Check DTX for mols
        mols = get_logd_agg_data(last_week_fmt)
        update_logd_model_data()
        sdf_path = f"/tmp/{last_week_fmt}_logd_dtx.sdf"

        if mols:
            self.assertTrue(os.path.exists(sdf_path))
            for mol in Chem.SDMolSupplier(sdf_path):
                props_dict = mol.GetPropsAsDict()
                self.assertEqual(len(props_dict), 4)
                self.assertTrue("logd_avg" in props_dict.keys())

        else:
            self.assertFalse(os.path.exists(sdf_path))

    def test_monitor_toklat_scoring(self):
        """
        Test that `monitor_toklat_scoring` correctly runs rDock for new Collections and
        summarizes results in an email to the Django admins.
        """
        last_week = datetime.datetime.now(
            pytz.timezone("US/Pacific")
        ) - datetime.timedelta(days=8)
        with self.subTest("No new collections - no email"):
            Collection.objects.all().delete()
            co1 = CompoundOccurrenceFactory()
            col1 = CollectionFactory()
            col1.created_on = last_week
            col1.save()
            col1.compound_occurrences.add(co1)

            monitor_toklat_scoring()
            self.assertEqual(mail.outbox, [])

        with self.subTest("New collection, no errors - email sent"):
            mail.outbox = []
            self.compound_a.series = self.series1
            self.compound_a.save()
            co2 = CompoundOccurrenceFactory(compound=self.compound_a)
            col2 = CollectionFactory()
            col2.compound_occurrences.add(co2)

            monitor_toklat_scoring()
            self.assertEqual(len(mail.outbox), 1)
            email = mail.outbox[0]
            self.assertIn("Toklat vs rDock Scoring Comparison", email.subject)

            # Check file
            self.assertEqual(len(email.attachments), 1)
            file = email.attachments[0]
            num_diff = file[0].split("_")[2]
            self.assertEqual(file[0], f"toklat_compare_{num_diff}_diff.zip")
            zf = ZipFile(BytesIO(file[1]))
            self.assertEqual(
                zf.namelist(), [f"{co2.pk}_{self.series1.dn_id}_rdock_out_scored.sdf"]
            )

            # Check HTML
            html_content = email.alternatives[0][0]
            self.assertIn("Basechem docked 1 CompoundOccurrence ", html_content)
            self.assertIn(
                f"1 CompoundOccurrence had {num_diff} differences in the 3 top poses chosen:",
                html_content,
            )
            self.assertNotIn("Errors encountered:", html_content)
            self.assertNotIn("was skipped", html_content)

        with self.subTest("New collection, has errors - email sent"):
            mail.outbox = []
            co3 = CompoundOccurrenceFactory()
            col2.compound_occurrences.add(co3)

            monitor_toklat_scoring()
            self.assertEqual(len(mail.outbox), 1)
            email = mail.outbox[0]
            self.assertIn("Toklat vs rDock Scoring Comparison", email.subject)

            # Check file
            self.assertEqual(len(email.attachments), 1)
            file = email.attachments[0]
            num_diff = file[0].split("_")[2]
            self.assertEqual(file[0], f"toklat_compare_{num_diff}_diff.zip")
            zf = ZipFile(BytesIO(file[1]))
            self.assertEqual(
                zf.namelist(), [f"{co2.pk}_{self.series1.dn_id}_rdock_out_scored.sdf"]
            )

            # Check HTML
            html_content = email.alternatives[0][0]
            self.assertIn("Basechem docked 2 CompoundOccurrences ", html_content)
            self.assertIn(
                f"1 CompoundOccurrence had {num_diff} differences in the 3 top poses chosen:",
                html_content,
            )
            self.assertIn(f"Error running dock for {co3.compound.name}", html_content)
            self.assertIn(f"1 CompoundOccurrence was skipped", html_content)
