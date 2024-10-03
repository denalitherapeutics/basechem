import datetime
import os

import pandas as pd
from django.conf import settings
from django.core import mail
from django.test import tag

from basechem.common.dtx_utils import ASSAY_DATA_COLUMNS
from basechem.common.tests.base import BasechemTestCase
from basechem.common.tests.test_constants import MOLTEXT
from basechem.main.assay_emailer_utils import (
    _all_assays_in,
    _any_assays_in,
    _assays_needing_pings,
    _clean_assay_names,
    _create_assay_html_table,
    _create_curve_file,
    _create_mol_image_file,
    _delete_files,
    _impatient_assays,
    depreciate_old_assays,
    process_new_data,
    send_new_data_assay_email,
    send_new_data_ping,
    update_assay_results_from_new_data_df,
    update_series_data,
)
from basechem.main.constants import MAX_HOURS_TO_WAIT
from basechem.main.models.collection_models import Collection
from basechem.main.tests.factories import BasechemUserFactory, CompoundFactory

ASSAY1 = "assay one long name"
ASSAY1_SHORT = "assay one long"
ASSAY2 = "assay two"
DN_ID = "DN0017079"
SAMPLE_IC50 = "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\n11v9tepLSkqqvP/25bvtv719LUtE/yN9v/b+OS3DXm9/fWX7H2k7jWX7c6rXrK6o\naJ+8+uoB4+yc2v/yyzLZsKHEcX2th/3z0+0olt9LtP6tWLHCCKBOnb7/SNuX03xK\nCWB2dracOHHCOMBevXrF7QC9PP7Jzb837fu/fa+XoeJk7XY+/LBI69bO7cf70HbV\n2U2eXDWfERmgB4wYMcJ8MOvWrTN5YCIEkAeP6n7YxUl3OytzJkUsWd/vHiIhxwGS\nASaYuXPnSlZWlhkJLioqSgoHGC2D0QNDww0iJMMxXNV9HOBf/xrc9/bt0zsDVHGq\niu/Ay22IDNAHVNcosNf/3qHl8f57R/vnbtYsuG2dj7f/9rpEjwLH8/nH8x1wNRgc\nYFoIIBkaGWAyZWh+hwwQB5hUDjAR/U81BxitnzhAHGDKCmB1HwcYSx6SzscBxvPd\nxfMdcD1AMkAcoA/+vavbAbrtvxcOCgeIA0QAASBpIAPEASaVA0xE/3GAOEAcYIoI\nINcDdFfP9QDJAJ2WSQaIA8QB4gBxgDhAfwggAHgHGSAOEAeIA8QB4gD9IYBkgO7q\nyQDJAJ2WSQaIA8QB4gBxgDhAfwggAHgHGWCCePXVV82HwYMHj/R56O8eAQQAQAAB\nABBAAAAEEAAAAQQABJCPAAAQQAAABBAAUvLHzt2bEEAAP/H+++/LpZdeKk2bNpW3\n337bUcw+/vhjadKkidx8883y4osvSoMGDeTKK680rwn3aNWqldSpU0cyMjJk48aN\nCCAAJI5IghSJhg0byooVK2Tfvn0yefJkeeSRR+TkyZNB79O29VSyL774Qs477zwp\nLi6WH3/8US655BJp166dbNmyJUg0tc2ysjLzvLS0VK677rqo/V+9erXceuutUq9e\nPalZs6Zcfvnl0r59+yBRdvveeNpGAAHidFBVwdGjR2Xw4MGyYcOGSr1PxWDt2rWB\neX1/7969g/p/xRVXyPz5883z2bNnm6leGaVRo0ayefNmI4IzZsww7xk1apRcf/31\nsn37dvO6Xbt2yWWXXRa1H82bN5clS5bIb7/9JhUVFWb62WefycUXXxz3e+NpGwGE\nlBcPt+4pVgeVzOjura6/nblz50rjxo0D8yoeoS6uWbNmsmjRIvP8yJEjkpuba9Z5\n+vTp5vlVV11l6nT6ySefRO3HRRddJF9//bWcOnXKzB87dsw4ztq1a8f93njaRgAR\nkGrd/Upm9xSrg4LoTJkyRdq2bWtcZa1atcy2cfvtt8vEiRPjfm88bSOACEi17n6l\ngoOC6jUC4QZxdDcYAUwC0lVAkl084vkccYDJE0OEG8RRwUQAfUyyCIjbUbbqHJ2D\n1DABsRiBSIM4CGAS4HcBcTvKVp2jc159joh48huBaIM4CGA143cBcTvKVp2jc159\njsnyHUByggCmgYC4HWWrztE5rz7HZPkOAAH0LX4XkLBfvsuBgKoeQEjkIRiAAEIU\nGEQAQADTlnQdREDAAQGEtB1EYAABEEBI20EEBhAAAYTwH16KDyIwgAAIIAAAApi+\nMIgAgACmLQwiACCAaQuDCAAIYNrCIAIAAggAgAACACCAAAAIIAAAAggAgAACACCA\nAA5MmjTJ3PfhggsukJycHPn1118TurwdO3YEbgYejVatWsnUqVODynReywEBBIiL\nGTNmSOvWreW7776TsrIyGT16tPTv3z+hy5w+fbrk5ubG9qOoUUPatGkj5eXlZl77\nqGf1cCtNBBAgbm655RZZvnx5YP7w4cPy5ptvmuc7d+6UTp06Sd26dSUzM9PMW6I0\nbtw44xr1oPJ58+aZ8p9++kmysrLM2TUtWrSQdevWOS5TxU9FMFTonOZ1+vLLLweW\noXcwy8/PRwARQID4UbE6ePCgY13Pnj1lzJgxcujQISksLJTs7OyAKL3wwgvm/rQq\nTCqCitYXFBSY1+tFJlQEndDdX90NjlUAv//+e3n44YfNvE5VsBFABBAgbho0aCD7\n9+93rNMLRVjiqFM9n9oSJd0VDRUrzRCj5Yfh8r9IAqjce++95iIWLVu2NOdzI4AI\nIEDc9O7dW5YtWxZUtn79ejOtX79+QABV2FTgIomVimlpaWnE5YXL/+xt6jJDBXDC\nhAlm8GPw4MGOfQAEEKDSbN68WTp06CDbt283rk53bZ966ilT1717d7MLrOXPPPOM\nyfciCeD9998vzz33nHn9Bx984LgL7JT/KVaWeODAAenSpctZAqjOUS9goRkgAogA\nAnjG7NmzpXHjxuZisUOGDDHXSlR++OEHycjIMDlhx44dA7ldOAHUQZKuXbuaQRMd\nud24ceNZy3LK/xQrS1QXqYMwoQKoaF90kAYBRAABABBAAAAEEAAAAQQAQAABABBA\nAAAEEAAAAQQAQAABABBAAIDq5P8BWqsFoNlYa/UAAAAASUVORK5CYII="
TEST_DATA = {
    DN_ID: [
        (
            "mixture of enantiomers",
            "target",
            "1300",
            datetime.datetime(2022, 5, 6, 20, 6, 33),
            ASSAY1,
            "Analysis 1",
            "Good",
            "50",
            SAMPLE_IC50,
        ),
        (
            "single known stereoisomer",
            "target",
            "1400",
            datetime.datetime(2022, 5, 6, 20, 6, 33),
            ASSAY2,
            "Analysis 2",
            "Okay",
            "40",
            SAMPLE_IC50,
        ),
    ]
}
ANALYSES = {ASSAY1: ["Analysis 1"], ASSAY2: ["Analysis 2"]}


class AssayEmailerUtilsTestCase(BasechemTestCase):
    def setUp(self):
        super().setUp()

        # Update DN ID of series to match the test data
        self.series1.dn_id = DN_ID
        self.series1.save()

        # Create users
        self.project_lead = BasechemUserFactory()
        self.project_sub = BasechemUserFactory()

        # Initialize project data
        self.project = self.series1.project
        now = datetime.datetime.now()
        last_ping = now - datetime.timedelta(days=7)
        self.last_week = last_ping.replace(microsecond=0)
        self.few_days_ago = (now - datetime.timedelta(days=3)).replace(microsecond=0)
        self.many_months_ago = (now - datetime.timedelta(days=120)).replace(
            microsecond=0
        )
        self.project.assays = {
            "data_exp_id": 1000,
            "last_ping_sent": last_ping,
            "control": DN_ID,
            "assay_exp_info": {
                ASSAY1: [1200, self.few_days_ago],
                ASSAY2: [1100, self.few_days_ago],
            },
        }
        self.project.leads.add(self.project_lead)
        self.project.subscribers.add(self.project_sub)

        # Generate test df from TEST_DATA
        assay_data_as_series = [
            (k, *tup) for k, val in TEST_DATA.items() for tup in val
        ]
        self.test_data_df = pd.DataFrame(
            assay_data_as_series, columns=ASSAY_DATA_COLUMNS
        )
        self.test_data_df["Moltext"] = MOLTEXT

    @tag("local", "dtx")
    def test_process_new_data(self):
        """
        Test that `process_new_data` sends the correct emails and updates the
        project and series models accordingly
        """
        # Create a filtered test data that only includes 1 of the 2 assays
        assay_1_data_df = self.test_data_df[
            self.test_data_df["Assay_Name"] == ASSAY1
        ].copy()
        self.project.assays["assay_exp_info"][ASSAY1] = [1000, self.last_week]
        self.project.assays["assay_exp_info"][ASSAY2] = [900, self.last_week]

        with self.subTest("Ping only"):
            mail.outbox = []
            # Set last data email to less than MAX_HOURS_TO_WAIT ago
            minutes_ago = MAX_HOURS_TO_WAIT * 60 - 2
            time = datetime.datetime.now() - datetime.timedelta(minutes=minutes_ago)
            self.project.assays["last_ping_sent"] = time

            process_new_data(self.project, assay_1_data_df, ANALYSES)
            self.assertEqual(len(mail.outbox), 1)
            # Check ping
            self.assertEqual(
                mail.outbox[0].subject, f"Dotmatics got {self.project.code} data!"
            )
            message = f"New data for {self.project.code}: {ASSAY1} has been added to Dotmatics"
            self.assertIn(message, mail.outbox[0].body)
            self.assertEqual(mail.outbox[0].to, [self.project_lead.email])

        with self.subTest("Data only"):
            mail.outbox = []
            num_collections = len(Collection.objects.all())
            # Set last data email to more than MAX_HOURS_TO_WAIT ago
            minutes_ago = MAX_HOURS_TO_WAIT * 60 + 2
            time = datetime.datetime.now() - datetime.timedelta(minutes=minutes_ago)
            self.project.assays["last_ping_sent"] = time

            process_new_data(self.project, assay_1_data_df, ANALYSES)
            self.assertEqual(len(mail.outbox), 1)
            # Check data
            self.assertEqual(
                mail.outbox[0].subject, f"You've got {self.project} assay data!"
            )
            self.assertEqual(mail.outbox[0].to, [self.project_sub.email])
            # Check collection with assayed compounds created
            self.assertEqual(len(Collection.objects.all()), num_collections + 1)
            collection = Collection.objects.all().order_by("pk").last()
            self.assertEqual(collection.compound_occurrences.all().count(), 1)
            self.assertEqual(
                collection.compound_occurrences.first().compound.dn_id, "DN0017079"
            )
            self.assertEqual(self.project, collection.project)
            self.assertEqual(self.admin_user, collection.owner)
            self.assertEqual(collection.metadata["assays"], ANALYSES)

        with self.subTest("Ping and data"):
            mail.outbox = []
            # Reset the project to before we sent the impatient data email
            self.project.assays["data_exp_id"] = 1000
            num_collections = len(Collection.objects.all())
            collection = Collection.objects.all().order_by("pk").last()

            process_new_data(self.project, self.test_data_df, ANALYSES)
            self.assertEqual(collection.compound_occurrences.all().count(), 1)
            self.assertEqual(
                collection.compound_occurrences.first().compound.dn_id, "DN0017079"
            )
            self.assertEqual(len(mail.outbox), 2)

            # Check ping
            self.assertEqual(
                mail.outbox[0].subject, f"Dotmatics got {self.project.code} data!"
            )
            self.assertEqual(mail.outbox[0].to, [self.project_lead.email])

            # Check data
            self.assertEqual(
                mail.outbox[1].subject, f"You've got {self.project} assay data!"
            )
            self.assertEqual(mail.outbox[1].to, [self.project_sub.email])
            # Check collection with assayed compounds created
            self.assertEqual(len(Collection.objects.all()), num_collections + 1)
            created_collection = Collection.objects.all().order_by("pk").last()
            self.assertEqual(len(created_collection.compound_occurrences.all()), 1)
            self.assertEqual(collection.metadata["assays"], ANALYSES)

        with self.subTest("No new data"):
            mail.outbox = []
            num_collections = len(Collection.objects.all())
            process_new_data(self.project, self.test_data_df, ANALYSES)
            self.assertEqual(len(mail.outbox), 0)

            # Check no new collection has been added
            self.assertEqual(len(Collection.objects.all()), num_collections)

    @tag("local", "dtx")
    def test_send_new_data_assay_email(self):
        """
        Test that `send_new_data_assay_email` sends an email with assay data to the correct people
        """
        # Update the project model to reflect the new data
        # (`_assays_needing_pings` runs before `test_send_new_data_assay_email`)
        _ = _assays_needing_pings(self.project, self.test_data_df)

        self.assertEqual(len(mail.outbox), 0)
        send_new_data_assay_email(self.project, self.test_data_df, self.collection.pk)

        self.assertEqual(len(mail.outbox), 1)
        # Test subject and recipients
        self.assertEqual(
            mail.outbox[0].subject, f"You've got {self.project} assay data!"
        )
        self.assertEqual(mail.outbox[0].to, [self.project_sub.email])
        self.assertEqual(mail.outbox[0].from_email, settings.DEFAULT_FROM_EMAIL)

        # Test body
        self.assertEqual(mail.outbox[0].alternatives[0][1], "text/html")
        html = f"latest data for {self.project.code}"
        self.assertIn(html, mail.outbox[0].alternatives[0][0])

        # Test attachment
        today = datetime.datetime.now().strftime("%m_%d_%Y")
        attached_file_name = f"{today}_{self.project.code}_assayreport.pdf"
        self.assertEqual(mail.outbox[0].attachments[0][0], attached_file_name)

        # Test project updated
        self.project.refresh_from_db()
        self.assertEqual(self.project.assays["data_exp_id"], 1400)

    @tag("local", "dtx")
    def test__create_assay_html_table(self):
        """
        Test that `_create_assay_html_table` returns a string with an html table
        """
        html_table = _create_assay_html_table(self.test_data_df, self.project)

        self.assertIn('<table border="1" class="dataframe assay_table">', html_table)
        self.assertIn("<td>DN0017079</td>", html_table)
        self.assertIn(
            f"/media/project_{self.project.code}/images/{DN_ID}_structure.png",
            html_table,
        )
        self.assertIn(
            f"/media/project_{self.project.code}/images/{DN_ID}_{ASSAY2.replace(' ', '_')}_curve_1.png",
            html_table,
        )
        self.assertIn(f"<td>{ASSAY1_SHORT}</td>", html_table)
        self.assertIn(f"<td>{ASSAY2}</td>", html_table)

    def test__create_mol_image_file(self):
        """
        Test that `_create_mol_image_file` returns an html string with the correct image
        """
        img_tag = _create_mol_image_file(DN_ID, MOLTEXT, self.project)

        expected_img_tag = (
            f"<img src=/media/project_{self.project.code}/images/{DN_ID}_structure.png"
        )
        self.assertIn(expected_img_tag, img_tag)

    def test__create_curve_file(self):
        """
        Test that `_create_curve_file` returns an html string with the correct image
        """
        img_tag = _create_curve_file(DN_ID, SAMPLE_IC50, "assay_name", 1, self.project)

        expected_img_tag = f"<img src=/media/project_{self.project.code}/images/{DN_ID}_assay_name_curve_1.png"
        self.assertIn(expected_img_tag, img_tag)

    def test__delete_files(self):
        """
        Test that `_delete_files` removes all local assay emailer files for the project
        """
        images_dir = f"{settings.MEDIA_ROOT}/project_{self.project}/images"
        # Fake pdf path, no file should exist here
        pdf_path = f"{settings.MEDIA_ROOT}/project_{self.project}/*.pdf"

        image_tag = _create_curve_file(
            DN_ID, SAMPLE_IC50, "assay_name", 1, self.project
        )
        # Check image generated
        self.assertTrue(len(os.listdir(images_dir)) > 0)
        _delete_files(self.project, pdf_path)
        # Check image deleted
        self.assertEqual(os.listdir(images_dir), [])

    def test_send_new_data_ping(self):
        """
        Test that `send_new_data_ping` sends an email for the correct assays
        """
        with self.subTest("New ping"):
            self.assertEqual(len(mail.outbox), 0)
            send_new_data_ping(self.project, self.test_data_df)

            # Check email sent
            message = f"New data for {self.project.code}: {ASSAY1}, {ASSAY2} has been added to Dotmatics"
            self.assertEqual(len(mail.outbox), 1)
            self.assertIn(message, mail.outbox[0].body)
            self.assertEqual(
                mail.outbox[0].subject, f"Dotmatics got {self.project.code} data!"
            )
            self.assertEqual(mail.outbox[0].to, [self.project_lead.email])
            self.assertEqual(mail.outbox[0].from_email, settings.DEFAULT_FROM_EMAIL)

        with self.subTest("No new ping"):
            self.assertEqual(len(mail.outbox), 1)
            send_new_data_ping(self.project, self.test_data_df)
            self.assertEqual(len(mail.outbox), 1)

    def test__assays_needing_pings(self):
        """
        Test that `_assays_needing_pings` returns the set of assays that need pings
        """
        self.assertEqual(self.project.assays["assay_exp_info"][ASSAY1][0], 1200)
        self.assertEqual(self.project.assays["assay_exp_info"][ASSAY2][0], 1100)
        with self.subTest("New pings"):
            pings = _assays_needing_pings(self.project, self.test_data_df)

            # Check pings
            self.assertEqual(sorted(list(pings)), [ASSAY1, ASSAY2])
            # Check project model updated
            self.assertEqual(self.project.assays["assay_exp_info"][ASSAY1][0], 1300)
            self.assertEqual(self.project.assays["assay_exp_info"][ASSAY2][0], 1400)

        with self.subTest("No new pings"):
            pings = _assays_needing_pings(self.project, self.test_data_df)

            # Check pings
            self.assertEqual(sorted(list(pings)), [])
            # Check project model not updated
            self.assertEqual(self.project.assays["assay_exp_info"][ASSAY1][0], 1300)
            self.assertEqual(self.project.assays["assay_exp_info"][ASSAY2][0], 1400)

    def test_update_series_data(self):
        """
        Test that `update_series_data` updates the assay data for any series that have new data
        """
        with self.subTest("no new data"):
            series_df = self.test_data_df[self.test_data_df["DN_ID"] == ""]
            update_series_data(series_df)

            # Check series
            self.series1.refresh_from_db()
            self.assertEqual(self.series1.assay_data, {})

        with self.subTest("new data"):
            series_df = self.test_data_df[
                self.test_data_df["DN_ID"] == self.series1.dn_id
            ]
            update_series_data(series_df)

            # Check series
            expected_assay_data = {ASSAY1: "50", ASSAY2: "40"}
            self.series1.refresh_from_db()
            self.assertEqual(self.series1.assay_data, expected_assay_data)

    def test__all_assays_in(self):
        """
        Test `_all_assays_in` accurately reports if all the assays are in
        """
        with self.subTest("all assays in"):
            assays_in = _all_assays_in(self.project)
            self.assertTrue(assays_in)

        with self.subTest("not all assays in"):
            # Set ASSAY1 to not have new data
            self.project.assays["assay_exp_info"][ASSAY1][0] = 1000

            assays_in = _all_assays_in(self.project)
            self.assertFalse(assays_in)

    def test__any_assays_in(self):
        """
        Test `_any_assays_in` accurately reports if any of the assays are in
        """
        with self.subTest("some assays in"):
            assays_in = _any_assays_in(self.project)
            self.assertTrue(assays_in)

        with self.subTest("no assays in"):
            # Set both assays to not have new data
            self.project.assays["assay_exp_info"][ASSAY1][0] = 900
            self.project.assays["assay_exp_info"][ASSAY2][0] = 800

            assays_in = _any_assays_in(self.project)
            self.assertFalse(assays_in)

    def test__impatient_assays(self):
        """
        Test `_impatient_assays` returns True only when both of the following conditions are met:
        (1) enough time has elapsed since the last ping
        (2) not all assays have new data
        (3) at least one assay has new data
        """
        with self.subTest("time elapsed, all assays in"):
            impatient = _impatient_assays(self.project)
            self.assertFalse(impatient)

        with self.subTest("time elapsed, some new data"):
            # Set one assay to not have new data
            self.project.assays["assay_exp_info"][ASSAY1][0] = 900

            impatient = _impatient_assays(self.project)
            self.assertTrue(impatient)

        with self.subTest("time elapsed, no new data"):
            # Set both assays to not have new data
            self.project.assays["assay_exp_info"][ASSAY2][0] = 800

            impatient = _impatient_assays(self.project)
            self.assertFalse(impatient)

        # Set time to less than MAX_HOURS_TO_WAIT ago
        minutes_ago = MAX_HOURS_TO_WAIT * 60 - 1
        time = datetime.datetime.now() - datetime.timedelta(minutes=minutes_ago)
        self.project.assays["last_ping_sent"] = time

        with self.subTest("time not elapsed, no new data"):
            impatient = _impatient_assays(self.project)
            self.assertFalse(impatient)

        with self.subTest("time not elapsed, some new data"):
            # Set one assay to have new data
            self.project.assays["assay_exp_info"][ASSAY1][0] = 1200

            impatient = _impatient_assays(self.project)
            self.assertFalse(impatient)

        with self.subTest("time not elapsed, all assays in"):
            # Set both assays to have new data
            self.project.assays["assay_exp_info"][ASSAY2][0] = 1100
            impatient = _impatient_assays(self.project)
            self.assertFalse(impatient)

    def test_depreciate_old_assays(self):
        """
        Test that `depreciate_old_assays` removes old assays from the project
        """
        self.project.assays["assay_exp_info"]["assay3"] = [600, self.many_months_ago]
        self.assertEqual(len(self.project.assays["assay_exp_info"]), 3)
        depreciate_old_assays(self.project)
        self.assertEqual(len(self.project.assays["assay_exp_info"]), 2)
        self.assertNotIn("assay3", list(self.project.assays["assay_exp_info"].keys()))

    def test__clean_assay_names(self):
        """
        Test that `clean_assay_names` reduces assay names to three words
        """
        assay_names = self.test_data_df["Assay_Name"].unique().tolist()
        self.assertEqual(sorted(assay_names), [ASSAY1, ASSAY2])

        self.test_data_df = _clean_assay_names(self.test_data_df)
        assay_names = self.test_data_df["Assay_Name"].unique().tolist()
        self.assertEqual(sorted(assay_names), [ASSAY1_SHORT, ASSAY2])

    def test_update_assay_results_from_new_data_df(self):
        """
        Test `update_assay_results_from_new_data_df` updates measured_data.assay_results for Compounds
        in the given dataframe of assay data
        """
        with self.subTest("Empty dataframe"):
            empty_df = pd.DataFrame(columns=ASSAY_DATA_COLUMNS)
            update_assay_results_from_new_data_df(empty_df)

        with self.subTest("Compound updated"):
            c = CompoundFactory(dn_id="DN0017079")
            self.assertEqual(c.measured_data, {})
            update_assay_results_from_new_data_df(self.test_data_df)
            c.refresh_from_db()
            expected_mesasured_data = {
                "assay_results": {
                    ASSAY1: {"Analysis 1": "50"},
                    ASSAY2: {"Analysis 2": "40"},
                }
            }
            self.assertEqual(c.measured_data, expected_mesasured_data)
