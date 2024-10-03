import datetime
import os

from django.conf import settings
from django.core import mail
from django.test import tag
from rdkit import Chem

from basechem.common.tests.base import BasechemNoMockTestCase
from basechem.main.constants import ALIGN, DOCK, MAX_N_CONFS_DISPLAY, TORSION
from basechem.main.management.commands.bulk_analysis import run_bulk_analysis
from basechem.main.models.collection_models import Collection, collection_files_path


@tag("local")
class TasksTestCase(BasechemNoMockTestCase):
    def test_run_bulk_analysis(self):
        """
        Test bulk analysis returns the expected results
        """
        existing_project_code = self.project1.code

        with self.subTest("Test bulk align"):
            self.assertEqual(len(mail.outbox), 0)
            expected_file_prefix = f"bulk_{ALIGN}_{existing_project_code}_{datetime.datetime.now().strftime('%m%d%Y')}"
            run_bulk_analysis(
                f"{settings.PROJECT_ROOT}/basechem/main/tests/testdata/7eri_ligand.sdf",
                existing_project_code,
                ALIGN,
            )

            # Check bulk collection is created
            collection = Collection.objects.latest("created_on")
            self.assertEqual(len(collection.compounds()), 1)
            self.assertEqual(expected_file_prefix, collection.name)

            # Check mol results from the job exist
            expected_file_path = collection_files_path(
                collection, f"{expected_file_prefix}.sdf", local=True
            )
            self.assertTrue(os.path.exists(expected_file_path))
            suppl = Chem.SDMolSupplier(expected_file_path)
            self.assertTrue(4 <= len([mol for mol in suppl]) <= MAX_N_CONFS_DISPLAY)

            # Check admins are notified via email of completion
            self.assertEqual(len(mail.outbox), 1)
            self.assertIn("successfully completed", mail.outbox[-1].body)

        with self.subTest("Test bulk dock"):
            self.assertEqual(len(mail.outbox), 1)
            expected_file_prefix = f"bulk_{DOCK}_{existing_project_code}_{datetime.datetime.now().strftime('%m%d%Y')}"
            run_bulk_analysis(
                f"{settings.PROJECT_ROOT}/basechem/main/tests/testdata/7eri_ligand.sdf",
                existing_project_code,
                DOCK,
            )

            # Check bulk collection is created
            collection = Collection.objects.latest("created_on")
            self.assertEqual(len(collection.compounds()), 1)
            self.assertEqual(expected_file_prefix, collection.name)

            # Check mol results from the job exist
            expected_file_path = collection_files_path(
                collection, f"{expected_file_prefix}.sdf", local=True
            )
            self.assertTrue(os.path.exists(expected_file_path))
            suppl = Chem.SDMolSupplier(expected_file_path)
            self.assertTrue(4 <= len([mol for mol in suppl]) <= 8)

            # Check admins are notified via email of completion
            self.assertEqual(len(mail.outbox), 2)
            self.assertIn("successfully completed", mail.outbox[-1].body)

        with self.subTest("Test bad analysis type"):
            self.assertEqual(len(mail.outbox), 2)
            expected_file_prefix = f"bulk_{TORSION}_{existing_project_code}_{datetime.datetime.now().strftime('%m%d%Y')}"
            run_bulk_analysis(
                "basechem/main/tests/testdata/7eri_ligand.sdf",
                existing_project_code,
                TORSION,
            )

            # Check new bulk collection is not created
            collection = Collection.objects.latest("created_on")
            self.assertEqual(len(collection.compounds()), 1)
            self.assertNotEqual(expected_file_prefix, collection.name)

            # Check mol results from the job exist
            expected_file_path = collection_files_path(
                collection, f"{expected_file_prefix}.sdf", local=True
            )
            self.assertFalse(os.path.exists(expected_file_path))

            # Check no email is sent
            self.assertEqual(len(mail.outbox), 2)
