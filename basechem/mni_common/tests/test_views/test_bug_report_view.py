from unittest.mock import patch

from django.conf import settings
from django.core import mail
from django.core.files.uploadedfile import SimpleUploadedFile
from django.urls import reverse

from ..base import MniCommonTestCase


class SendBugReportTestCase(MniCommonTestCase):
    def setUp(self):
        super().setUp()
        self.url = reverse("send_bug_report")
        self.form_data = {
            "description": "I found a bug!",
            "url": "https://testurl.com/basechem",
        }
        with open(
            f"{settings.PROJECT_NAME}/mni_common/static/mni_common/images/company-logo.png",
            "rb",
        ) as fp:
            self.file_a = SimpleUploadedFile("file_a.xlsx", fp.read())
        with open(
            f"{settings.PROJECT_NAME}/mni_common/static/mni_common/images/user-white.png",
            "rb",
        ) as fp:
            self.file_b = SimpleUploadedFile("file_b.xlsx", fp.read())

    def test_unauthenticated(self):
        """
        User must be logged in to use this view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    def test_no_files(self):
        """
        Test sending a bug report with no files succeeds
        """
        response = self.client.post(self.url, data=self.form_data)
        self.assertEqual(response.json()["errors"], {})
        self.assertEqual(len(mail.outbox), 1)
        email = mail.outbox[-1]
        # Check email contents
        self.assertEqual(len(email.to), len(settings.ADMINS))
        self.assertIn("https://testurl.com/basechem", email.body)
        self.assertIn("I found a bug!", email.body)
        self.assertEqual(email.attachments, [])

    def test_has_files(self):
        """
        Test sending a bug report with files succeeds
        """
        self.form_data["files"] = [self.file_a, self.file_b]
        response = self.client.post(self.url, data=self.form_data)
        self.assertEqual(response.json()["errors"], {})
        self.assertEqual(len(mail.outbox), 1)
        # Check email contents
        email = mail.outbox[-1]
        self.assertEqual(len(email.to), len(settings.ADMINS))
        self.assertIn("https://testurl.com/basechem", email.body)
        self.assertIn("I found a bug!", email.body)
        self.assertEqual(len(email.attachments), 2)
        self.assertEqual(email.attachments[0][0], self.file_a.name)
        self.assertEqual(email.attachments[1][0], self.file_b.name)

    def test_failed_bug_report(self):
        """
        Test that the user receives an error message and the admins receive an email if
        the bug report fails to send.
        """
        self.form_data["files"] = [self.file_a]

        # Enforce an unexpected failure by adding a side effect to a random function in the try block
        with patch("mimetypes.guess_type", autospec=True) as mock_obj:
            mock_obj.side_effect = Exception("Unexpected Error!!!")
            response = self.client.post(self.url, data=self.form_data)

        expected_errors = {
            "__all__": [
                f"Something went wrong - our development team will not receive this report. Please email {settings.ADMINS[0][1]} with your bug."
            ]
        }
        self.assertEqual(response.json()["errors"], expected_errors)
        self.assertEqual(len(mail.outbox), 1)
        # Check email contents
        email = mail.outbox[-1]
        self.assertEqual(len(email.to), len(settings.ADMINS))
        self.assertIn("tried to file a bug report", email.body)
        self.assertIn("Unexpected Error!!!", email.body)
        self.assertEqual(len(email.attachments), 0)
