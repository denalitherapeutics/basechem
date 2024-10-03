from django.conf import settings
from django.test import tag
from django.urls import reverse

from ...storage import save_media_file
from ..base import MniCommonTestCase


class DownloadFileTestCase(MniCommonTestCase):
    """
    Test the `DownloadFileView` for both static and media files
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.url = reverse(
            "download", args=["static", "mni_common/images/company-logo.png"]
        )

    def test_unauthenticated(self):
        """
        User must be logged in to use the view
        """
        self.client.logout()
        response = self.client.get(self.url)
        self.assertRedirects(response, f"{reverse('login')}?next={self.url}")

    @tag("local")
    def test_download_static(self):
        """
        Test downloading a static file. This test must be run locally because `collectstatic`
        is not run when running django tests
        """
        response = self.client.get(self.url)
        self.assertEqual(response.status_code, 200)
        self.assertIn("company-logo.png", response["Content-Disposition"])

    def test_download_media(self):
        """
        Test downloading a media file
        """
        filepath = "a_dir/a_file.png"
        # Save a file to the MEDIA directory to download
        with open(
            f"{settings.PROJECT_NAME}/mni_common/static/mni_common/images/company-logo.png",
            "rb",
        ) as fp:
            save_media_file(filepath, fp)

        response = self.client.get(reverse("download", args=["media", filepath]))
        self.assertEqual(response.status_code, 200)
        self.assertIn("a_file.png", response["Content-Disposition"])
