import tempfile

from django.test import TestCase, override_settings

from .factories import TEST_PASSWORD, UserFactory


@override_settings(MEDIA_ROOT=tempfile.mkdtemp())
class MniCommonTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.TEST_PASSWORD = TEST_PASSWORD
        cls.user = UserFactory()

    def setUp(self):
        super().setUp()
        self.client.login(username=self.user.username, password=self.TEST_PASSWORD)
