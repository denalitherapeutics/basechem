from django.core import mail
from django.urls import reverse

from basechem.common.tests.base import BasechemTestCase, BasechemViewTestMixin
from basechem.main.constants import ALIGN


class BasechemUserTestCase(BasechemViewTestMixin, BasechemTestCase):
    def test_str(self):
        """Test for string representation."""
        self.assertEqual(str(self.user), self.user.get_full_name())

    def test_get_full_name(self):
        with self.subTest("No first_name or last_name"):
            self.user.first_name = ""
            self.user.last_name = ""
            self.user.save()
            self.assertEqual(self.user.get_full_name(), "")

        with self.subTest("Has a first_name, but no last_name"):
            self.user.first_name = "Alex"
            self.user.last_name = ""
            self.user.save()
            self.assertEqual(self.user.get_full_name(), "Alex")

        with self.subTest("No first_name, but has a last_name"):
            self.user.first_name = ""
            self.user.last_name = "Smith"
            self.user.save()
            self.assertEqual(self.user.get_full_name(), "Smith")

        with self.subTest("Has a first_name and a last_name"):
            self.user.first_name = "Alex"
            self.user.last_name = "Smith"
            self.user.save()
            self.assertEqual(self.user.get_full_name(), "Alex Smith")

    def test_get_short_name(self):
        with self.subTest("No first_name or last_name"):
            self.user.first_name = ""
            self.user.last_name = ""
            self.user.save()
            self.assertEqual(self.user.get_short_name(), "")

        with self.subTest("Has a first_name, but no last_name"):
            self.user.first_name = "Alex"
            self.user.last_name = ""
            self.user.save()
            self.assertEqual(self.user.get_short_name(), "Alex")

        with self.subTest("No first_name, but has a last_name"):
            self.user.first_name = ""
            self.user.last_name = "Smith"
            self.user.save()
            self.assertEqual(self.user.get_short_name(), "Smith")

        with self.subTest("Has a first_name and a last_name"):
            self.user.first_name = "Alex"
            self.user.last_name = "Smith"
            self.user.save()
            self.assertEqual(self.user.get_short_name(), "Alex")

    def test_get_email_name(self):
        with self.subTest("Has Denali email address"):
            self.user.email = "anemployee@dnli.com"
            self.assertEqual(self.user.get_email_name(), "anemployee")

        with self.subTest("Has other email address"):
            self.user.email = "outsourced@example.com"
            self.assertEqual(self.user.get_email_name(), "outsourced")

        with self.subTest("User has blank email"):
            self.user.email = ""
            self.assertEqual(self.user.get_email_name(), "")

        with self.subTest("User has meaningless email address"):
            self.user.email = "this is not an email"
            self.assertEqual(self.user.get_email_name(), "this is not an email")

    def test_add_easter_egg_point(self):
        http_request = self.client.get(
            reverse(ALIGN, kwargs={"collection_id": self.collection_one_cpd.id})
        ).wsgi_request
        egg_name = "Hiking"

        with self.subTest("First point"):
            self.assertEqual(self.user.easter_egg_points, {})

            self.user.add_easter_egg_point(http_request, egg_name)

            self.user.refresh_from_db()
            self.assertEqual(self.user.easter_egg_points, {egg_name: 1})
            self.assertEqual(len(mail.outbox), 0)  # Check no email sent

        with self.subTest("Not enough to win"):
            self.user.add_easter_egg_point(http_request, egg_name)

            self.user.refresh_from_db()
            self.assertEqual(self.user.easter_egg_points, {egg_name: 2})
            self.assertEqual(len(mail.outbox), 0)  # Check no email sent

        with self.subTest("Discover egg!"):
            self.user.add_easter_egg_point(http_request, egg_name)

            self.user.refresh_from_db()
            self.assertEqual(self.user.easter_egg_points, {egg_name: 3})
            self.assertEqual(len(mail.outbox), 1)  # Check email sent

        with self.subTest("After egg discovery"):
            self.user.add_easter_egg_point(http_request, egg_name)

            self.user.refresh_from_db()
            self.assertEqual(self.user.easter_egg_points, {egg_name: 4})
            self.assertEqual(len(mail.outbox), 1)  # Check no email sent
