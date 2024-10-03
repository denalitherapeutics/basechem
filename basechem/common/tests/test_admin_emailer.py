from django.conf import settings
from django.core import mail

from basechem.common.admin_emailer import send_easter_egg_email
from basechem.common.tests.base import BasechemTestCase


class EmailTestCase(BasechemTestCase):
    def test_send_easter_egg_email(self):
        """
        Test that `send_easter_egg_email` sends an email
        """
        self.assertEqual(len(mail.outbox), 0)

        send_easter_egg_email(self.owner, "Test Egg")

        message = f"{self.owner.first_name} has found the 'Test Egg' Easter egg!"
        self.assertEqual(len(mail.outbox), 1)
        self.assertEqual(mail.outbox[0].body, message)
        self.assertEqual(mail.outbox[0].subject, "Basechem Easter Egg Found!")
        self.assertEqual(mail.outbox[0].to, [settings.ADMIN_EMAIL, self.owner.email])
        self.assertEqual(mail.outbox[0].from_email, settings.DEFAULT_FROM_EMAIL)
