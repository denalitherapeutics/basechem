from django.conf import settings
from django.test import TestCase, override_settings

from basechem.common.context_processors import sentry_dsn


class ContextProcessorsTest(TestCase):
    """
    The Javascript Sentry SDK is initialized in base.html based on the value of the SENTRY_DSN
    context variable. These tests check that the SENTRY_DSN context variable is set correctly
    based on whether SENTRY_DSN is specified in django settings.
    """

    @override_settings(SENTRY_DSN="abc123")
    def test_with_sentry_dsn(self):
        env = getattr(settings, "ENVIRONMENT")
        result = sentry_dsn(None)
        self.assertEqual(result, {"ENVIRONMENT": env, "SENTRY_DSN": "abc123"})

    @override_settings(SENTRY_DSN=None)
    def test_without_sentry_dsn(self):
        env = getattr(settings, "ENVIRONMENT")
        delattr(settings, "SENTRY_DSN")
        result = sentry_dsn(None)
        self.assertEqual(result, {"ENVIRONMENT": env, "SENTRY_DSN": False})
