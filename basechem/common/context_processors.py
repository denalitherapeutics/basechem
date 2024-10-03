from django.conf import settings


def sentry_dsn(request):
    """
    Provide SENTRY_DSN and ENVIRONMENT in templates so we can integrate Sentry for Javascript errors
    """
    return {
        "SENTRY_DSN": getattr(settings, "SENTRY_DSN", False),
        "ENVIRONMENT": getattr(settings, "ENVIRONMENT", False),
    }
