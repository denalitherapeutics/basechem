# flake8: noqa
import os

from boto3.session import Session

from .base import *  # noqa

# must set .env ENVIRONMENT=something other than 'local' to use this file
DEBUG = os.environ.get("DEBUG", False)
DOMAIN = os.environ.get("DOMAIN", "*")
WEBSERVER_ROOT = os.environ.get("WEBSERVER_ROOT", "/home/app/web/")

if ENVIRONMENT == "prod":
    Q_CLUSTER["workers"] = 8  # default queue
    Q_CLUSTER["ALT_CLUSTERS"]["fast"][
        "workers"
    ] = 8  # propcalc and new_collection only rn
    Q_CLUSTER["ALT_CLUSTERS"]["slow"]["workers"] = 4  # torsion only rn
    INDUCTIVE_VERSION = "latest"

MIDDLEWARE += ("basechem.common.middleware.ExceptionMiddleware",)

# Set up for django-storages S3 integration
AWS_ACCESS_KEY_ID = os.environ.get("AWS_ACCESS_KEY_ID", "")
AWS_SECRET_ACCESS_KEY = os.environ.get("AWS_SECRET_ACCESS_KEY", "")
AWS_STORAGE_BUCKET_NAME = os.environ.get("AWS_STORAGE_BUCKET_NAME", "")
AWS_SESSION_TOKEN = os.environ.get("AWS_SESSION_TOKEN", "")

DEFAULT_FILE_STORAGE = "basechem.mni_common.storage_backends.MediaStorage"
STATICFILES_STORAGE = "basechem.mni_common.storage_backends.StaticStorage"
AWS_S3_REGION_NAME = os.environ.get("AWS_S3_REGION_NAME", "us-west-2")
AWS_S3_CUSTOM_DOMAIN = "%s.s3.%s.amazonaws.com" % (
    AWS_STORAGE_BUCKET_NAME,
    AWS_S3_REGION_NAME,
)
AWS_DEFAULT_ACL = "public-read"

##########################################################################################
#                                    Logger Settings                                     #
##########################################################################################

LOGGING_DIR = os.path.join(WEBSERVER_ROOT, "log")

# Configure AWS Logging, only on the actual "test" and "prod" sites (not staging containers)
if ENVIRONMENT in ["test", "prod"]:
    # Get a boto3 session to use for AWS logging
    logger_boto3_session = Session(
        aws_access_key_id=AWS_ACCESS_KEY_ID,
        aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
        region_name=AWS_S3_REGION_NAME,
    )
    # AWS log group that contains the "analytics" and "django" log streams
    log_group = "/basechem/mni-test/basechem"
    if ENVIRONMENT == "prod":
        log_group = "/basechem/mni-prod/basechem"

    # Configure django error logging
    LOGGING["formatters"]["aws_error_formatter"] = {
        "format": "%(asctime)s [%(levelname)s] %(message)s [%(pathname)s:%(lineno)d]",
        "datefmt": "%Y-%m-%d %H:%M:%S",
    }
    LOGGING["handlers"]["django_aws"] = {
        "level": "INFO",
        "class": "watchtower.CloudWatchLogHandler",
        "boto3_client": logger_boto3_session.client("logs"),
        "log_group": log_group,
        "stream_name": "django",
        "formatter": "aws_error_formatter",
    }
    LOGGING["loggers"]["django"]["handlers"] = ["django_aws"]

    # Configure analytics logging
    LOGGING["formatters"]["analytics_formatter"] = {
        "format": "%(message)s",
        "datefmt": "%Y-%m-%d %H:%M:%S",
    }
    LOGGING["handlers"]["analytics_aws"] = {
        "level": "INFO",
        "class": "watchtower.CloudWatchLogHandler",
        "boto3_client": logger_boto3_session.client("logs"),
        "log_group": log_group,
        "stream_name": "analytics",
        "formatter": "analytics_formatter",
    }
    LOGGING["loggers"]["analytics"]["handlers"] = ["analytics_aws"]

##########################################################################################
#                                     Email settings                                     #
##########################################################################################

EMAIL_HOST = os.environ.get("EMAIL_HOST", "mail.smtp2go.com")
EMAIL_HOST_USER = os.environ.get("EMAIL_HOST_USER", "")
EMAIL_HOST_PASSWORD = os.environ.get("EMAIL_HOST_PASSWORD", "")
EMAIL_USE_TLS = os.environ.get("EMAIL_USE_TLS", True)
EMAIL_USE_SSL = os.environ.get("EMAIL_USE_SSL", False)
# use TLS or SSL, not both:
assert not (EMAIL_USE_TLS and EMAIL_USE_SSL)
if EMAIL_USE_TLS:
    default_smtp_port = 587
elif EMAIL_USE_SSL:
    default_smtp_port = 465
else:
    default_smtp_port = 25
EMAIL_PORT = os.environ.get("EMAIL_PORT", default_smtp_port)
EMAIL_SUBJECT_PREFIX = "[Basechem %s] " % ENVIRONMENT.title()
DEFAULT_FROM_EMAIL = os.environ.get("DEFAULT_FROM_EMAIL", "")
SERVER_EMAIL = DEFAULT_FROM_EMAIL

##########################################################################################
#                                     Other settings                                     #
##########################################################################################

CSRF_COOKIE_SECURE = bool(os.environ.get("CSRF_COOKIE_SECURE", True))
SESSION_COOKIE_SECURE = bool(os.environ.get("SESSION_COOKIE_SECURE", True))
SESSION_COOKIE_HTTPONLY = bool(os.environ.get("SESSION_COOKIE_HTTPONLY", True))

ALLOWED_HOSTS = [DOMAIN]

# Use template caching on deployed servers
for backend in TEMPLATES:
    if backend["BACKEND"] == "django.template.backends.django.DjangoTemplates":
        default_loaders = ["django.template.loaders.filesystem.Loader"]
        if backend.get("APP_DIRS", False):
            default_loaders.append("django.template.loaders.app_directories.Loader")
            # Django gets annoyed if you both set APP_DIRS True and specify your own loaders
            backend["APP_DIRS"] = False
        loaders = backend["OPTIONS"].get("loaders", default_loaders)
        for loader in loaders:
            if (
                len(loader) == 2
                and loader[0] == "django.template.loaders.cached.Loader"
            ):
                # We're already caching our templates
                break
        else:
            backend["OPTIONS"]["loaders"] = [
                ("django.template.loaders.cached.Loader", loaders)
            ]
