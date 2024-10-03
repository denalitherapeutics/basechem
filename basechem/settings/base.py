import os
import sys

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, os.pardir))
PROJECT_NAME = BASE_DIR.split("/")[-1]

ENVIRONMENT = os.environ.get("ENVIRONMENT", "local")

SAML_IP_METADATA_URL = os.environ.get("SAML_IP_METADATA_URL", "")
SECRET_KEY = os.environ.get("SECRET_KEY", "SECRET_KEY")

BASE_URL = os.environ.get("BASE_URL", "https://localhost:8000")

ADMIN_USER = os.environ.get("ADMIN_USER", "admin")
ADMIN_PASSWORD = os.environ.get("ADMIN_PASSWORD", "admin")
ADMIN_EMAIL = os.environ.get("ADMIN_EMAIL", "")

# Set necessary variables from .env for local development
if ENVIRONMENT == "local":
    DEBUG = True
    CRISPY_FAIL_SILENTLY = False

    EMAIL_BACKEND = "django.core.mail.backends.console.EmailBackend"
    ALLOWED_HOSTS = ["*"]
    AWS_STORAGE_BUCKET_NAME = None

ADMINS = [
    (os.environ.get("DJANGO_ADMIN_NAME", ""), os.environ.get("DJANGO_ADMIN_EMAIL", ""))
]

# Application definition

INSTALLED_APPS = [
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "django.contrib.admin",
    "django.contrib.humanize",
    "django.contrib.sitemaps",
    "crispy_forms",
    "formtools",
    "storages",
    "django_q",
    "basechem.sso",
    "basechem.users",
    "basechem.main",
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.locale.LocaleMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

if ENVIRONMENT == "local":
    INSTALLED_APPS += ("debug_toolbar",)
    MIDDLEWARE += ("debug_toolbar.middleware.DebugToolbarMiddleware",)

    INTERNAL_IPS = ["127.0.0.1"]
    import socket

    # tricks to have debug toolbar when developing with docker
    ip = socket.gethostbyname(socket.gethostname())
    INTERNAL_IPS += [ip[:-1] + "1"]

ROOT_URLCONF = "basechem.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [
            os.path.join(BASE_DIR, "templates"),
            os.path.join(BASE_DIR, "mni_common", "templates"),
        ],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.contrib.auth.context_processors.auth",
                "django.template.context_processors.debug",
                "django.template.context_processors.i18n",
                "django.template.context_processors.tz",
                "django.template.context_processors.request",
                "django.contrib.messages.context_processors.messages",
                "dealer.contrib.django.context_processor",
            ],
            "libraries": {
                "mni_common_tags": "basechem.mni_common.templatetags.mni_common_tags"
            },
        },
    }
]

WSGI_APPLICATION = "basechem.wsgi.application"

# Database
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": os.environ.get("DB_NAME", "basechem"),
        "USER": os.environ.get("DB_USER", "basechem"),
        "PASSWORD": os.environ.get("DB_PASSWORD", "basechem"),
        "HOST": os.environ.get("DB_HOST", "db"),
        "PORT": os.environ.get("DB_PORT", "5432"),
    }
}

# Absolute filesystem path to the directory that will hold user-uploaded files.
MEDIA_ROOT = os.path.join(PROJECT_ROOT, "public", "media")
MEDIA_URL = "/media/"

LOGGING_DIR = os.path.join(PROJECT_ROOT, "log")

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "filters": {"require_debug_false": {"()": "django.utils.log.RequireDebugFalse"}},
    "formatters": {
        "basic": {"format": "%(asctime)s %(name)-20s %(levelname)-8s %(message)s"}
    },
    "handlers": {
        "mail_admins": {
            "level": "ERROR",
            "filters": ["require_debug_false"],
            "class": "django.utils.log.AdminEmailHandler",
        },
        "console": {
            "level": "INFO",
            "class": "logging.StreamHandler",
            "formatter": "basic",
        },
        "debug_file": {
            "level": "INFO",
            "class": "logging.FileHandler",
            "filename": os.path.join(LOGGING_DIR, "main_debug.log"),
        },
        "analytics_file": {
            "level": "INFO",
            "class": "logging.FileHandler",
            "filename": os.path.join(LOGGING_DIR, "analytics.log"),
        },
    },
    "loggers": {
        "django.request": {
            "handlers": ["mail_admins"],
            "level": "ERROR",
            "propagate": True,
        },
        "django.security": {
            "handlers": ["mail_admins"],
            "level": "ERROR",
            "propagate": True,
        },
        "django": {
            "handlers": ["debug_file"],
            "level": "DEBUG",
            "propagate": True,
        },
        "analytics": {
            "handlers": ["analytics_file"],
            "level": "INFO",
            "propagate": True,
        },
    },
    "root": {"handlers": ["console"], "level": "INFO"},
}

Q_CLUSTER = {
    "name": "default",  # name of the queue from Q_CLUSTER_CONFIG_LIST that should be used by default
    "orm": "default",
    "workers": 2,
    "save_limit": 5000,
    "ack_failures": True,
    "has_replica": True,
    "timeout": 43200,  # timeout after 12 hours
    "retry": 43201,  # Never retry (because retry > timeout and ack_failures is True)
    "max_attempts": 1,
    "ALT_CLUSTERS": {
        "slow": {
            "timeout": 604800,  # timeout after 7 days
            "retry": 604801,
            "workers": 1,
        },
        "fast": {
            "timeout": 3600,  # timeout after 1 hour
            "retry": 3601,
            "workers": 2,
            "cpu_affinity": 3,
        },
    },
}


# Internationalization
LANGUAGE_CODE = "en-us"
LOCALE_PATHS = (os.path.join(PROJECT_ROOT, "locale"),)

TIME_ZONE = "America/Los_Angeles"

USE_I18N = True
USE_L10N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.9/howto/static-files/
# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/home/media/media.lawrence.com/static/"
STATIC_ROOT = os.path.join(PROJECT_ROOT, "public", "static")

# URL prefix for static files.
# Example: "http://media.lawrence.com/static/"
STATIC_URL = "basechem/static/"

# Additional locations of static files
STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "static"),
    os.path.join(BASE_DIR, "mni_common", "static"),
    os.path.join(BASE_DIR, "common", "ketcher"),
)

# https://docs.djangoproject.com/en/1.9/topics/auth/passwords/#password-validation
AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator"
    },
    {"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator"},
    {"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator"},
    {"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator"},
]

AUTH_USER_MODEL = "users.BasechemUser"

CRISPY_TEMPLATE_PACK = "bootstrap4"

LOGIN_URL = "login"
LOGIN_REDIRECT_URL = "homepage"
LOGOUT_REDIRECT_URL = "login"

# Default pagination is set to 25 (result per page).
DEFAULT_PAGINATION = 25

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"


SENTRY_DSN = os.environ.get("SENTRY_DSN", None)

if SENTRY_DSN:
    # Sentry config
    # https://docs.sentry.io/platforms/python/django/
    # Note: this integrates with logging without our having to do anything else.
    import sentry_sdk
    from sentry_sdk.integrations.django import DjangoIntegration

    sentry_sdk.init(
        dsn=SENTRY_DSN,
        integrations=[DjangoIntegration()],
        environment=ENVIRONMENT,
        send_default_pii=True,
        in_app_include="basechem",
    )

    # Put public DSN into context for templates to access
    TEMPLATES[0]["OPTIONS"]["context_processors"] += (
        "basechem.common.context_processors.sentry_dsn",
    )


# If the django test suite is running, set password hashers
if "test" in sys.argv:
    PASSWORD_HASHERS = (
        "django.contrib.auth.hashers.SHA1PasswordHasher",
        "django.contrib.auth.hashers.MD5PasswordHasher",
    )
    LOGGING["root"]["handlers"] = []
    Q_CLUSTER["sync"] = True

if os.environ.get("GITHUB_WORKFLOW"):
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.postgresql",
            "NAME": "basechem",
            "USER": "basechem",
            "PASSWORD": "basechem",
            "HOST": "postgres",
            "PORT": "5432",
        }
    }

#######################
# "External" settings #
#######################

DTX_HOST = os.environ.get("DTX_HOST", "")
DTX_USER = os.environ.get("DTX_USER", "")
DTX_PASSWORD = os.environ.get("DTX_PASSWORD", "")

ESP_DIR = f"{BASE_DIR}/common/ESP_DNN"
TOKLAT_DIR = os.environ.get("TOKLAT_DIR", "/opt/toklat")
MMPDB_DIR = os.environ.get("MMPDB_DIR", "/opt/mmpdb-3.1")

INDUCTIVE_CUSTOMER_ID = os.environ.get("INDUCTIVE_CUSTOMER_ID", "")
INDUCTIVE_API_KEY = os.environ.get("INDUCTIVE_API_KEY", "")
INDUCTIVE_VERSION = os.environ.get("INDUCTIVE_VERSION", "dev")
INDUCTIVE_BIO_ENABLED = (
    INDUCTIVE_CUSTOMER_ID and INDUCTIVE_API_KEY and INDUCTIVE_VERSION
)
MAYACHEMTOOLS_DIR = os.environ.get("MAYACHEMTOOLS_DIR", "/opt/mayachemtools")
MAYACHEMTOOLS_ENV = os.environ.get("MAYACHEMTOOLS_ENV", "/usr/local/envs/cenv")
MAYACHEMTOOLS_CONDA_EXEC_PATH = os.environ.get(
    "MAYACHEMTOOLS_CONDA_EXEC_PATH", "/usr/local/bin/conda"
)

SLURM_REST_API_HOST = os.environ.get("SLURM_REST_API_HOST", "")
SLURM_SSH_KEY_PATH = os.environ.get("SLURM_SSH_KEY_PATH", "")
if SLURM_REST_API_HOST or "test" in sys.argv:
    SLURM_REST_API_BASE_URL = f"http://{SLURM_REST_API_HOST}:6899/slurm/v0.0.36"
    SLURM_SHARED_FILES_TMP_DIR = os.environ.get("SLURM_SHARED_FILES_TMP_DIR", "")
    SHARED_FILES_TMP_DIR_FROM_CONTAINER = (
        "/opt/shared_files/local"
        if ENVIRONMENT == "local"
        else SLURM_SHARED_FILES_TMP_DIR
    )
