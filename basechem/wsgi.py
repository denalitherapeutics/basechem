"""
WSGI config for ProteinRequestSystem project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.9/howto/deployment/wsgi/
"""

from django.core.wsgi import get_wsgi_application

from . import load_env

load_env.load_env()
application = get_wsgi_application()
