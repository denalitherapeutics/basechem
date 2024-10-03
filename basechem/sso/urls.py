from django.urls import path

from . import views

urlpatterns = [
    path(r"", views.homepage_gateway, name="login"),
    path(r"saml_login/", views.saml_login, name="sso_login"),
    path(r"saml_success/", views.saml_success, name="sso_success"),
]
