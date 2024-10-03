from unittest import mock

from basechem.main.tests.factories import BasechemUserFactory


class BasechemSSOMixin:
    """
    A mixin that sets SSO-related mock values:
      - self.mock_get_saml_auth_return_value for the expected return value of
        basechem.sso.views.get_saml_auth()
      - self.mock_get_saml_settings_return_value for the expected return value of
        basechem.sso.utils.get_saml_settings()
    """

    def setUp(self):
        super().setUp()
        self.user1 = BasechemUserFactory(
            username="requestor1@example.com",
            email="requestor1@example.com",
            first_name="requestor",
            last_name="one",
        )
        self.client.login(username=self.user1.username, password="testpassword")

        # Mock the get_saml_auth() function in basechem.sso.views.py
        self.patch_saml_auth_object = mock.patch("basechem.sso.views.get_saml_auth")
        self.mock_get_saml_auth = self.patch_saml_auth_object.start()
        self.mock_idp_url = (
            "https://dev-1.okta.com/app/appidhere/something/sso/saml?"
            "SAMLRequest=samlhere&RelayState=http%3A%2F%2Flocalhost%3A8000%2Fsaml_login%2F"
        )
        mock_saml_auth_object = mock.Mock()
        # Calling .login() on the mock_saml_auth_object returns the self.mock_idp_url.
        mock_saml_auth_object.login.return_value = self.mock_idp_url
        # Calling .get_errors() on the mock_saml_auth_object returns an empty array.
        mock_saml_auth_object.get_errors.return_value = []
        # Calling .is_authenticated() on the mock_saml_auth_object returns True.
        mock_saml_auth_object.is_authenticated.return_value = True
        self.mock_get_saml_auth.return_value = mock_saml_auth_object

        # Mock the get_saml_settings() function in basechem.sso.views.py
        self.patcher_get_saml_settings = mock.patch(
            "basechem.sso.utils.get_saml_settings"
        )
        self.mock_get_saml_settings = self.patcher_get_saml_settings.start()
        self.mock_get_saml_settings.return_value = {
            "idp": {
                "entityId": "http://www.okta.com/entityIdHere",
                "singleSignOnService": {
                    "url": "https://dev-1.okta.com/app/appidhere/something/sso/saml",
                    "binding": "urn:oasis:names:tc:SAML:2.0:bindings:HTTP-Redirect",
                },
                "x509cert": "thex509certgoeshere",
            },
            "sp": {
                "entityId": "http://localhost:8000/saml_login/",
                "assertionConsumerService": {
                    "url": "http://localhost:8000/saml_login/",
                    "binding": "urn:oasis:names:tc:SAML:2.0:bindings:HTTP-POST",
                },
            },
        }

        # Clean up our patchers
        self.addCleanup(self.patch_saml_auth_object.stop)
        self.addCleanup(self.patcher_get_saml_settings.stop)
