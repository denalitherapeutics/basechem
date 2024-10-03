import random
from unittest import mock

from django.test import TestCase
from faker import Faker

from basechem.sso.utils import get_request_data, get_saml_auth, get_saml_settings


class TestGetRequestData(TestCase):
    """Test case for the get_request_data() utility function."""

    def test_get_data(self):
        """Given a request object, the function returns a dictionary of data."""
        mock_request_object = mock.MagicMock()
        mock_request_object.META = {
            "HTTP_HOST": "localhost:8000",
            "PATH_INFO": "/saml_success/",
            "SERVER_PORT": None,
        }
        mock_request_object.GET = {}
        mock_request_object.POST = {}

        subtests = (
            # is_request_secure, http_host, description
            (True, "localhost:8000", "secure request for localhost"),
            (True, "somesite.com:8000", "secure request for somesite.com"),
            (False, "localhost:8000", "insecure request for localhost"),
            (False, "somesite.com:8000", "insecure request for somesite.com"),
        )

        for is_request_secure, http_host, description in subtests:
            with self.subTest(description):
                mock_request_object.is_secure.return_value = is_request_secure
                mock_request_object.META["HTTP_HOST"] = http_host

                # Create the expected_result
                expected_result = {
                    "https": "on" if is_request_secure else "off",
                    "http_host": mock_request_object.META["HTTP_HOST"],
                    "script_name": mock_request_object.META["PATH_INFO"],
                    "get_data": mock_request_object.GET,
                    "post_data": mock_request_object.POST,
                }
                # If the http_host includes 'localhost', then the server port
                # is included in the request data.
                if "localhost" in http_host:
                    expected_result["server_port"] = mock_request_object.META[
                        "SERVER_PORT"
                    ]
                else:
                    expected_result["server_port"] = None

                self.assertEqual(get_request_data(mock_request_object), expected_result)


@mock.patch("basechem.sso.utils.get_saml_settings")
@mock.patch("basechem.sso.utils.OneLogin_Saml2_Auth")
class TestGetSamlAuth(TestCase):
    """Test case for the get_saml_auth() utility function."""

    def test_get_saml_auth(self, mock_saml2_auth_obj, mock_get_saml_settings):
        """Given a dictionary of request data, function returns an OneLogin_Saml2_Auth object."""

        request_data = {
            "https": "on",
            "http_host": "localhost:8000",
            "script_name": "/saml_success/",
            "server_port": "8000",
            "get_data": {},
            "post_data": {},
        }

        result = get_saml_auth(request_data)

        self.assertEqual(mock_saml2_auth_obj.call_count, 1)
        mock_saml2_auth_obj.assert_called_with(request_data, mock_get_saml_settings())
        self.assertEqual(result, mock_saml2_auth_obj())


@mock.patch("basechem.sso.utils.OneLogin_Saml2_IdPMetadataParser")
class TestGetSamlSettings(TestCase):
    """Test case for the get_saml_settings() utility function."""

    def test_get_saml_settings(self, mock_saml2_metadata_parser):
        """Given a dictionary of request data, function returns an OneLogin_Saml2_Auth object."""

        mock_saml2_metadata_parser.parse_remote.return_value = {
            "idp": {
                "entityId": "http://www.okta.com/entityIdHere",
                "singleSignOnService": {
                    "url": "https://dev-1.okta.com/app/appidhere/something/sso/saml",
                    "binding": "urn:oasis:names:tc:SAML:2.0:bindings:HTTP-Redirect",
                },
                "x509cert": "thex509certgoeshere",
            },
            "sp": {
                "NameIDFormat": "urn:oasis:names:tc:SAML:1.1:nameid-format:unspecified"
            },
        }

        http_host = Faker().domain_name
        protocol = random.choice(["http", "https"])

        expected_result = {
            "idp": mock_saml2_metadata_parser.parse_remote.return_value["idp"],
            "sp": {
                "entityId": f"{protocol}://{http_host}/saml_login/",
                "assertionConsumerService": {
                    "url": f"{protocol}://{http_host}/saml_login/",
                    "binding": "urn:oasis:names:tc:SAML:2.0:bindings:HTTP-POST",
                },
            },
        }

        self.assertEqual(get_saml_settings(http_host, protocol), expected_result)
