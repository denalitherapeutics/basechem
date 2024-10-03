from urllib.parse import urlencode

from django.conf import settings
from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from basechem.sso.tests.base import BasechemSSOMixin


class TestLoginPage(BasechemSSOMixin, TestCase):
    """Test case for the login view."""

    def test_reverse_of_login(self):
        self.assertEqual("/", reverse("login"))

    def test_get_login_page_unauthenticated(self):
        """Login page should load"""
        self.client.logout()
        response = self.client.get(reverse("login"))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "sso/login.html")

    def test_authenticated(self):
        self.client.login(username=self.user1.username, password="testpassword")
        response = self.client.get(reverse("login"))
        self.assertRedirects(response, reverse(settings.LOGIN_REDIRECT_URL))

    def test_authenticated_with_next(self):
        # If already authenticated, and ?next provided, we redirect there.
        next_url = reverse("homepage") + "?q=einstein"
        response = self.client.get(
            reverse("login") + "?" + urlencode({"next": next_url})
        )
        self.assertRedirects(response, next_url)


class TestSSOLoginView(BasechemSSOMixin, TestCase):
    """Test case for the sachml_login view."""

    def test_get_unauthenticated(self):
        """An unauthenticated user is redirected to Okta to log in."""
        self.client.logout()
        response = self.client.get(reverse("sso_login"))
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response["Location"], self.mock_idp_url)

    def test_get_authenticated(self):
        """An authenticated user is redirected to the 'sso_success' view."""
        self.client.login(username=self.user1.username, password="testpassword")
        response = self.client.get(reverse("sso_login"))
        self.assertEqual(response.status_code, 302)
        self.assertEqual(response["Location"], reverse("sso_success"))

    def test_post(self):
        """POSTing to the 'sso_login' view is not allowed."""
        response = self.client.post(reverse("sso_login"))
        self.assertEqual(response.status_code, 405)


class TestSSOLoginSuccessView(BasechemSSOMixin, TestCase):
    """Test case for the saml_success view."""

    def test_unauthenticated(self):
        self.client.logout()

        with self.subTest("GET"):
            response = self.client.get(reverse("sso_success"))
            self.assertEqual(response.status_code, 302)
            self.assertEqual(response["Location"], reverse(settings.LOGIN_URL))

        with self.subTest("POST"):
            # Calling .get_errors() on the auth object returns errors.
            self.mock_get_saml_auth.return_value.get_errors.return_value = ["errors"]
            self.mock_get_saml_auth.return_value.get_last_error_reason.return_value = (
                "something was invalid"
            )

            response = self.client.post(reverse("sso_success"))

            self.assertEqual(response.status_code, 200)
            # The errors from the POST are shown to the user.
            self.assertEqual(
                response.context["errors"],
                [
                    self.mock_get_saml_auth.return_value.get_last_error_reason.return_value
                ],
            )

    def test_get(self):
        """GETting the saml_success view redirects the user to the homepage."""
        # We assume the  user is logged in (unauthenticated users are tested in
        # test_unauthenticated()).
        self.client.login(username=self.user1.username, password="testpassword")

        response = self.client.get(reverse("sso_success"))
        self.assertRedirects(response, reverse(settings.LOGIN_REDIRECT_URL))

    def test_post_sso_errors(self):
        """SSO login errors are shown to the user in the template."""
        # Calling .get_errors() on the auth object returns errors.
        self.mock_get_saml_auth.return_value.get_errors.return_value = ["errors"]
        self.mock_get_saml_auth.return_value.get_last_error_reason.return_value = (
            "something was invalid"
        )

        response = self.client.post(reverse("sso_success"))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            response.context["errors"],
            [self.mock_get_saml_auth.return_value.get_last_error_reason.return_value],
        )

    def test_post_sso_not_authenticated(self):
        """If SSO login says the user is not authenticated, errors are shown to the user."""
        # Calling .is_authenticated() on the auth object returns False.
        self.mock_get_saml_auth.return_value.is_authenticated.return_value = False

        response = self.client.post(reverse("sso_success"))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["errors"], ["Authentication error"])

    def test_post_sso_success_local_user_found(self):
        """If SSO login succeeds, and finds a user in basechem, local user is logged in."""
        # Calling .get_attributes() on the auth object returns the self.user's email.
        self.mock_get_saml_auth.return_value.get_attributes.return_value = {
            "email": [self.user1.email],
            "username": [self.user1.email],
        }

        response = self.client.post(reverse("sso_success"))

        self.assertRedirects(response, reverse(settings.LOGIN_REDIRECT_URL))
        self.assertEqual(response.wsgi_request.user, self.user1)

    def test_post_sso_success_local_user_found_with_next(self):
        """If SSO login succeeds, and finds a user in basechem, local user is logged in."""
        # If we provided a "next" url, it gets passed to success as RelayState and we
        # redirect to it.
        self.mock_get_saml_auth.return_value.get_attributes.return_value = {
            "email": [self.user1.email],
            "username": [self.user1.email],
        }
        next_url = "/foo/bar"
        response = self.client.post(
            reverse("sso_success"), data={"RelayState": next_url}
        )
        self.assertRedirects(response, next_url, fetch_redirect_response=False)

    def test_post_sso_success_no_local_user(self):
        """
        Handle the case that user is logged in with SSO, but no matching user exists in basechem.

        If the SSO login succeeds, but no matching user exists in basechem, the
        user is created, then logged in.
        """
        # Calling .get_attributes() on the auth object returns data that does not
        # match any user in basechem.
        new_user_email = "someemail@example.com"
        self.mock_get_saml_auth.return_value.get_attributes.return_value = {
            "username": [new_user_email],
            "email": [new_user_email],
            "first_name": ["New"],
            "last_name": ["User"],
        }
        self.assertFalse(get_user_model().objects.filter(email=new_user_email).exists())

        response = self.client.post(reverse("sso_success"))

        self.assertRedirects(response, reverse(settings.LOGIN_REDIRECT_URL))
        self.assertEqual(response.wsgi_request.user.email, new_user_email)
        self.assertTrue(
            get_user_model()
            .objects.filter(email=new_user_email, first_name="New", last_name="User")
            .exists()
        )

    def test_post_sso_success_multiple_users(self):
        """
        Handle case that user is logged in with SSO, but multiple matching users exist in basechem.

        If the SSO login succeeds, but multiple matching user exist in basechem,
        the user receives an error.
        """
        # Calling .get_attributes() on the auth object returns data that matches
        # two users in basechem.
        user1 = get_user_model().objects.create(
            username="user1@example.com",
            email="user1@example.com",
            first_name="User",
            last_name="One",
        )
        get_user_model().objects.create(
            username="USER1@example.com",
            email="USER1@example.com",
            first_name="User",
            last_name="One",
        )
        self.mock_get_saml_auth.return_value.get_attributes.return_value = {
            "username": [user1.email],
            "email": [user1.email],
            "first_name": [user1.first_name],
            "last_name": [user1.last_name],
        }

        response = self.client.post(reverse("sso_success"))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["errors"], ["multiple users found!"])
