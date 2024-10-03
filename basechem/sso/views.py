import datetime

from django.conf import settings
from django.contrib.auth import authenticate, get_user_model, login
from django.contrib.auth.forms import AuthenticationForm
from django.http import HttpResponseRedirect
from django.shortcuts import render, reverse
from django.utils.http import url_has_allowed_host_and_scheme
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods

from basechem.common.analytic import LOGIN, Analytic
from basechem.sso.utils import get_request_data, get_saml_auth


@require_http_methods(["GET", "POST"])
def homepage_gateway(request):
    """
    Unauthenticated users are directed to the login page while authenticated users are directed to the homepage.
    This view also handles all authentication but for SSO it redirects to another view
    """
    if request.user.is_authenticated:
        # Log this login attempt to AWS
        day_of_week = datetime.datetime.now().strftime("%A")
        Analytic(LOGIN, None, request.user, day_of_week=day_of_week)
        if "next" in request.GET:
            next_url = request.GET["next"]
            if url_has_allowed_host_and_scheme(next_url, allowed_hosts=None):
                return HttpResponseRedirect(next_url)
        return HttpResponseRedirect(reverse(settings.LOGIN_REDIRECT_URL))
    else:
        template = "sso/login.html"
        sso_available = settings.SAML_IP_METADATA_URL
        form = AuthenticationForm()

        if "next" in request.GET:
            # Save "next" in session for use in later request
            request.session["next"] = request.GET["next"]

        if request.method == "POST":
            # Used when logging in with basic authentication instead of SSO
            form = AuthenticationForm(data=request.POST)
            if form.is_valid():
                username = form.cleaned_data.get("username")
                password = form.cleaned_data.get("password")
                user = authenticate(username=username, password=password)
                login(request, user)
                day_of_week = datetime.datetime.now().strftime("%A")
                Analytic(LOGIN, None, request.user, day_of_week=day_of_week)
                if "next" in request.session:
                    next_url = request.session["next"]
                    if url_has_allowed_host_and_scheme(next_url, allowed_hosts=None):
                        return HttpResponseRedirect(next_url)
                return HttpResponseRedirect(reverse(settings.LOGIN_REDIRECT_URL))
            else:
                return render(
                    request,
                    "sso/login.html",
                    context={
                        "sso_available": sso_available,
                        "form": form,
                        "errors": form.errors["__all__"],
                    },
                )
        return render(
            request,
            template,
            context={"sso_available": sso_available, "form": form, "errors": []},
        )


@require_http_methods(["GET"])
def saml_login(request):
    """Take a GET request, and redirect user to login with the Identity Provider (Okta)."""
    # If the request's user is already logged in, then redirect to the 'saml_success' view.
    if request.user.is_authenticated:
        return HttpResponseRedirect(reverse("sso_success"))

    request_dict = get_request_data(request)
    auth = get_saml_auth(request_dict)
    # If we had a "next" parameter, tell SAML to go there after login.
    idp_login_url = auth.login(return_to=request.session.get("next", None))
    return HttpResponseRedirect(idp_login_url)


@csrf_exempt
@require_http_methods(["GET", "POST"])
def saml_success(request):
    """
    Handle result of SSO authentication.

    If the user has successfully authenticated through SSO, then log the user in.
    Otherwise, show the user an error.
    """
    # Successful SSO authentication results in a POST request to this view.
    # We check things out, then if all looks good, redirect somewhere.
    login_redirect = reverse(settings.LOGIN_REDIRECT_URL)
    day_of_week = datetime.datetime.now().strftime("%A")
    if request.method == "POST":
        request_dict = get_request_data(request)

        auth = get_saml_auth(request_dict)
        auth.process_response()
        errors = auth.get_errors()
        if not errors:
            if auth.is_authenticated():
                # Authenticate the user to basechem using the 'username' field.
                matching_users = get_user_model().objects.filter(
                    username__iexact=auth.get_attributes()["username"][0]
                )
                errors = []
                if not matching_users:
                    # There is no matching user in basechem for this username,
                    # so create one, and log that user in.
                    new_user = get_user_model().objects.create(
                        username=auth.get_attributes()["username"][0],
                        email=auth.get_attributes()["email"][0],
                        first_name=auth.get_attributes()["first_name"][0],
                        last_name=auth.get_attributes()["last_name"][0],
                    )
                    login(request, new_user)
                    # Log the successful login to AWS
                    Analytic(LOGIN, None, request.user, day_of_week=day_of_week)
                elif matching_users.count() > 1:
                    # Multiple users in basechem have this username. This shouldn't
                    # happen, so give the user an error.
                    errors.append("multiple users found!")
                    return render(
                        request, "sso/login.html", {"errors": ["multiple users found!"]}
                    )
                else:
                    # Get the matching user in basechem, and log this user in.
                    the_user = matching_users.first()
                    login(request, the_user)
                    # Log the successful login to AWS
                    Analytic(LOGIN, None, request.user, day_of_week=day_of_week)

                if request_dict["post_data"].get("RelayState"):
                    login_redirect = request_dict["post_data"]["RelayState"]

                return HttpResponseRedirect(login_redirect)
            else:
                # The user was not authenticated, so give the user an error.
                return render(
                    request, "sso/login.html", {"errors": ["Authentication error"]}
                )
        else:
            # There were errors with authentication, so show them to the user.
            human_readable_error = auth.get_last_error_reason()
            return render(request, "sso/login.html", {"errors": [human_readable_error]})

    else:
        # Non-POST requests are redirected to the LOGIN_REDIRECT_URL or LOGIN_URL.
        if request.user.is_authenticated:
            return HttpResponseRedirect(login_redirect)
        return HttpResponseRedirect(reverse(settings.LOGIN_URL))
