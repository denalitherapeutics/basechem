from django.conf import settings
from onelogin.saml2.auth import OneLogin_Saml2_Auth
from onelogin.saml2.idp_metadata_parser import OneLogin_Saml2_IdPMetadataParser


def get_saml_settings(http_host, protocol):
    """Return the SAML settings for the Identity Provider and the Service Provider."""
    # The URL in basechem that receives the Identity Provider (Okta) response.
    acs_url = f"{protocol}://{http_host}/saml_login/"
    # Get the Identity Provider (Okta) metadata.
    idp_data = OneLogin_Saml2_IdPMetadataParser.parse_remote(
        settings.SAML_IP_METADATA_URL
    )
    return {
        "idp": idp_data["idp"],
        "sp": {
            "entityId": acs_url,
            "assertionConsumerService": {
                "url": acs_url,
                "binding": "urn:oasis:names:tc:SAML:2.0:bindings:HTTP-POST",
            },
        },
    }


def get_saml_auth(req):
    """Get and return a SAML auth object."""
    protocol = "https" if req["https"] == "on" else "http"
    saml_settings = get_saml_settings(http_host=req["http_host"], protocol=protocol)
    auth = OneLogin_Saml2_Auth(req, saml_settings)
    return auth


def get_request_data(request):
    """Return a dictionary of data about the request."""
    # Only include the server port in the request data if the HTTP_HOST includes 'localhost'.
    server_port = None
    if "localhost" in request.META.get("HTTP_HOST", ""):
        if "8000" in request.META.get("HTTP_HOST", ""):
            server_port = request.META.get("SERVER_PORT")
        else:
            server_port = "1337"
            request.META["HTTP_HOST"] = "localhost:1337"

    return {
        "https": "on" if request.is_secure() else "off",
        "http_host": request.META.get("HTTP_HOST", ""),
        "script_name": request.META.get("PATH_INFO"),
        "server_port": server_port,
        "get_data": request.GET.copy(),
        "post_data": request.POST.copy(),
    }
