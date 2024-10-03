import logging
from urllib.parse import urlparse

from django.contrib import messages
from django.shortcuts import redirect
from django.urls import resolve


def show_2d_warning(request, collection, analysis):
    """
    Add message to request if 2D, not 3D, structures exist in collection compounds
    :param request: a request object
    :param collection: a Collection object
    :param analysis: a string, the analysis this is called from (ex: "torsion")
    """
    collection.get_cos_for_analysis(analysis)
    num_2d = sum(
        [co.molblock == "" for co in collection.get_cos_for_analysis(analysis)]
    )
    if num_2d:
        verb = "is" if num_2d == 1 else "are"
        messages.info(
            request,
            f"{num_2d} of your uploaded structures {verb} 2D. If you meant to do this, great! Basechem will generate a 3D starting conformation for you. If not, please create a new collection with your intended structures.",
        )


class ExceptionMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        return self.get_response(request)

    @staticmethod
    def process_exception(request, exception):
        """
        Redirects user to the previous page they were on (HTTP REFERER). Also, calls an optional handle_exception
        function from the destination view class that raised the exception to display a message for the user. Uses
        the django messages framework (https://docs.djangoproject.com/en/3.1/ref/contrib/messages/) to show errors to
        user.
        :param request: HTTPRequest object for the request made
        :param exception: Exception object that was caught by the middleware
        """
        http_referer = urlparse(request.META["HTTP_REFERER"])

        # get the view class of the URL that raised the exception
        destination_class = resolve(request.path).func.view_class
        handle_exception = getattr(destination_class, "handle_exception", None)
        error_msg = ""

        # check to make sure handle exception is implemented for the view
        if callable(handle_exception):
            error_msg = handle_exception(exception)
        if not error_msg:
            error_msg = "Oh no! An unexpected error has occurred. We've notified our administrators and they are on the case!"

        # Log message to sentry
        logging_message = f"{request.user} encountered an unexpected error. {repr(type(exception))} raised in {repr(destination_class)} view. Error:\n{exception}"
        logging.error(logging_message)

        messages.error(request, error_msg)  # show message to user
        return redirect(http_referer.path)
