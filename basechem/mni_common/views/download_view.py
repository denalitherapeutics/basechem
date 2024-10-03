import os
from wsgiref.util import FileWrapper

from django.contrib.auth.mixins import LoginRequiredMixin
from django.http import HttpResponse
from django.views.generic import View

from ..storage import select_media_storage, select_static_storage


class DownloadFileView(LoginRequiredMixin, View):
    """
    View to download static or media files
    """

    @staticmethod
    def get(request, **kwargs):
        storage = select_static_storage()
        if kwargs.get("storage") == "media":
            storage = select_media_storage()
        filepath = kwargs.get("filepath")
        wrapper = FileWrapper(storage.open(filepath, "rb"))
        response = HttpResponse(wrapper, content_type="application/force-download")
        response[
            "Content-Disposition"
        ] = f"attachment; filename={os.path.basename(filepath)}"
        return response
