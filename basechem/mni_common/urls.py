from django.urls import path

from .views.bug_report_view import SendBugReportView
from .views.download_view import DownloadFileView

urlpatterns = [
    path(
        r"download/<str:storage>/<path:filepath>",
        DownloadFileView.as_view(),
        name="download",
    ),
    path(r"send-bug-report/", SendBugReportView.as_view(), name="send_bug_report"),
]
