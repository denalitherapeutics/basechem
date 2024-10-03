import mimetypes

from django.conf import settings
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.files.uploadedfile import InMemoryUploadedFile, TemporaryUploadedFile
from django.core.mail import EmailMultiAlternatives, mail_admins
from django.http import JsonResponse
from django.urls import reverse_lazy
from django.views.generic import FormView

from ..forms.forms import BugReportForm


class SendBugReportView(LoginRequiredMixin, FormView):
    form_class = BugReportForm
    # The success URL never gets used (because of the JSON Responses), but is a requirement for FormViews
    success_url = reverse_lazy("login")

    def form_valid(self, form):
        super().form_valid(form)
        if form.is_valid():
            subject = f"Bug reported by {self.request.user}"
            full_url = f"{settings.BASE_URL}{form.cleaned_data['url']}"
            message = f"{self.request.user} reported a bug at {full_url} with the following description:\n{form.cleaned_data['description']}"
            files = []
            if form.cleaned_data.get("files"):
                files = self.request.FILES.getlist("files")

            try:
                # Attach files to email and send to admins
                recipients = [email for admin, email in settings.ADMINS]
                msg = EmailMultiAlternatives(
                    subject, message, settings.DEFAULT_FROM_EMAIL, recipients
                )
                for file in files:
                    # Django keeps files in memory if they're small enough, so we have to accommodate
                    # both InMemoryUploadedFile and TemporaryUploadedFile:
                    if isinstance(file, TemporaryUploadedFile):
                        msg.attach_file(file.temporary_file_path())
                    elif isinstance(file, InMemoryUploadedFile):
                        msg.attach(
                            file.name,
                            file.file.getvalue(),
                            mimetypes.guess_type(file.name)[0],
                        )
                msg.send()
            except Exception as e:
                # Catch this exception in case there are unexpected errors due to the attachments
                mail_admins(
                    "Bug report failed to send",
                    f'{self.request.user} tried to file a bug report, but PropheSeq encountered an "{e}" error.',
                )
                return JsonResponse(
                    {
                        "errors": {
                            "__all__": [
                                f"Something went wrong - our development team will not receive this report. Please email {recipients[0]} with your bug."
                            ]
                        }
                    }
                )

        return JsonResponse({"errors": form.errors})

    def form_invalid(self, form):
        super().form_invalid(form)
        return JsonResponse({"errors": form.errors})
