from django import forms
from django.core.exceptions import ValidationError

from .widgets import MultipleFileInput


class MultipleFileField(forms.FileField):
    """
    FileField for uploading multiple files in a single field. This code is taken directly from
    the Django Documentation. The docs imply that this will be in the Django codebase in a future release (>4.2)
    (https://docs.djangoproject.com/en/4.2/topics/http/file-uploads/#uploading-multiple-files).
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("widget", MultipleFileInput())
        super().__init__(*args, **kwargs)

    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            result = [single_file_clean(d, initial) for d in data]
        else:
            result = single_file_clean(data, initial)
        return result


class NoValidateMultipleChoiceField(forms.MultipleChoiceField):
    """
    A custom multiple choice field that does not check that the selected values are included
    in the original choices. This is used in cases where users may add their own values
    in addition to the pre-set choice list when filling out the form.
    """

    def validate(self, value):
        if self.required and not value:
            raise ValidationError("This is a required field.")


class NoValidateChoiceField(forms.ChoiceField):
    """
    A custom choice field that does not check that the selected values are included
    in the original choices. This is used in cases when the choices are added dynamically
    as the user fills out the form
    """

    def validate(self, value):
        return
