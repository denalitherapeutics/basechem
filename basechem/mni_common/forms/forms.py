from crispy_forms.helper import FormHelper
from django import forms

from .fields import MultipleFileField


class CrispyModelForm(forms.ModelForm):
    """
    A ModelForm that has a FormHelper already attached
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.form_tag = False
        self.helper.label_class = "form-label"


class CrispyForm(forms.Form):
    """
    A Form that has a FormHelper already attached
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.form_tag = False
        self.helper.label_class = "form-label"


class BugReportForm(CrispyForm):
    description = forms.CharField(
        max_length=1000,
        required=True,
        help_text="Describe the bug/suggestion",
        widget=forms.Textarea(attrs={"style": "height: 100px", "required": False}),
    )

    url = forms.CharField(max_length=200, widget=forms.HiddenInput())
    files = MultipleFileField(
        help_text="Please attach any relevant files or screenshots", required=False
    )

    def __init__(self, url=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if url:
            self.initial["url"] = url
