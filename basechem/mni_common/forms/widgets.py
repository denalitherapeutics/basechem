from django import forms


class DateWidget(forms.DateInput):
    input_type = "date"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.attrs = {"style": "width:200px", "type": "date"}


class MultipleFileInput(forms.ClearableFileInput):
    """
    ClearableFileInput that allows multiple files to be attached
    """

    allow_multiple_selected = True
