def include_initial_in_choices(form, field_name, initial, multiple=False):
    """
    Used with `NoValidate` fields, this function is called in a form's `init` method to make
    sure that all initial values are available options (to avoid the field being blank when
    there is an initial value). This function adds the given `initial` value to the field's
    choices
    :param form: a Django Form object, the form being initialized
    :param field_name: a string, the name of the field being processed
    :param initial: a string, the initial value
    :param multiple: optional, a boolean, is this a multiple choice field?
    """
    field = form.fields[field_name]
    if multiple:
        initial = initial.split(", ")

    # Add initial value to available choices
    if field.choices:
        current_choices = [value for value, _ in field.choices]
        initial_vals = initial if isinstance(initial, list) else [initial]
        for initial_val in initial_vals:
            if initial_val not in current_choices:
                field.choices += [(initial_val, initial_val)]

    form.initial[field_name] = initial
