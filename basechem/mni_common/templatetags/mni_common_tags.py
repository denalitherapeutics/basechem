import os

from django import template
from django.conf import settings

register = template.Library()


@register.simple_tag
def get_setting(setting_name):
    """
    Returns the current value of a Django setting
    :param setting_name: the name of an attribute in the Django settings (ex: 'ENVIRONMENT')
    """
    return getattr(settings, setting_name, "")


@register.filter(name="dict_get")
def dict_get(dictionary, key):
    """
    Returns the value of `key` in `dictionary`
    :param dictionary: a dictionary
    :param key: a string, the key whose value should be retrieved
    """
    if isinstance(dictionary, dict):
        return dictionary.get(key)


@register.filter(name="no_break_hyphen")
def no_break_hyphen(value):
    """
    When a single word (string w/ no spaces) is longer than the container it's given, the rendered
    html will split the word at hyphens if it can. Hyphens appear in ASO sequences, which leads to
    a sequence like this:
    "-ddsfhsdfjkshdfkjshdfkjfs"
    being rendered like this:
    "-
    dfjsdflksdflkjsdflkjsd"
    . This function replaces all the hyphens with no-breaking hyphens (&#8209;).
    """
    return value.replace("-", "&#8209;")


@register.simple_tag
def call_object_method(obj, method_name, *args):
    """
    This template tag can be used to call a method that requires arguments on an object.
    For example, to call `protein.can_user_cancel(user)`, write:
    {% call_object_method protein can_user_cancel user %}
    :param obj: an instance of a Python class
    :param method_name: a string, the name of a method on that class
    :param *args: the arguments that method should be called with
    :returns: the result of the method called on the object with the given arguments
    """
    method = getattr(obj, method_name)
    return method(*args)


@register.filter(name="filename")
def filename(filepath):
    """
    :param filepath: a string, the filepath for a file
    :returns: a string, just the filename
    """
    return os.path.basename(filepath)


@register.filter(name="bool_as_int")
def bool_as_int(value):
    """
    This function is useful when using a python boolean value in javascript. Without this,
    the values `True` and `False` become the strings "True" and "False"
    :param value: a boolean value
    :returns: 1 if true, 0 if false
    """
    if value:
        return 1
    return 0


@register.filter(name="sorted_dict")
def sorted_dict(data):
    """
    :param data: a dictionary
    :returns: the dictionary in alphabetical order by key
    """
    sorted_data = dict(sorted(data.items(), key=lambda x: x[0].lower()))
    return sorted_data


@register.filter(name="get_attribute")
def get_attribute(obj, attribute):
    """
    Returns the obj's value for an attribute. If the attribute has a double underscore (meaning
    it is a nested field, as with a ForeignKeyField or JSONField), then this function calls
    itself recursively to get the nested value.
    :param obj: a python object (can be a Django object or not)
    :param attribute: a string, the attribute whose value should be retrieved
    :returns: the obj's value for `attribute`, or the empty string if it does not exist.
    """
    empty_value = ""
    # Try to return the obj's attribute, if possible.
    if hasattr(obj, str(attribute)):
        return getattr(obj, attribute)
    elif hasattr(obj, "keys") and attribute in obj.keys():
        # If the attribute is a key for the object (meaning that the object is a
        # dictionary), then get that key.
        return obj.get(attribute, empty_value)
    elif "__" in attribute:
        # The attribute has a nested field (like if the field is a JSONField), so
        # return the value of the nested field by calling getattribute() on it.
        field_name, nested_field_name = attribute.split("__", 1)

        return get_attribute(get_attribute(obj, field_name), nested_field_name)

    return empty_value


@register.filter(name="class_name")
def class_name(obj):
    """
    :param obj: a python object
    :returns: a string, the class name of the object
    """
    try:
        return obj.__class__.__name__
    except:
        return str(type(obj))


@register.filter(name="is_instance")
def is_instance(val, instance_type):
    """
    Given a value and an instance type, returns True if the value is of type instance_type
    :param val: variable who's instance type needs to be determined
    :return: True if value is of type instance_type, False otherwise
    """
    return isinstance(val, eval(instance_type))


@register.filter(name="get_range")
def get_range(value):
    """
    Given a number, return a list of integers from 0 to that number (exclusive)
    :param value: the top of the range
    """
    return range(int(value))
