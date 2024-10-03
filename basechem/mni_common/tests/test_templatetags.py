from django.template import Context, Template
from django.test import TestCase, override_settings


class TemplatetagTestCase(TestCase):
    @override_settings(FAKE_SETTING="Hello World")
    def test_get_setting(self):
        """
        Test that `get_setting` retrieves the value of the given setting when possible
        """
        with self.subTest("Setting exists"):
            out = Template(
                "{% load mni_common_tags %}" "{% get_setting 'FAKE_SETTING' %}"
            ).render(Context({}))
            self.assertEqual(out, "Hello World")

        with self.subTest("Setting does not exist"):
            out = Template(
                "{% load mni_common_tags %}" "{% get_setting 'FAKE_SETTING_2' %}"
            ).render(Context({}))
            self.assertEqual(out, "")

    def test_dict_get(self):
        """
        Test that `dict_get` returns the value of the given key in the given dictionary
        """
        data = {"a": "The letter a", "b": "The letter b"}
        with self.subTest("Key exists"):
            out = Template(
                "{% load mni_common_tags %}" "{{ data|dict_get:'a' }}"
            ).render(Context({"data": data}))
            self.assertEqual(out, "The letter a")

        with self.subTest("Key does not exist"):
            out = Template(
                "{% load mni_common_tags %}" "{{ data|dict_get:'c' }}"
            ).render(Context({"data": data}))
            self.assertEqual(out, "None")

        with self.subTest("Not a dictionary"):
            out = Template("{% load mni_common_tags %}" "{{ ''|dict_get:'c' }}").render(
                Context({"data": data})
            )
            self.assertEqual(out, "None")

    def test_no_break_hyphen(self):
        """
        Test that `no_break_hyphen` replaces hyphens with the no-break hyphen
        """
        with self.subTest("One hyphen"):
            out = Template(
                "{% load mni_common_tags %}" "{{ value|no_break_hyphen }}"
            ).render(Context({"value": "-hi there"}))
            self.assertEqual(out, "&amp;#8209;hi there")

        with self.subTest("Multiple hyphens"):
            out = Template(
                "{% load mni_common_tags %}" "{{ value|no_break_hyphen }}"
            ).render(Context({"value": "-hi-there"}))
            self.assertEqual(out, "&amp;#8209;hi&amp;#8209;there")

        with self.subTest("No hyphens"):
            out = Template(
                "{% load mni_common_tags %}" "{{ value|no_break_hyphen }}"
            ).render(Context({"value": "hi there"}))
            self.assertEqual(out, "hi there")

    def test_call_object_method(self):
        """
        Test that `call_object_method` calls an object's method w/ any number of arguments
        """

        class TestObject:
            def test_func(self, a=1, b=2):
                return a + b

        obj = TestObject()
        with self.subTest("No arguments"):
            out = Template(
                "{% load mni_common_tags %}" "{% call_object_method obj 'test_func' %}"
            ).render(Context({"obj": obj}))
            self.assertEqual(out, "3")

        with self.subTest("One argument"):
            out = Template(
                "{% load mni_common_tags %}"
                "{% call_object_method obj 'test_func' 3 %}"
            ).render(Context({"obj": obj}))
            self.assertEqual(out, "5")

        with self.subTest("Multiple arguments"):
            out = Template(
                "{% load mni_common_tags %}"
                "{% call_object_method obj 'test_func' 3 7 %}"
            ).render(Context({"obj": obj}))
            self.assertEqual(out, "10")

    def test_filename(self):
        """
        Test that `filename` returns the filename when given a path
        """
        with self.subTest("Filename only"):
            out = Template(
                "{% load mni_common_tags %}" "{{ filepath|filename }}"
            ).render(Context({"filepath": "abc.pdf"}))
            self.assertEqual(out, "abc.pdf")

        with self.subTest("Filepath"):
            out = Template(
                "{% load mni_common_tags %}" "{{ filepath|filename }}"
            ).render(Context({"filepath": "x/y/z/abc.pdf"}))
            self.assertEqual(out, "abc.pdf")

        with self.subTest("Not a filename"):
            out = Template(
                "{% load mni_common_tags %}" "{{ filepath|filename }}"
            ).render(Context({"filepath": "hello world"}))
            self.assertEqual(out, "hello world")

    def test_bool_as_int(self):
        """
        Test that `bool_as_int` converts booleans to integers
        """
        with self.subTest("True"):
            out = Template(
                "{% load mni_common_tags %}" "{{ bool|bool_as_int }}"
            ).render(Context({"bool": True}))
            self.assertEqual(out, "1")

        with self.subTest("False"):
            out = Template(
                "{% load mni_common_tags %}" "{{ bool|bool_as_int }}"
            ).render(Context({"bool": False}))
            self.assertEqual(out, "0")

    def test_sorted_dict(self):
        """
        Test that `sorted_dict` returns the dictionary in order by key
        """
        data = {"b": "B", "c": "C", "a": "A"}
        out = Template(
            "{% load mni_common_tags %}"
            "{% for k in data|sorted_dict %}{{ k }}{% endfor %}"
        ).render(Context({"data": data}))
        self.assertEqual(out, "abc")

    def test_get_attribute(self):
        """
        Test that `get_attribute` retrieves the value of an object's attribute when possible
        """

        class TestObject:
            a = None
            b = None

            def __init__(self, a=None, b=None):
                self.a = a
                self.b = b

        obj_b = TestObject(a="A", b={"hello": "world"})
        obj_a = TestObject(a="Z", b=obj_b)

        with self.subTest("Simple attribute"):
            out = Template(
                "{% load mni_common_tags %}" "{{ obj|get_attribute:'a' }}"
            ).render(Context({"obj": obj_b}))
            self.assertEqual(out, "A")

        with self.subTest("Nested attribute"):
            out = Template(
                "{% load mni_common_tags %}" "{{ obj|get_attribute:'b__a' }}"
            ).render(Context({"obj": obj_a}))
            self.assertEqual(out, "A")

        with self.subTest("Nested dictionary"):
            out = Template(
                "{% load mni_common_tags %}" "{{ obj|get_attribute:'b__b__hello' }}"
            ).render(Context({"obj": obj_a}))
            self.assertEqual(out, "world")

        with self.subTest("Attribute does not exist"):
            out = Template(
                "{% load mni_common_tags %}" "{{ obj|get_attribute:'c' }}"
            ).render(Context({"obj": obj_a}))
            self.assertEqual(out, "")

    def test_class_name(self):
        """
        Test that `class_name` returns the class name of the object
        """
        with self.subTest("Nonetype"):
            out = Template("{% load mni_common_tags %}" "{{ obj|class_name }}").render(
                Context({"obj": None})
            )
            self.assertEqual(out, "NoneType")

        with self.subTest("Built-in types"):
            # String
            out = Template("{% load mni_common_tags %}" "{{ obj|class_name }}").render(
                Context({"obj": "hello"})
            )
            self.assertEqual(out, "str")
            # Int
            out = Template("{% load mni_common_tags %}" "{{ obj|class_name }}").render(
                Context({"obj": 1})
            )
            self.assertEqual(out, "int")

        with self.subTest("Custom types"):

            class ACustomClass(object):
                def __init__(self, *args, **kwargs):
                    pass

            out = Template("{% load mni_common_tags %}" "{{ obj|class_name }}").render(
                Context({"obj": ACustomClass()})
            )
            self.assertEqual(out, "ACustomClass")

    def test_is_instance(self):
        """
        Test that the `is_instance` tag returns expected boolean
        """
        text = "No MMP found."
        out = Template(
            "{% load mni_common_tags %}" '{{ text|is_instance:"str" }}'
        ).render(Context({"text": text}))
        self.assertEqual("True", out)
        out = Template(
            "{% load mni_common_tags %}" '{{ text|is_instance:"int" }}'
        ).render(Context({"text": text}))
        self.assertEqual("False", out)

    def test_get_range(self):
        """
        Test that `get_range` returns a range of numbers
        """
        with self.subTest("Range from 0 -> 0"):
            out = Template(
                "{% load mni_common_tags %}"
                "{% for i in 0|get_range %}{{ i }}{% endfor %}"
            ).render(Context({}))
            self.assertEqual(out, "")

        with self.subTest("Range from 0 -> 5"):
            out = Template(
                "{% load mni_common_tags %}"
                "{% for i in 5|get_range %}{{ i }}{% endfor %}"
            ).render(Context({}))
            self.assertEqual(out, "01234")

        with self.subTest("Range from 0 -> 5 (string)"):
            out = Template(
                "{% load mni_common_tags %}"
                "{% for i in '5'|get_range %}{{ i }}{% endfor %}"
            ).render(Context({}))
            self.assertEqual(out, "01234")
