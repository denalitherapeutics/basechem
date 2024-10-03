from django.test import TestCase

from ..str_utils import (
    capitalize_first,
    get_A1_column_index_from_name,
    get_A1_column_name_from_index,
    natsort,
)


class StrUtilsTestCase(TestCase):
    def test_capitalize_first(self):
        """
        Test `capitalize_first` capitalizes the first character of a string without changing the
        case of the rest of the string
        """
        with self.subTest("All lower"):
            self.assertEqual(capitalize_first("hello"), "Hello")

        with self.subTest("All upper"):
            self.assertEqual(capitalize_first("HELLO"), "HELLO")

        with self.subTest("Mixed case"):
            self.assertEqual(capitalize_first("heLlo"), "HeLlo")

    def test_natsort(self):
        """
        Test `natsort` can be used to natural sort a list of strings
        """
        for l, expected_l in [
            (["z", "d", "a", "f"], ["a", "d", "f", "z"]),
            (
                ["11", "20", "21", "1", "2", "9", "10"],
                ["1", "2", "9", "10", "11", "20", "21"],
            ),
            (["a", "b", "c", "A", "B", "C"], ["a", "A", "b", "B", "c", "C"]),
        ]:
            self.assertEqual(sorted(l, key=natsort), expected_l)

    def test_get_A1_column_name_from_index(self):
        """
        Test `get_A1_column_name_from_index` can convert a column index to a column name using the A1-reference style
        """
        for col_idx, expected_col_name in [
            (1, "A"),
            (26, "Z"),
            (27, "AA"),
            (52, "AZ"),
            (53, "BA"),
            (702, "ZZ"),
            (703, "AAA"),
            (18278, "ZZZ"),
            (18279, "AAAA"),
        ]:
            with self.subTest(f"{col_idx} = {expected_col_name}"):
                self.assertEqual(
                    get_A1_column_name_from_index(col_idx), expected_col_name
                )

        with self.subTest("Non-integer"):
            with self.assertRaises(ValueError):
                get_A1_column_name_from_index(1.5)

        with self.subTest("Negative integer"):
            with self.assertRaises(ValueError):
                get_A1_column_name_from_index(-1)

    def test_get_A1_column_index_from_name(self):
        """
        Test `get_A1_column_index_from_name` can convert a column name to a column index using the A1-reference style
        """
        for col_name, expected_col_idx in [
            ("a", 1),
            ("A", 1),
            ("Z", 26),
            ("AA", 27),
            ("AZ", 52),
            ("BA", 53),
            ("ZZ", 702),
            ("AAA", 703),
            ("aaa ", 703),
            ("ZZZ", 18278),
            ("AAAA", 18279),
        ]:
            with self.subTest(f"{col_name} = {expected_col_idx}"):
                self.assertEqual(
                    get_A1_column_index_from_name(col_name), expected_col_idx
                )

        with self.subTest("Non alpha"):
            with self.assertRaises(ValueError):
                get_A1_column_index_from_name("AK6")
