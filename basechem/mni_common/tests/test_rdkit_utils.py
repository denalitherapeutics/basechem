from django.test import TestCase
from rdkit import Chem

from ..rdkit_utils import moltext_to_svg


class RdkitUtilsTestCase(TestCase):
    def test_moltext_to_svg(self):
        """
        Test `moltext_to_svg` converts a string of moltext to an svg string
        """
        with self.subTest("Invalid moltext"):
            result = moltext_to_svg("this is not moltext")
            self.assertEqual(result, "")

        with self.subTest("Valid moltext"):
            moltext = Chem.MolToMolBlock(Chem.MolFromSmiles("C1CCCCC1"))
            result = moltext_to_svg(moltext)
            self.assertEqual(result[:5], "<?xml")
            self.assertEqual(result[-7:], "</svg>\n")
