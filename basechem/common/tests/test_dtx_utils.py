from django.test import TestCase, tag

from basechem.common.dtx_utils import check_dtx_for_inchi


@tag("local", "dtx", "external")
class DTXUtilTest(TestCase):
    def test_check_dtx_for_inchi(self):
        with self.subTest("single match"):
            dn153 = "InChI=1/C14H13ClF5N5/c15-12-10(6-25(24-12)8(3-16)4-17)22-13-21-5-9(14(18,19)20)11(23-13)7-1-2-7/h5-8H,1-4H2,(H,21,22,23)"
            self.assertEqual("DN0000153", check_dtx_for_inchi(dn153))

        with self.subTest("no match"):
            bad_inchi = "InChI=1/C14H13ClF5N5/c15-12-10(6-25(24-12)8(3-16)4-17)22-13-21-5-9(14(18,19)20)11(23-13)7-1-2-7/h5-8H,1-4H2,(H,21,22,23)C12"
            self.assertEqual("", check_dtx_for_inchi(bad_inchi))

        with self.subTest("multiple matches, single known & single unknown"):
            # DN0021363 is single unknown stereoisomer, DN0021769 is single known stereoisomer
            inchi = "InChI=1/C15H14F2N2O2/c1-8-3-4-10(19-13(20)12-6-15(12,16)17)5-11(8)14-18-9(2)7-21-14/h3-5,7,12H,6H2,1-2H3,(H,19,20)/t12-/m0/s1"
            # No mixtures, only one single unknown -> keep the single unknown
            self.assertEqual("DN0021363", check_dtx_for_inchi(inchi))

        with self.subTest("multiple matches, single unknown & mixture"):
            # DN0018288 is single unknown stereoisomer,  DN0017477 is mixture of enantiomers
            inchi = "InChI=1/C17H18N4O/c1-13(14-7-9-19-10-8-14)21(2)17(22)20-12-16-6-4-3-5-15(16)11-18/h3-10,13H,12H2,1-2H3,(H,20,22)/t13-/m0/s1"
            # Keep the mixture
            self.assertEqual("DN0017477", check_dtx_for_inchi(inchi))

        with self.subTest(
            "multiple matches, mixture of enantiomers & mixture of diastereomers"
        ):
            # DN0023392 is mixture of enantiomers, DN0022385 is mixture of diastereomers
            inchi = "InChI=1/C13H16N2O3S/c1-10-9-18-11(2)8-15(10)19(16,17)13-5-3-4-12(6-13)7-14/h3-6,10-11H,8-9H2,1-2H3/t10?,11?"
            # Both are mixtures -> keep both
            self.assertEqual("DN0022385,DN0023392", check_dtx_for_inchi(inchi))
