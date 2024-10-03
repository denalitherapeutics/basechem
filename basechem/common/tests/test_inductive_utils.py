from django.core import mail
from django.test import tag

from basechem.common.inductive_utils import (
    run_inductive_alogd_predict,
    run_inductive_lm_predict,
    update_inductive_logd_data,
)
from basechem.common.mocks.mock_inductive_utils import (
    mock_run_inductive_alogd_predict,
    mock_run_inductive_lm_predict,
)
from basechem.common.tests.base import BasechemNoMockTestCase


@tag("local", "inductive", "external")
class InductiveBioUtilsTestCase(BasechemNoMockTestCase):
    def test_run_inductive_lm_predict(self):
        """
        Tests `run_inductive_lm_predict` and tests that output of `mock_run_inductive_lm_predict`
        tracks with `run_inductive_lm_predict`
        """
        input_file = "basechem/main/tests/testdata/test_onecomp.sdf"
        with self.subTest("Rat, no image"):
            result = run_inductive_lm_predict(input_file, "R", False)
            mock_result = mock_run_inductive_lm_predict(input_file, "R", False)

            # Check that lists are the same length, then compare the first element
            self.assertEqual(len(result), len(mock_result))
            result = result[0]
            mock_result = mock_result[0]

            # "prediction" can be cast integer
            self.assertEqual(int, type(int(result["prediction"])))
            # "out_of_domain" can be cast as a bool
            self.assertEqual(bool, type(bool(result["out_of_domain"])))
            # "interp_image" should be the empty string
            self.assertEqual(result["interp_image"], "")
            # "probs_image" should be the empty string
            self.assertEqual(result["probs_image"], "")

            # Check that the mock result has the correct format
            self.assertEqual(mock_result.keys(), result.keys())
            for key in mock_result.keys():
                self.assertEqual(type(mock_result[key]), type(result[key]))
            self.assertEqual(mock_result["interp_image"], "")
            self.assertEqual(mock_result["probs_image"], "")

        with self.subTest("Rat, image"):
            result = run_inductive_lm_predict(input_file, "R", True)
            mock_result = mock_run_inductive_lm_predict(input_file, "R", True)

            # Check that lists are the same length, then compare the first element
            self.assertEqual(len(result), len(mock_result))
            result = result[0]
            mock_result = mock_result[0]

            # "prediction" can be cast integer
            self.assertEqual(int, type(int(result["prediction"])))
            # "out_of_domain" can be cast as a bool
            self.assertEqual(bool, type(bool(result["out_of_domain"])))
            # images should be base64 images of at least 50 chars
            self.assertEqual(str, type(result["interp_image"]))
            self.assertTrue(len(result["interp_image"]) > 50)
            self.assertEqual(str, type(result["probs_image"]))
            self.assertTrue(len(result["probs_image"]) > 50)

            # Check that the mock result has the correct format
            self.assertEqual(mock_result.keys(), result.keys())
            for key in mock_result.keys():
                self.assertEqual(type(mock_result[key]), type(result[key]))
            self.assertTrue(len(mock_result["interp_image"]) > 50)
            self.assertTrue(len(mock_result["probs_image"]) > 50)

        with self.subTest("Human, no image"):
            result = run_inductive_lm_predict(input_file, "H", False)
            mock_result = mock_run_inductive_lm_predict(input_file, "H", False)

            # Check that lists are the same length, then compare the first element
            self.assertEqual(len(result), len(mock_result))
            result = result[0]
            mock_result = mock_result[0]

            # "prediction" can be cast integer
            self.assertEqual(int, type(int(result["prediction"])))
            # "out_of_domain" can be cast as a bool
            self.assertEqual(bool, type(bool(result["out_of_domain"])))
            # "interp_image" should be the empty string
            self.assertEqual(result["interp_image"], "")
            # "probs_image" should be the empty string
            self.assertEqual(result["probs_image"], "")

            # Check that the mock result has the correct format
            self.assertEqual(mock_result.keys(), result.keys())
            for key in mock_result.keys():
                self.assertEqual(type(mock_result[key]), type(result[key]))
            self.assertEqual(mock_result["interp_image"], "")
            self.assertEqual(mock_result["probs_image"], "")

        with self.subTest("Human, image"):
            result = run_inductive_lm_predict(input_file, "H", True)
            mock_result = mock_run_inductive_lm_predict(input_file, "H", True)

            # Check that lists are the same length, then compare the first element
            self.assertEqual(len(result), len(mock_result))
            result = result[0]
            mock_result = mock_result[0]

            # "prediction" can be cast integer
            self.assertEqual(int, type(int(result["prediction"])))
            # "out_of_domain" can be cast as a bool
            self.assertEqual(bool, type(bool(result["out_of_domain"])))
            # images should be base64 images of at least 50 chars
            self.assertEqual(str, type(result["interp_image"]))
            self.assertTrue(len(result["interp_image"]) > 50)
            self.assertEqual(str, type(result["probs_image"]))
            self.assertTrue(len(result["probs_image"]) > 50)

            # Check that the mock result has the correct format
            self.assertEqual(mock_result.keys(), result.keys())
            for key in mock_result.keys():
                self.assertEqual(type(mock_result[key]), type(result[key]))
            self.assertTrue(len(mock_result["interp_image"]) > 50)
            self.assertTrue(len(mock_result["probs_image"]) > 50)

    def test_run_inductive_alogd_predict(self):
        """
        Tests `run_inductive_alogd_predict` and tests that output of `mock_run_inductive_alogd_predict`
        tracks with `run_inductive_alogd_predict`
        """
        input_file = "basechem/main/tests/testdata/test_onecomp.sdf"

        result = run_inductive_alogd_predict(input_file)
        mock_result = mock_run_inductive_alogd_predict(input_file)

        self.assertEqual(len(result), len(mock_result))
        result = result[0]
        mock_result = mock_result[0]

        self.assertEqual(float, type(float(result["prediction"])))
        self.assertEqual(mock_result.keys(), result.keys())
        for key in mock_result.keys():
            self.assertEqual(type(mock_result[key]), type(result[key]))

    def test_update_inductive_logd_data(self):
        """
        Tests `update_inductive_logd_data` to make sure no email is sent with
        empty data
        """
        self.assertEqual(len(mail.outbox), 0)
        input_file = "/tmp/empty.sdf"
        with open(input_file, "w") as f:
            pass

        update_inductive_logd_data(input_file)

        self.assertEqual(len(mail.outbox), 0)
