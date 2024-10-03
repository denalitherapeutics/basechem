import os
from unittest.mock import Mock, patch

from basechem.common.constants import FAILED, RUNNING, UNKNOWN
from basechem.common.slurm_utils import get_job_state, run_on_slurm_node
from basechem.common.tests.base import BasechemTestCase


@patch("basechem.common.slurm_utils.generate_slurm_api_token")
class TestSlurmUtils(BasechemTestCase):
    @patch("requests.post")
    def test_job_submission_returns_job_id(
        self, mock_requests_post, mock_generate_slurm_api_token
    ):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_requests_post.return_value = mock_response
        mock_response.json.return_value = {"job_id": 1}
        job_id = run_on_slurm_node("test", ["echo 'test sleep job'"], os.getcwd())
        self.assertEqual(mock_requests_post.call_count, 1)
        self.assertEqual(job_id, 1)

    @patch("requests.get")
    def test_get_job_state(self, mock_requests_get, mock_generate_slurm_api_token):
        job_id = 1
        mock_response = Mock()
        mock_response.status_code = 200
        mock_requests_get.return_value = mock_response

        with self.subTest("running"):
            mock_response.json.return_value = {
                "jobs": [{"job_id": job_id, "job_state": RUNNING}]
            }
            mock_requests_get.status_code = 200
            job_state = get_job_state(job_id)
            self.assertEqual(mock_requests_get.call_count, 1)
            self.assertEqual(job_state, RUNNING)

        with self.subTest("failed"):
            mock_response.json.return_value = {
                "jobs": [{"job_id": job_id, "job_state": FAILED}]
            }
            mock_requests_get.status_code = 200
            job_state = get_job_state(job_id)
            self.assertEqual(job_state, FAILED)

        with self.subTest("unknown"):
            mock_response.status_code = 500
            mock_requests_get.status_code = 200
            job_state = get_job_state(job_id)
            self.assertEqual(job_state, UNKNOWN)
