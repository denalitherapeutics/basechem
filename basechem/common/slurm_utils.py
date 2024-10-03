import logging
import os
import subprocess

import requests
from django.conf import settings
from django.core.mail import mail_admins

from basechem.common.constants import (
    ADMIN_FAILURE,
    COMPLETING,
    CONFIGURING,
    LAUNCH_FAILED_DESC,
    LAUNCH_FAILED_REASON,
    PENDING,
    RUNNING,
    SLURM_REST_API_USER,
    UNKNOWN,
)

logger = logging.getLogger("django")


def generate_slurm_api_token():
    """
    Generates a JWT token for use with the Slurm Rest API and adds it to os.environ.
    :return: a JWT token (str)
    """
    output = subprocess.check_output(
        [
            "ssh",
            "-i",
            settings.SLURM_SSH_KEY_PATH,
            "-o",
            "StrictHostKeyChecking=no",
            f"{SLURM_REST_API_USER}@{settings.SLURM_REST_API_HOST}",
            "/opt/slurm/bin/scontrol token lifespan=999999",
        ]
    )
    # output is in the format SLURM_JWT=TOKEN
    token = output.decode().split("=")[1].strip("\n")
    os.environ["SLURM_API_TOKEN"] = token
    return token


def get_or_create_slurm_api_token():
    """
    Retrieves the JWT token to be used with the Slurm Rest API or generates a new
    one if it is not found in os.environ.
    :return: a JWT token (str)
    """
    token = os.environ.get("SLURM_API_TOKEN")
    if not token:
        token = generate_slurm_api_token()

    return token


def get_bash_script_from_commands(commands):
    """
    Turns the specified list of bash commands into a bash script.
    :param commands: a list of bash commands to run
    :return: a bash script (str) that encodes the list of bash commands to run
    """
    return f"#!/bin/bash\n" + "\n".join(commands)


def run_on_slurm_node(job_name, commands, working_dir):
    """
    :param job_name: name (str) assigned to the Slurm job
    :param commands: a list of bash commands to run
    :param working_dir: path to the directory (str) on the Slurm node on which the job will be run
    :return job_id: id (str) of the job submitted to the Slurm cluster
    """
    # The path_to_create tmp directory is owned by the root user since it is created from
    # the container so we need to change permissions such that it can be written to from within
    # a Slurm cluster node. SLURM_SHARED_FILES_TMP_DIR should have write
    # permissions to the Slurm user such that we can match the permissions on working_dir_relative_to_slurm
    # to those on the SLURM_SHARED_FILES_TMP_DIR.
    setup_commands = [
        f"sudo chown * -Rv --reference={settings.SLURM_SHARED_FILES_TMP_DIR}",
        f"cd {working_dir}",
    ]
    script = get_bash_script_from_commands(setup_commands + commands)
    response = requests.post(
        f"{settings.SLURM_REST_API_BASE_URL}/job/submit",
        json={
            "job": {
                "name": job_name,
                "ntasks": 1,
                "nodes": 1,
                "partition": "cpu-q",
                "current_working_directory": settings.SLURM_SHARED_FILES_TMP_DIR,
                "standard_input": "/dev/null",
                "standard_output": os.path.join(
                    settings.SLURM_SHARED_FILES_TMP_DIR, f"{job_name}.log"
                ),
                "standard_error": os.path.join(
                    settings.SLURM_SHARED_FILES_TMP_DIR, f"{job_name}_error.log"
                ),
                "environment": {
                    "PATH": f"/bin:/usr/bin/:/usr/local/bin/:{settings.MAYACHEMTOOLS_DIR}/bin",
                    "LD_LIBRARY_PATH": "/lib/:/lib64/:/usr/local/lib",
                },
            },
            "script": script,
        },
        headers={
            "X-SLURM-USER-NAME": SLURM_REST_API_USER,
            "X-SLURM-USER-TOKEN": get_or_create_slurm_api_token(),
        },
    )

    if response.status_code == 200:
        response_json = response.json()
        job_id = response_json["job_id"]
        return job_id
    else:
        message = f"Failed to submit the job {job_name} to the Slurm cluster.\n{script}"
        logger.error(message)
        mail_admins(ADMIN_FAILURE, message)


def get_job_data(job_id):
    """
    Retrieves the job data from the Slurm Rest API.
    :param job_id: id of the job submitted to the Slurm cluster
    :return: a dictionary of job information (includes keys 'start_time', 'job_state', 'submit_time', etc.)
    """
    response = requests.get(
        f"{settings.SLURM_REST_API_BASE_URL}/job/{job_id}",
        headers={
            "X-SLURM-USER-NAME": "basechem",
            "X-SLURM-USER-TOKEN": get_or_create_slurm_api_token(),
        },
    )
    if response.status_code == 200:
        jobs = response.json()["jobs"]
        for job in jobs:
            if job["job_id"] == job_id:
                return job
    return {}


def get_job_state(job_id):
    """
    Retrieves the job state from the Slurm Rest API.
    :param job_id: id of the job submitted to the Slurm cluster
    :return: one of "RUNNING", "FAILED", "COMPLETED", or "UNKNOWN" based on the Slurm
        API response
    """
    return get_job_data(job_id).get("job_state", UNKNOWN)


def has_job_completed(job_id):
    """
    Determines if the job is still running or scheduled to run at a later time.
    :param job_id: id of the job submitted to the Slurm cluster
    :return: True if a job is no longer in a "PENDING", "RUNNING", "CONFIGURING", or "COMPLETING", state
    """
    return get_job_state(job_id) not in [PENDING, RUNNING, CONFIGURING, COMPLETING]


def has_job_stalled(job_id):
    """
    Determines if the job has stalled and emails the django admins if so. A stalled job is one
    that is pending due to launch failure.
    :param job_id: id of the job submitted to the Slurm cluster
    :return: True if a job is stalled
    """
    job_data = get_job_data(job_id)
    stalled = (
        job_data["job_state"] == PENDING
        and job_data.get("state_reason") == LAUNCH_FAILED_REASON
        and job_data.get("state_description") == LAUNCH_FAILED_DESC
    )
    if stalled:
        message = f"The slurm job {job_id} is pending due to launch failure. This job will remain in the slurm queue until it is released by an admin. The Basechem task that started this job will not wait for its result."
        mail_admins(ADMIN_FAILURE, message)
    return stalled
