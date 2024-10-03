import os
import subprocess

from django.conf import settings
from rdkit import Chem

from basechem.common.analysis_utils import (
    get_torsion_commands,
    process_torsion_job_result,
)
from basechem.common.slurm_utils import get_bash_script_from_commands


def mock_submit_torsion_job_to_slurm(
    job_name,
    dihedral,
    dihedral_smarts,
    dihedral_atoms,
    input_filename,
    output_filename,
    slurm_working_dir,
):
    """
    Intended to patch basechem.common.analysis_utils.submit_torsion_job_to_slurm and so it
    should match the signature of that function. Runs the Mayachemtools torsion script
    locally instead of submitting a job to the Slurm cluster.
    """
    container_working_dir = f"/opt/shared_files/local_testsuite/{'/'.join(slurm_working_dir.split('/')[-3:])}/"
    commands = get_torsion_commands(
        dihedral_smarts, dihedral_atoms, input_filename, output_filename, dihedral
    )
    script = get_bash_script_from_commands(
        [f"export PATH=$PATH:{settings.MAYACHEMTOOLS_DIR}/bin"] + commands
    )
    tmp_script_file = os.path.join(container_working_dir, f"run_torsion_{dihedral}.sh")
    with open(tmp_script_file, "w+") as script_file:
        script_file.write(script)
    subprocess.call(
        ["bash", tmp_script_file],
        cwd=container_working_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
    # Return dihedral or 1 so that the dihedral 0 doesn't return a falsey job id
    return dihedral or 1


def mock_submit_torsion_job_to_slurm_USE_CACHE(
    job_name,
    dihedral,
    dihedral_smarts,
    dihedral_atoms,
    input_filename,
    output_filename,
    slurm_working_dir,
):
    """
    Intended to patch basechem.common.analysis_utils.submit_torsion_job_to_slurm and so it
    should match the signature of that function. Uses cached results instead of running
    the Mayachemtools torsion script so that some tests can run significantly faster.
    """
    cached_output_path = "basechem/main/tests/testdata/test_torsion_results.sdf"
    for mol in Chem.SDMolSupplier(cached_output_path):
        if int(mol.GetProp("Torsion_Angle")) == dihedral:
            with Chem.SDWriter(output_filename) as w:
                w.write(mol)
    # Return dihedral or 1 so that the dihedral 0 doesn't return a falsey job id
    return dihedral or 1


def mock_submit_torsion_job_to_slurm_fail_slurm_api(
    job_name,
    dihedral,
    dihedral_smarts,
    dihedral_atoms,
    input_filename,
    output_filename,
    slurm_working_dir,
):
    """
    Intended to patch basechem.common.analysis_utils.submit_torsion_job_to_slurm and so it
    should match the signature of that function. Runs the Mayachemtools torsion script
    locally instead of submitting a job to the Slurm cluster, failing to start jobs for dihedrals that are multiples of 20.
    This is used to test torsion jobs that only return a fraction of the results.
    """
    if not dihedral % 20:
        # Simulate job failing to start by returning a falsey value for job ID
        return
    else:
        return mock_submit_torsion_job_to_slurm_USE_CACHE(
            job_name,
            dihedral,
            dihedral_smarts,
            dihedral_atoms,
            input_filename,
            output_filename,
            slurm_working_dir,
        )


def mock_process_torsion_job_result_fail_psi4(
    input_filename, output_filename, job_name, container_working_dir
):
    """
    Intended to patch basechem.common.analysis_utils.process_torsion_job_result and so it
    should match the signature of that function. Simulates Psi4 failures for all dihedrals that
    are not multiples of 90.
    """
    dihedral = int(os.path.splitext(input_filename)[0].split("_")[-1])
    if dihedral % 90:
        raise Exception("Job failed")
    else:
        return process_torsion_job_result(
            input_filename, output_filename, job_name, container_working_dir
        )


def mock_has_job_completed(job_id):
    """
    This is intended to be used alongside mock_submit_torsion_job_to_slurm to patch
    basechem.common.analysis_utils.has_job_completed.
    """
    # Always return True because the jobs were run synchronously locally instead of on slurm
    return True


def mock_generate_torsion_alerts_NOOP(input_filepath):
    """
    No-op mock for `generate_torsion_alerts` to use in all test cases except those that directly test the mayachem script.
    This mock simply returns the input file path (without adding torsion alerts) because the mayachem code cannot be run in GitHub Actions.
    :param input_path: a string, sdf file that contains conformers to generate torsion alerts for
    :return: a string, filepath to the output file with torsion alert properties added (when applicable)
    """
    return input_filepath
