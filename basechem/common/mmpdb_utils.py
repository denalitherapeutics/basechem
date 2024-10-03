import datetime
import os
import shutil
import subprocess

from django.conf import settings
from django.core.mail import mail_admins

from basechem.common.constants import ADMIN_FAILURE


def get_mmpdb_postgres_db():
    """
    Return mmpdb postgres host url
    """
    pw = settings.DATABASES["default"]["PASSWORD"]
    host = settings.DATABASES["default"]["HOST"]
    port = settings.DATABASES["default"]["PORT"]
    user = settings.DATABASES["default"]["USER"]
    env = "prod" if settings.ENVIRONMENT == "prod" else "test"

    return f"postgresql://{user}:{pw}@{host}:{port}/mmpdb_bc_{env}"


def initialize_mmpdb(db_env):
    """
    Set up a new MMPDB database with the provided test data if db_env is "test"
    """
    if db_env == "test":
        test_smi = "basechem/common/mmpdb_data/test_data.smi"
        frag_db = _fragment_mmpdb(test_smi)
        if frag_db:
            indexed = _index_mmpdb(frag_db)
            return indexed


def _fragment_mmpdb(smi_path):
    """
    Helper to fragment the given smiles file in MMPDB, requires cache.fragdb already be set up
    :param smi_path: path to smiles file to fragment
    :return: path to file with newly fragmented data if it exists
    """
    date_str = datetime.datetime.today().strftime("%Y%m%d")
    # We use the matcher_alpha cut-smarts which is custom from matcher but duplicated in our MMPDB install
    # it doesn't need to be added here since the cache.fragdb file already has information about the cut-smarts
    cache_frag = f"basechem/common/mmpdb_data/cache.fragdb"
    new_frag = f"/tmp/{date_str}_new.fragdb"
    if os.path.exists(cache_frag):
        frag_cmd = [
            "mmpdb",
            "fragment",
            smi_path,
            "--cache",
            cache_frag,
            "-o",
            new_frag,
        ]
    else:
        frag_cmd = [
            "mmpdb",
            "fragment",
            smi_path,
            "-o",
            cache_frag,
            "--cut-smarts",
            "matcher_alpha",
        ]

    process = subprocess.run(
        args=frag_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if process.returncode == 0 and os.path.exists(new_frag):
        shutil.move(new_frag, cache_frag)
        return cache_frag

    mail_admins(
        ADMIN_FAILURE,
        f"MMPDB Fragmenting failed for {new_frag}.\nSTDOUT:\n{process.stdout.decode()}\nSTDERR:\n{process.stderr.decode()}",
    )


def _index_mmpdb(frag_db):
    """
    Helper to index the given fragment db in MMPDB, mmpdb will drop all tables if they exist already
    and recreate the databse from scratch
    :param frag_db: the fragments to index
    :return: True if the process completes successfully
    """
    if _db_exists():
        # This is a weird bug in MMPDB where it deletes the db but doesn't remake if it exists already
        # so if the DB exists already, run twice to repopulate
        index_cmd = ["mmpdb", "index", frag_db, "-o", get_mmpdb_postgres_db()]
        process = subprocess.run(
            args=index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

    index_cmd = ["mmpdb", "index", frag_db, "-o", get_mmpdb_postgres_db()]
    process = subprocess.run(
        args=index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if process.returncode == 0:
        return True

    mail_admins(
        ADMIN_FAILURE,
        f"MMPDB Indexing failed.\nSTDOUT:\n{process.stdout.decode()}\nSTDERR:\n{process.stderr.decode()}",
    )


def _db_exists():
    """
    Return true if a current MMPDB exists
    """
    db = get_mmpdb_postgres_db()
    list_cmd = ["mmpdb", "list", db]
    process = subprocess.run(
        args=list_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if db in process.stdout.decode():
        return True
    return False


def generate_mmpdb(smiles, output_filepath, radius=0):
    """
    Generate MMPs for the given smiles string using mmpdb
    :param smiles: the smiles of the compound to generate MMPs for
    :param radius: the radius argument for generate. Each available radius has a corresponding rule environment. Larger radii have more specific environments. Default 0.
    :output_filepath: the filepath where the generate file will be stored
    :return: True if generate command was successful
    """
    generate_cmd = [
        "mmpdb",
        "generate",
        get_mmpdb_postgres_db(),
        "--smiles",
        smiles,
        "-o",
        output_filepath,
        "--radius",
        str(radius),
    ]
    process = subprocess.run(
        args=generate_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    if process.returncode == 0:
        return True

    mail_admins(
        ADMIN_FAILURE,
        f"MMPDB Generate failed for smiles '{smiles}'.\nSTDOUT:\n{process.stdout.decode()}\nSTDERR:\n{process.stderr.decode()}",
    )


def loadprops_mmpdb(propfile_path):
    """
    Load properties into mmpdb for the IDs/props in the given file
    If the given property already exists, the value will be updated
    :param propfile_path: path to tsv file containing properties to load
    :return: True if generate command was successful
    """
    loadprops_cmd = [
        "mmpdb",
        "loadprops",
        "--properties",
        propfile_path,
        get_mmpdb_postgres_db(),
    ]
    process = subprocess.run(
        args=loadprops_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if process.returncode == 0:
        return True

    mail_admins(
        ADMIN_FAILURE,
        f"MMPDB loadprops failed for file {propfile_path}.\nSTDOUT:\n{process.stdout.decode()}\nSTDERR:\n{process.stderr.decode()}",
    )
