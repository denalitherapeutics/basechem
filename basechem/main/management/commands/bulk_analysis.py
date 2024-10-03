import datetime
import os

from django.core.files.base import File
from django_q.tasks import async_task
from rdkit import Chem

from basechem.common.rdkit_utils import RDKitWrappers
from basechem.main.constants import ALIGN, DOCK
from basechem.main.models.collection_models import Collection, collection_files_path
from basechem.main.models.project_models import Project
from basechem.mni_common.storage import save_media_file
from basechem.users.models import BasechemUser


def run_bulk_analysis(infile, project_code, analysis):
    """
    Run provided analysis for all the compounds in the given infile
    :param infile: sdf file containing compounds to analyze
    :param project_code: project code for a project already created in basechem
    :param analysis: the type of analysis to run; dock or align
    """
    name = (
        f"bulk_{analysis}_{project_code}_{datetime.datetime.now().strftime('%m%d%Y')}"
    )

    if analysis == DOCK:
        collection = _bulk_collection_setup(infile, project_code, name)
        async_task(
            _bulk_dock,
            collection=collection,
            name=name,
            task_name=name,
            hook="basechem.main.tasks.simple_email_task_hook",
        )
    elif analysis == ALIGN:
        collection = _bulk_collection_setup(infile, project_code, name)
        async_task(
            _bulk_align,
            collection=collection,
            name=name,
            task_name=name,
            hook="basechem.main.tasks.simple_email_task_hook",
        )
    else:
        print(f"{analysis} is not an option, please use '{DOCK}' or '{ALIGN}'")


def _bulk_dock(collection, name):
    """
    Dock the compounds in the given collection
    :param collection: the collection with compounds to dock
    :param name: the name of the task
    :return: path to sdf file with docked poses, dict of {DN_ID: error}
    """
    out_filepath = collection_files_path(collection, f"{name}.sdf", local=True)
    error_dict = {}

    w = Chem.SDWriter(open(out_filepath, "a+"))

    for co in collection.get_cos_for_analysis(DOCK):
        try:
            poses = co.dock_to_receptor(reference=None)
            for conf_id, pose in poses.items():
                mol = Chem.MolFromMolBlock(pose["moltext"])
                mol.SetProp("conf_id", conf_id)
                del pose["moltext"]
                for k, v in pose.items():
                    mol.SetProp(k, str(v))
                w.write(mol)

        except Exception as e:
            error_dict[co.compound.name] = e
            continue
    w.close()

    if os.path.exists(out_filepath):
        media_path = collection_files_path(collection, f"{name}.sdf", local=False)
        # Save file to S3
        with open(out_filepath, "rb") as fp:
            save_media_file(media_path, fp)

    return out_filepath, error_dict


def _bulk_align(collection, name):
    """
    Align the compounds in the given collection
    :param collection: the collection with compounds to align
    :param name: the name of the task
    :return: path to sdf file with aligned poses, dict of {DN_ID: error}
    """
    out_filepath = collection_files_path(collection, f"{name}.sdf", local=True)
    error_dict = {}

    w = Chem.SDWriter(open(out_filepath, "a+"))

    for co in collection.get_cos_for_analysis(ALIGN):
        try:
            poses = co.superimpose_to_ref(reference=None)
            for conf_id, pose in poses.items():
                mol = Chem.MolFromMolBlock(pose["moltext"])
                mol.SetProp("conf_id", conf_id)
                del pose["moltext"]
                for k, v in pose.items():
                    mol.SetProp(k, str(v))
                w.write(mol)

        except Exception as e:
            error_dict[co.compound.name] = e
            continue
    w.close()

    if os.path.exists(out_filepath):
        media_path = collection_files_path(collection, f"{name}.sdf", local=False)
        # Save file to S3
        with open(out_filepath, "rb") as fp:
            save_media_file(media_path, fp)

    return out_filepath, error_dict


def _bulk_collection_setup(infile, project_code, name):
    """
    Process compounds from the given sdf file into Compounds and make a collection
    :param infile: sdf file with compounds to run analysis on
    :param project_code: project code for the compounds in the file
    :param name: name of the collection
    :return: new collection
    """
    collection = Collection(
        name=name,
        project=Project.objects.get(code=project_code),
        owner=BasechemUser.objects.get(first_name="ADMIN"),
    )
    collection.save()
    with open(infile, "rb") as f:
        collection.sdf_file.save(os.path.basename(infile), File(f))
        collection.save()

    suppl = Chem.SDMolSupplier(infile, removeHs=False)
    all_mols = [m for m in suppl]
    romols = [(m, RDKitWrappers.is_3d(m)) for m in all_mols if m is not None]

    collection.handle_romols(romols)
    return collection
