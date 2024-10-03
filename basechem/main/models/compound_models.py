import base64
import logging
import os
import re
from collections import Counter, OrderedDict
from pathlib import Path

from django.conf import settings
from django.contrib.postgres.indexes import GistIndex
from django.core.exceptions import ObjectDoesNotExist, ValidationError
from django.core.files.base import File
from django.core.files.storage import FileSystemStorage
from django.core.mail import mail_admins
from django.core.validators import FileExtensionValidator
from django.db import connection, models
from django.utils.safestring import mark_safe
from django.utils.timezone import now
from django_q.models import OrmQ, Task
from django_q.tasks import async_task
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdMolTransforms

from basechem.common.analysis_utils import (
    copy_rdock_ligand_prm,
    generate_rdock_grid,
    generate_rdock_system_prm,
    generate_toklat_scores,
    generate_torsion_alerts,
    mol_to_toklat_annotations,
    mol_to_torsion_alerts,
    run_apbs,
    run_esp_predict,
    run_lsalign,
    run_mc_torsion_scan,
    run_rdock,
)
from basechem.common.constants import ADMIN_NOTIFICATION
from basechem.common.dtx_utils import get_agg_ic50_data, get_registered_structures
from basechem.common.file_utils import get_tmp_file
from basechem.common.inductive_utils import (
    run_inductive_alogd_predict,
    run_inductive_lm_predict,
)
from basechem.common.rdkit_utils import RDKitWrappers
from basechem.main.constants import *
from basechem.main.db_utils import MolField
from basechem.main.mmp_utils import generate_mmps
from basechem.main.utils import strip_series_from_conf_id
from basechem.mni_common.rdkit_utils import moltext_to_svg
from basechem.mni_common.storage import save_media_file, select_media_storage
from basechem.users.models import BasechemUser

logger = logging.getLogger("django")


def series_files_path(instance, filename, local=False):
    """
    Directory for placing project level files such as series files, from MediaRoot
    """
    # Check if the filename contains extra directory information to split off
    if len(filename.split("/")) > 1:
        diradd = os.path.dirname(filename)
        filename = filename.split("/")[-1]
        dir_path = f"project_{instance.project.code}/{instance.dn_id}/{diradd}"
    else:
        dir_path = f"project_{instance.project.code}/{instance.dn_id}"
    if local:
        dir_path = f"{settings.MEDIA_ROOT}/{dir_path}"
        os.makedirs(dir_path, exist_ok=True)

    return f"{dir_path}/{filename}"


def compound_files_path(instance, filename, co=None, local=False):
    """
    Directory for placing generated files for compounds, from MediaRoot
    :param instance: the Compound instance
    :param filename: the name of the file
    :param co: True if a CompoundOccurrence subdirectory is needed
    :param local: True if the local filepath should be returned/created;
        False if the relative path should be returned
    :return: relative or complete filepath with filename
    """
    if co:
        dir_path = f"compounds/{instance.virtual_id}/co_{co.pk}"
    else:
        dir_path = f"compounds/{instance.virtual_id}"
    if local:
        dir_path = f"{settings.MEDIA_ROOT}/{dir_path}"
        os.makedirs(dir_path, exist_ok=True)

    return f"{dir_path}/{filename}"


class Series(models.Model):
    """
    Class to create series and save key information for use in analysis
    Instances of this model will primarily be set and updated in
    the admin panel
    """

    dn_id = models.CharField(max_length=11)
    name = models.CharField(max_length=50)

    project = models.ForeignKey("main.Project", on_delete=models.PROTECT, null=True)
    # Assay data of the form {ASSAY_NAME: IC50_VALUE}
    assay_data = models.JSONField(default=dict, blank=True, null=True)
    ground_state_file = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["sdf"])],
    )
    bound_state_file = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["sdf"])],
    )
    receptor_file = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["pdb"])],
    )

    receptor_residues = models.CharField(max_length=500, blank=True)
    active = models.BooleanField(default=True)
    smiles = models.CharField(max_length=500)

    # These fields are required for docking jobs using rdock - the rdock_ fields are auto set but can be updated manually
    receptor_file_mol2 = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["mol2"])],
    )
    rdock_prm = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["prm"])],
        help_text="Optional: This will be auto-generated on the first rDock run for this series",
    )
    rdock_grid = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["grd"])],
        help_text="Optional: This will be auto-generated on the first rDock run for this series",
    )
    rdock_as = models.FileField(
        upload_to=series_files_path,
        storage=select_media_storage,
        blank=True,
        validators=[FileExtensionValidator(allowed_extensions=["as"])],
        help_text="Optional: This will be auto-generated on the first rDock run for this series",
    )

    class Meta:
        verbose_name_plural = "Series"

    def save(self, *args, **kwargs):
        """
        If the Series is not active, re-assign Series to all the Compound
        objects currently assigned to this Series
        """
        super().save(*args, **kwargs)
        if not self.active:
            task_name = f"reassign_series_{self.id}_{self.name}"
            queued = [
                oq for oq in OrmQ.objects.all() if oq.task.get("name") == task_name
            ]
            if not queued:
                async_task(self.reassign_compounds, task_name=task_name)

    def reassign_compounds(self):
        """
        For each Compound assigned to this Series, pick a new Series. This function is called
        asynchronously from the Series save method when a Series is marked as inactive. Additionally,
        sends an admin email with stats about the reassigned Compounds
        """
        series_counter = Counter()
        unassigned = []
        for c in self.compound_set.all():
            try:
                new_series = c.pick_series()
                if new_series:
                    series_counter[new_series] += 1
                    c.series = new_series
                    c.save()
                else:
                    unassigned.append(c)
            except:
                unassigned.append(c)
        # Send an email to Django admins with information about reassigned Compounds
        num_reassigned = sum(series_counter.values())
        message = f"Reassigned Compounds assigned to Series {self.name} (PK: {self.pk}). {num_reassigned}/{num_reassigned + len(unassigned)} Compounds were assigned to new Series:\n"
        for series, count in series_counter.most_common():
            message += f"- {count} assigned to {series.name} (PK: {series.pk})\n"
        if unassigned:
            message += "\n\nThe following Compounds (PKs) could not have a new Series assigned:\n"
            message += ", ".join([str(c.pk) for c in unassigned])

        mail_admins(ADMIN_NOTIFICATION, message)

    def mol(self):
        """
        Creates a mol object using the conformation in the `default_file` (either the bound state
        or ground state file).
        :returns: an RDKit mol object
        """
        localpath = self.default_file(path=True)
        mol = Chem.MolFromMolFile(localpath)
        mol.SetProp("_Name", self.dn_id)
        return mol

    def pdb_block(self):
        """
        Return PDB block of text of receptor file
        """
        if self.receptor_file:
            localpath = self.path_to_local_file(self.receptor_file)
            mol = Chem.MolFromPDBFile(localpath)
            return Chem.MolToPDBBlock(mol)
        return ""

    def default_file(self, path=False):
        """
        Use bound_state_file as a default if it exists
        :param path: if True, return path to localfile
        :return: file field or path to localfile
        """
        default = None
        if self.bound_state_file:
            default = self.bound_state_file
        elif self.ground_state_file:
            default = self.ground_state_file

        if default and path:
            fname = os.path.split(default.name)[-1]
            localpath = series_files_path(self, fname, local=True)
            if not os.path.exists(localpath):
                content = default.read()
                with open(localpath, "wb") as fw:
                    fw.write(content)

            return localpath
        return default

    def path_to_local_file(self, field):
        """
        Returns path to local file for the given field, will
        save files from model fields if they don't exist locally
        :param field: one of this series' FileFields (ex: series.receptor_file)
        :return: localpath of the given fields file
        """
        fname = os.path.split(field.name)[-1]
        localpath = series_files_path(self, fname, local=True)

        if not os.path.exists(localpath):
            content = field.read()
            with open(localpath, "wb") as fw:
                fw.write(content)

        return localpath

    ######################
    #### DOCK METHODS ####
    ######################

    def get_rdock_prm(self):
        """
        Generates a parameter file and saves it to the series if it doesn't exist
        :return: path to local rDock prm grid parameter file
        """
        expected_prm_path = series_files_path(
            self, f"{self.dn_id}_rdock.prm", local=True
        )

        # Save needed receptor files locally
        local_mol2 = self.path_to_local_file(self.receptor_file_mol2)
        local_bound = self.path_to_local_file(self.bound_state_file)

        # If prm file doesn't exist anywhere, generate a local copy and save to model
        if not self.rdock_prm and not os.path.exists(expected_prm_path):
            _ = generate_rdock_system_prm(local_mol2, local_bound, self.dn_id)
            path = Path(expected_prm_path)
            with path.open(mode="rb") as f:
                self.rdock_prm = File(f, name=path.name)
                self.save()

        # If prm file doesn't exist locally, save local copy
        if not os.path.exists(expected_prm_path):
            content = self.rdock_prm.read()
            with open(expected_prm_path, "wb") as fw:
                fw.write(content)

        self._get_rdock_grid(expected_prm_path)

        return expected_prm_path

    def _get_rdock_grid(self, rdock_prm_file):
        """
        Generates .grid and .as files for rdock jobs and saves them to the series if they
        doesn't exist
        :rdock_prm_file: localpath to rdock parameter file
        """
        rdock_filename = os.path.splitext(os.path.basename(rdock_prm_file))[0]
        grid_filepath = series_files_path(
            self, f"{rdock_filename}_cav1.grd", local=True
        )
        as_filepath = series_files_path(self, f"{rdock_filename}.as", local=True)

        if not self.rdock_grid:
            generate_rdock_grid(rdock_prm_file)
            if os.path.exists(grid_filepath) and os.path.exists(as_filepath):
                f = open(grid_filepath)
                path = Path(grid_filepath)
                with path.open(mode="rb") as f:
                    self.rdock_grid = File(f, name=path.name)
                    self.save()

                path = Path(as_filepath)
                with path.open(mode="rb") as f:
                    self.rdock_as = File(f, name=path.name)
                    self.save()

        if self.rdock_grid and not os.path.exists(grid_filepath):
            content = self.rdock_grid.read()
            with open(grid_filepath, "wb") as fw:
                fw.write(content)
        if self.rdock_as and not os.path.exists(as_filepath):
            content = self.rdock_as.read()
            with open(as_filepath, "wb") as fw:
                fw.write(content)

    #####################
    #### ESP METHODS ####
    #####################

    def _path_to_ligand_esp_directory(self):
        """
        Create a pdb file in a subdirectory for ESP-DNN to use if it doesn't exist
        :return: path to esp directory
        """
        pdb_file = series_files_path(
            self, f"espdnn/ligand/{self.id}_ligand_espdnn.pdb", local=True
        )
        if not os.path.exists(pdb_file):
            mols = [
                mol
                for mol in Chem.SDMolSupplier(
                    self.default_file(path=True), removeHs=False
                )
            ]
            Chem.MolToPDBFile(mols[0], pdb_file)
        return os.path.dirname(pdb_file)

    def _path_to_receptor_esp_directory(self):
        """
        Create a pdb file in a subdirectory for ESP-DNN to use if it doesn't exist
        :return: path to esp directory
        """
        pdb_file = series_files_path(
            self, f"espdnn/receptor/{self.id}_receptor_espdnn.pdb", local=True
        )
        if not os.path.exists(pdb_file):
            with open(pdb_file, "w") as fw:
                # Sometimes pdbs from Maestro don't work with the ESP DNN, use rdkit to convert
                # the saved pdb file to rdkit's format before saving it to the esp_dnn directory
                new_pdb = self.pdb_block()
                fw.write(new_pdb)

        return os.path.dirname(pdb_file)

    def generate_ligand_esp_map(self):
        """
        Generate esp map file for the series ligand
        :return: pqr file as a text string if the map was generated
        """
        pdb_dir = self._path_to_ligand_esp_directory()
        if run_esp_predict(pdb_dir):
            espdnn_filepath = f"{pdb_dir}/{self.id}_ligand_espdnn.pdb.pqr"
            with open(espdnn_filepath) as f:
                pqr_text = f.read()
            dx_filepath = run_apbs(espdnn_filepath)
            dx_text = ""
            if os.path.exists(dx_filepath):
                with open(dx_filepath) as f:
                    dx_text = f.read()
            return {"pqr": pqr_text, "dx": dx_text}

        return {"pqr": "", "dx": ""}

    def generate_receptor_esp_map(self):
        """
        Generate esp map file for the series receptor
        :return: pqr file as a text string if the map was generated
        """
        pdb_dir = self._path_to_receptor_esp_directory()
        if run_esp_predict(pdb_dir, "protein"):
            with open(f"{pdb_dir}/{self.id}_receptor_espdnn.pdb.pqr") as f:
                esp_text = f.read()
                return esp_text

        return ""


class CompoundQuerySet(models.QuerySet):
    def has_substruct(self, smiles):
        """
        Returns QuerySet with Compounds where the smiles is a substructure of the Compound,
        stereo is not taken into consideration; cannot be chained with other filters
        """
        return self.raw(
            "select * from public.main_compound where rdkit_mol@>%s", params=[smiles]
        )

    def is_substruct(self, smiles):
        """
        Returns QuerySet with Compounds where the Compound is a substructure of the given smiles,
        stereo is not taken into consideration; cannot be chained with other filters
        """
        return self.raw(
            "select * from public.main_compound where rdkit_mol<@%s", params=[smiles]
        )

    def is_equal(self, smiles):
        """
        Returns QuerySet with Compounds that are identical to the given smiles,
        stereo is not taken into consideration; cannot be chained with other filters
        """
        return self.raw(
            "select * from public.main_compound where rdkit_mol@=%s", params=[smiles]
        )


class Compound(models.Model):
    """
    Class describing a chemical compound:
    - project: FK to Project
    - series: FK to Series, can be null

    - virtual_id: virtual generated upon first time compound is seen
    - inchi: rdkit inchi generated from ROMol object w/ flag '/suu'
    - smiles: SMILES representation of compound
    - rdkit_mol: Postgres RDKit Mol representation

    - dn_id: populated when compound has a DN number in DTX
    - modified: datetime when parent compound was last modified
    - measured_data: data pulled from experiments in Dotmatics
    - metadata: additional information related to this compound
    - mmps: ManyToMany field to other Compounds that are classified as a matched pair
    - is_mmpdb_analogue: a boolean, True if this Compound was generated by mmpdb and has never
        appeared in Dotmatics or been used by a real user

    - sdf_file: location of generic sdf file
    - ligprep_file: location of ligprepped structure file
    - confgen_file: location of conformer file
    """

    project = models.ForeignKey("main.Project", on_delete=models.PROTECT, null=True)
    series = models.ForeignKey(Series, on_delete=models.PROTECT, null=True)

    virtual_id = models.CharField(max_length=15, blank=True)
    inchi = models.CharField(max_length=700, blank=True)
    smiles = models.CharField(max_length=500, blank=True)
    rdkit_mol = MolField(null=True, blank=True)

    # `dn_id` can be up to 3 dn ids, separated by a comma
    dn_id = models.CharField(max_length=29, blank=True)
    modified = models.DateTimeField(auto_now_add=True)
    measured_data = models.JSONField(default=dict, blank=True)
    metadata = models.JSONField(default=dict, blank=True)
    mmps = models.ManyToManyField("main.Compound", blank=True)
    is_mmpdb_analogue = models.BooleanField(default=False)

    svg_file = models.FileField(
        blank=True, storage=select_media_storage, upload_to=compound_files_path
    )
    sdf_file = models.FileField(
        blank=True, storage=select_media_storage, upload_to=compound_files_path
    )
    ligprep_file = models.FileField(
        blank=True, storage=select_media_storage, upload_to=compound_files_path
    )
    confgen_file = models.FileField(
        blank=True, storage=select_media_storage, upload_to=compound_files_path
    )

    objects = CompoundQuerySet.as_manager()

    class Meta:
        # Provide index to enable performing SSS
        indexes = [GistIndex(fields=["rdkit_mol"])]

    def _update_rdkit_mol(self):
        """
        Custom SQL to update the RDKit Mol DB field (called in save); Compound does not have to be saved again after calling this function
        """
        update_molfield_sql = f'UPDATE public.main_compound SET "rdkit_mol"=mol_from_smiles(smiles::cstring) WHERE "id"=%s RETURNING *;'
        with connection.cursor() as cursor:
            cursor.execute(update_molfield_sql, [self.pk])

    def save(self, *args, **kwargs):
        """
        Set values for ID fields
        """
        self.modified = now()
        # When a new compound is seen, populate the 'virtual_id' field based on its 'id'
        is_initial_save = not self.id
        # Must save first to get self.id value
        super().save(*args, **kwargs)
        if is_initial_save:
            self.virtual_id = f"VSM-{self.id:011}"
            if not self.is_mmpdb_analogue and not self.series:
                self.series = self.pick_series()
            self.save()

        self._update_rdkit_mol()

    @property
    def name(self):
        if not self.dn_id:
            return self.virtual_id
        return self.dn_id

    def update_dn_id(self, dn_id):
        """
        Given a DN ID, updates the `dn_id` field on this compound and updates the names of all structures
        in the Compound's file fields
        :param dn_id: a string, the new DN ID
        """
        self.dn_id = dn_id
        if dn_id:
            self.is_mmpdb_analogue = False
        self.save()
        if not self.dn_id:
            # Don't bother updating the files to the VSM if there is no DN ID
            return
        file_fields = [
            field.name
            for field in self._meta.get_fields()
            if isinstance(field, models.FileField)
            and ".sdf" in getattr(self, field.name).name
        ]
        for file_field in file_fields:
            file = getattr(self, file_field)
            if not file:
                continue
            filename = os.path.split(file.name)[-1]
            old_filepath = compound_files_path(self, filename, local=True)
            if not os.path.exists(old_filepath):
                # Copy file from S3
                with open(old_filepath, "wb") as fw:
                    fw.write(file.read())

            suppl = Chem.SDMolSupplier(old_filepath, removeHs=False)
            writer = Chem.SDWriter(filename)
            for mol in suppl:
                mol.SetProp("_Name", self.dn_id)
                writer.write(mol)
            writer.close()

            # Delete the old file before saving the new file to avoid a filename clash
            file.delete(save=False)
            with open(filename, "rb") as file:
                setattr(self, file_field, File(file, name=filename))
                self.save()

            # Delete the temporary file
            os.remove(filename)

        # Regenerate the SDF file with the Dotmatics orientation
        dn_mol = Chem.MolFromMolBlock(get_registered_structures([dn_id])[dn_id])
        self.get_sdf_file(dn_mol)

    ########################
    ## FORMATTING METHODS ##
    ########################

    def mol(self, hydrogens=False):
        """
        Return this compound as an RDKit mol object from the compound's SDF file (not 3D)
        :param hydrogens: a boolean, should hydrogens be added explicitly
        :returns: an RDKit mol object
        """
        sdf_path, _ = self.get_sdf_file()
        mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
        mol.SetProp("_Name", self.name)
        if hydrogens:
            mol = Chem.AddHs(mol, addCoords=True)
        return mol

    def generate_svg_file(self):
        """
        Generates and saves an svg file to the model.
        """
        filename = f"{self.pk}.svg"
        localpath = compound_files_path(self, filename, local=True)
        if not self.svg_file:
            # Generate the file if it does not exist
            if not os.path.exists(localpath):
                with open(localpath, "w") as f:
                    f.write(self._svg())
            # Save the file to the model
            with open(localpath, "rb") as data:
                if isinstance(select_media_storage(), FileSystemStorage):
                    # If media storage is local storage, delete it so there is no filename clash
                    os.remove(localpath)
                self.svg_file.save(filename, data)

    def get_sdf_file(self, mol=None):
        """
        Generates and saves sdf_file to model
        :param mol: optional, a 2D RDkit mol object to write to the SDF file. If omitted or 3D, a mol
            will be generated from the Compound's SMILES. This argument is used to maintain
            user-uploaded or Dotmatics orientation for depictions
        :return: path to file (including name), filename for the sdf
        """
        filename = f"{self.pk}.sdf"
        localpath = compound_files_path(self, filename, local=True)

        # If a mol has been specified, delete the existing sdf file so that regenerates w/ the given molblock
        if mol:
            if RDKitWrappers.is_3d(mol):
                # 3D mols should not be saved as the Compound sdf
                mol = None
            else:
                if os.path.exists(localpath):
                    os.remove(localpath)
                self.sdf_file.delete()
                mol.SetProp("_Name", self.name)

        if not self.sdf_file:
            if not os.path.exists(localpath):
                if not mol:
                    # Use smiles to generate the mol since there was no mol supplied
                    mol_from_smiles = Chem.MolFromSmiles(self.smiles)
                    mol_from_smiles.SetProp("_Name", self.name)
                    return self.get_sdf_file(mol_from_smiles)
                writer = Chem.SDWriter(localpath)
                writer.write(mol)
                writer.close()

            with open(localpath, "rb") as data:
                if isinstance(select_media_storage(), FileSystemStorage):
                    # If media storage is local storage, delete it so there is no filename clash
                    os.remove(localpath)
                self.sdf_file.save(filename, data)

        elif self.sdf_file and not os.path.exists(localpath):
            content = self.sdf_file.read()
            if content:
                with open(localpath, "wb") as fw:
                    fw.write(content)
            else:
                # If the file is empty/gone missing try to regenerate it
                self.sdf_file.delete()
                self.get_sdf_file()

        return localpath, filename

    def _svg(self, size_x=450, size_y=150, transparent=False):
        """
        Draw this Compound as an SVG in the given canvas size
        :param size_x: the width of the svg in pixels
        :param size_y: the height of the svg in pixels
        :param transparent: a boolean, should the background be transparent? If False, the background is white
        :return: a string with SVG data
        """
        moltext = Chem.MolToMolBlock(self.mol())
        return moltext_to_svg(moltext, size_x, size_y, transparent)

    def get_inline_svg(self, size_x=450, size_y=150, transparent=False):
        """
        Creates a safe html rendering of the compound depiction as svg text
        :param size_x: the width of the svg in pixels
        :param size_y: the height of the svg in pixels
        :param transparent: a boolean, should the background be transparent? If False, the background is white
        :return: safe html svg tag
        """
        svg_prefix = re.compile(r"<\?xml .*\?>\s*<svg")
        result = re.sub(svg_prefix, "<svg", self._svg(size_x, size_y, transparent))
        result = mark_safe(result)
        return result

    #######################
    ### GENERAL METHODS ###
    #######################

    def pick_series(self):
        """
        Select the closest series match for this compound
        :return: series object or None
        """
        series_qs = Series.objects.filter(project=self.project, active=True).order_by(
            "dn_id"
        )
        if series_qs.exists():
            selected_series = RDKitWrappers.pick_series(
                self.mol(), [(s.mol(), s.dn_id) for s in series_qs]
            )
            if selected_series:
                return series_qs.get(dn_id=selected_series)
            else:
                logger.error(
                    f"Could not select a series for this compound: {self.virtual_id}"
                )

    def get_ligprep_path(self):
        """
        If ligprep has not been run, start job and return/save path to compound
        :return: filepath to ligprep file or None if file was not generated
        """
        filename = f"{self.pk}_ligprep.sdf"
        localpath = compound_files_path(self, filename, local=True)

        if not self.ligprep_file:
            if not os.path.exists(localpath):
                # TODO: GypsumDL is not a good alternative for running ligprep so for now we are just
                # taking the low energy conformation as a starting point
                mol = Chem.SDMolSupplier(self.get_confgen_path(), removeHs=False)[0]
                w = Chem.SDWriter(localpath)
                w.write(mol)
                w.close()
            with open(localpath, "rb") as data:
                if isinstance(select_media_storage(), FileSystemStorage):
                    # If media storage is local storage, delete it so there is no filename clash
                    os.remove(localpath)
                self.ligprep_file.save(filename, data)

        elif self.ligprep_file and not os.path.exists(localpath):
            content = self.ligprep_file.read()
            with open(localpath, "wb") as fw:
                fw.write(content)

        return localpath

    def get_first_ligprep_path(self):
        """
        Returns the first structure from ligprep in a new file
        """
        filename = f"{self.pk}_first_ligprep.sdf"
        localpath = compound_files_path(self, filename, local=True)

        if not os.path.exists(localpath):
            ligprep_file = self.get_ligprep_path()
            suppl = [m for m in Chem.SDMolSupplier(ligprep_file, removeHs=False)]
            with Chem.SDWriter(localpath) as w:
                w.write(suppl[0])

        return localpath

    def get_confgen_path(self):
        """
        If confgen has not been run, start job and return/save path to compound
        :return: filepath to file with conformers or None if file was not generated
        """
        filename = f"{self.pk}_conformers.sdf"
        localpath = compound_files_path(self, filename, local=True)

        if not self.confgen_file:
            if not os.path.exists(localpath):
                _ = RDKitWrappers.confgen(self.mol(), localpath)

            with open(localpath, "rb") as data:
                if isinstance(select_media_storage(), FileSystemStorage):
                    # If media storage is local storage, delete it so there is no filename clash
                    os.remove(localpath)
                self.confgen_file.save(filename, data)

        elif self.confgen_file and not os.path.exists(localpath):
            content = self.confgen_file.read()
            with open(localpath, "wb") as fw:
                fw.write(content)

        return localpath

    def get_lsalign_path(self, reference=None):
        """
        Returns the filepath to the lowest energy flexibly aligned structure to the given reference
        :param reference: parent object to superimpose to; Compound, Series or
            None which defaults to this Compounds' Series
        :return: filepath to lowest energy flexibly aligned file or None if file was not generated
        """
        ref_file, ref_id = self._pick_reference_file_and_id(reference)

        filename = f"{self.pk}_lsaligned_to_{ref_id}.sdf"
        localpath = compound_files_path(self, filename, local=True)

        if not os.path.exists(localpath):
            localpath = run_lsalign(self.get_confgen_path(), ref_file, localpath)

        return localpath

    ######################
    ## PROPCALC METHODS ##
    ######################

    def predict_logd(self):
        """
        Runs the InductiveBio logD model for this compound
        :return: {"ilogd": str, "measured": str, "latest_data_date": str}
        """
        sdf_path, _ = self.get_sdf_file()

        ib_logd_response = run_inductive_alogd_predict(sdf_path)[0]
        ib_logd = {}
        ib_logd["ilogd"] = ib_logd_response["prediction"]
        ib_logd["measured"] = ib_logd_response["measured"]
        ib_logd["latest_data_date"] = ib_logd_response["latest_data_date"]

        return ib_logd

    def predict_lms(self, species, image=True):
        """
        Generates predicted values for LMs from the Inductive API
        :param species: H or R for human or rat LM predictions
        :param image: boolean for if interpretation images should be generated
        :return: dictionary with "prediction, "measured", and "out_of_domain" values and images
        """
        infile, _ = self.get_sdf_file()
        # Only passed one compound, so use only element in list
        lm_prediction = run_inductive_lm_predict(infile, species, image)[0]

        # Save base64 images to files
        if image:
            lm_prediction["interp_img_path"] = self.inductive_img_path(
                species, "interp", lm_prediction["interp_image"]
            )
            lm_prediction["probs_img_path"] = self.inductive_img_path(
                species, "probs", lm_prediction["probs_image"]
            )

        return lm_prediction

    def inductive_img_path(self, species, name, img):
        """
        Saves given file and returns media path to the image
        :param species: the species of LMs being predicted
        :param name: the name of the image type ("interp" or "probs")
        :param img: the base64 image to save
        """
        filename = f"{self.name}_{species.lower()}lm_{name}_img.png"
        filepath = compound_files_path(self, filename)

        img = base64.decodebytes(str.encode(img))
        tmp_file = get_tmp_file(img, filename)
        with open(tmp_file, "rb") as fp:
            media_path = save_media_file(filepath, fp)

        return media_path

    def mol_w_properties(self, inductive=True, dn_id=True, hydrogens=False):
        """
        Helper to return rdkit calculated properties
        :param inductive: Bool for including InductiveBio based properties
        :param dn_id: Bool for including the DN ID
        :param hydrogens: a boolean, should hydrogens be added explicitly
        :return: RoMol that contains all properties
        """
        mol = self.mol(hydrogens=hydrogens)
        mol_w_props = RDKitWrappers.generate_properties(mol)
        if inductive:
            # Get LogD predictions
            ilogd_predict = self.predict_logd()
            mol_w_props.SetProp("logd_prediction", ilogd_predict["ilogd"])
            mol_w_props.SetProp("logd_measured", ilogd_predict["measured"])
            mol_w_props.SetProp(
                "latest_logd_data_date", ilogd_predict["latest_data_date"]
            )
            # Get Rat LM predictions
            rlm_predict = self.predict_lms(species="R")
            mol_w_props.SetProp("rlm_prediction", rlm_predict["prediction"])
            mol_w_props.SetProp("rlm_measured", rlm_predict["measured"])
            mol_w_props.SetProp("rlm_probabilities", rlm_predict["probs_img_path"])
            mol_w_props.SetProp("rlm_interpretation", rlm_predict["interp_img_path"])
            mol_w_props.SetProp("rlm_ood", rlm_predict["out_of_domain"])
            mol_w_props.SetProp("latest_lm_data_date", rlm_predict["latest_data_date"])
            mol_w_props.SetProp("model_version", rlm_predict["model_version"])
            # Get Human LM predictions
            hlm_predict = self.predict_lms(species="H")
            mol_w_props.SetProp("hlm_prediction", hlm_predict["prediction"])
            mol_w_props.SetProp("hlm_measured", hlm_predict["measured"])
            mol_w_props.SetProp("hlm_probabilities", hlm_predict["probs_img_path"])
            mol_w_props.SetProp("hlm_interpretation", hlm_predict["interp_img_path"])
            mol_w_props.SetProp("hlm_ood", hlm_predict["out_of_domain"])
            mol_w_props.SetProp("latest_lm_data_date", hlm_predict["latest_data_date"])
            mol_w_props.SetProp("model_version", hlm_predict["model_version"])

        if dn_id and self.dn_id:
            mol_w_props.SetProp("dn_id", self.dn_id)

        return mol_w_props

    def run_propcalc(self, inductive):
        """
        Async tasks that gets launched for the given compound
        Currently needed to include svg result in final result
        :param inductive: bool, should InductiveBio props be included
        :return: dictionary of properties
        """
        mol = self.mol_w_properties(inductive=inductive)
        return mol.GetPropsAsDict()

    #####################
    ### ALIGN METHODS ###
    #####################

    def _get_superimposed_path(self, reference):
        """
        Wrapper to run superimpose to reference with the given conformer file
        :param reference: parent object to superimpose to; Compound, Series or
            None which defaults to this Compounds' Series
        :return: a path to the file with superimposed conformers
        """
        _, ref_id = self._pick_reference_file_and_id(reference)
        lsalign_path = self.get_lsalign_path(reference)

        superimpose_output_filename = f"{self.id}_superimposed_to_{ref_id}.sdf"
        localpath = compound_files_path(self, superimpose_output_filename, local=True)

        if not os.path.exists(localpath):
            localpath = RDKitWrappers.superimpose(
                self.get_confgen_path(), lsalign_path, localpath
            )

        return localpath

    def _select_best_confs(self, confs):
        """
        Helper function to select the 4 lowest rmsd conformers and the 4
        lowest energy conformers. If those confs overlap, there will be
        less than 8 confs shown, no replacements made.
        :param confs: list of conformers sorted by relative_energy, ascending
        :return: Ordered dictionary of the top MAX_N_CONFS_DISPLAY conformers in the form
            {conf_id: {moltext: mol, r_mmff_rel_energy: float, r_bc_rmsd_to_lsalign:float }}
        """
        sorted_by_rmsd_confs = sorted(
            confs, key=lambda x: float(x.GetProp("r_bc_rmsd_to_lsalign"))
        )
        confs_w_energy = self._calculate_relative_energies(sorted_by_rmsd_confs)
        sorted_by_energy_confs = sorted(
            confs_w_energy, key=lambda x: float(x.GetProp("r_mmff_rel_energy"))
        )

        confs = (
            sorted_by_rmsd_confs[: int(MAX_N_CONFS_DISPLAY / 2)]
            + sorted_by_energy_confs[: int(MAX_N_CONFS_DISPLAY / 2)]
        )

        confs_dict = OrderedDict()
        for conf in confs:
            conf_id = conf.GetProp("s_conf_number")
            alerts, energy = mol_to_torsion_alerts(conf)
            confs_dict[conf_id] = {
                "moltext": Chem.MolToMolBlock(conf),
                "r_mmff_rel_energy": "%.3f" % float(conf.GetProp("r_mmff_rel_energy")),
                "r_bc_rmsd_to_lsalign": "%.3f"
                % float(conf.GetProp("r_bc_rmsd_to_lsalign")),
                "torsion_alerts": alerts,
                "torsion_alerts_total_energy": energy,
            }

        return confs_dict

    def _pick_reference_file_and_id(self, reference):
        """
        Selects which structure file and id should be used for each type of reference
        :param reference: a Compound, Series, or None
        :return: filepath to selected reference structure, ref_id
        """
        reference_file = None
        if not reference:
            reference = self.series

        if type(reference) == type(self):
            reference_file = reference.get_ligprep_path()
            ref_id = f"c-{reference.pk}"
        elif type(reference) == Series:
            reference_file = reference.default_file(path=True)
            ref_id = f"s-{reference.pk}"

        return reference_file, ref_id

    def _calculate_relative_energies(self, confs):
        """
        Calculate relative energy of each conformer to lowest energy in each
        set of conformers
        :param confs: list of conformers
        :return: list of conformers with energy values as a property
        """
        if "isTemplate" in confs[0].GetPropNames():
            # Add energy for just ls-aligned conformer since the other confs already have energy calculated
            lsaligned = confs[0]
            mp = AllChem.MMFFGetMoleculeProperties(lsaligned)
            ffm = AllChem.MMFFGetMoleculeForceField(lsaligned, mp)
            energy_value_M = ffm.CalcEnergy()
            lsaligned.SetProp("r_mmff_energy", str(energy_value_M))

        min_energy_in_set = min([float(c.GetProp("r_mmff_energy")) for c in confs])
        for conf in confs:
            relative_energy = float(conf.GetProp("r_mmff_energy")) - float(
                min_energy_in_set
            )
            conf.SetProp("r_mmff_rel_energy", str(relative_energy))

        return confs

    #####################
    #### MMP METHODS ####
    #####################

    def mmp_analysis(self, constant_smiles=None, variable_smiles=None):
        """
        Run MMP analysis on this compound and return a list of Compound objects that are MMPs
        for the given constant and variable regions
        :param constant_smiles: a string, the SMILES string of the user-defined constant region
        :param variable_smiles: a string, the SMILES string of the user-defined variable region
        :return: a list of PKs of Compound objects that are MMPs for this Compound for the given constant and variable regions
        """
        df = self.find_mmps()
        if constant_smiles:
            df["constant"] = df["constant"].apply(
                lambda constant: Chem.MolToSmiles(Chem.MolFromSmiles(constant))
            )
            df = df.loc[df["constant"] == constant_smiles]
        if variable_smiles:
            df["variable"] = df["from_smiles"].apply(
                lambda variable: Chem.MolToSmiles(
                    Chem.MolFromSmiles(re.sub(r"\[\*\:\d\]", "*", variable))
                )
            )
            df = df.loc[df["variable"] == variable_smiles]

        # Filter the mmps to only include Compound objects that are present in the filtered dataframe
        comp_pks = list(
            set([comp.pk for comp_list in df["vsms"].tolist() for comp in comp_list])
        )
        return list(self.mmps.filter(pk__in=comp_pks).values_list("pk", flat=True))

    def find_mmps(self, radius=0):
        """
        Find mmps and save them to the Compound's mmps field
        :param radius: the radius argument value for generate
        """
        generate_path = compound_files_path(self, "generate.csv", local=True)
        df = generate_mmps(self.smiles, generate_path, radius=radius)
        if df.empty:
            return

        df["vsms"] = df.apply(
            lambda x: list(Compound.objects.is_equal(x["final"])), axis=1
        )
        df["keep_row"] = df["meets_mmp_prop_req"] | df["vsms"]
        df = df.loc[df.keep_row, :]

        c_fp = Chem.RDKFingerprint(self.mol())
        for i, row in df.iterrows():
            comps = row["vsms"]
            mol = row["final_mol"]
            if not comps:
                try:
                    comps = [
                        Compound.objects.get(inchi=Chem.MolToInchi(mol, options="/suu"))
                    ]
                except ObjectDoesNotExist:
                    comps = [
                        Compound.objects.create(
                            inchi=Chem.MolToInchi(mol, options="/suu"),
                            smiles=Chem.MolToSmiles(mol),
                            is_mmpdb_analogue=True,
                            project=self.project,
                        )
                    ]
                df.loc[i, "vsms"] = comps
            for comp in comps:
                sim = DataStructs.TanimotoSimilarity(c_fp, Chem.RDKFingerprint(mol))
                comp.metadata.setdefault("similarity", {})
                comp.metadata["similarity"][self.pk] = sim
                comp.save()
            self.mmps.add(*comps)
        self.mmps.filter(project__code=AUTO).update(project=self.project)
        self.mmps.remove(self)
        return df

    def update_mmp_dtx_avg_assay_data(self, dtx_data=None, dns_to_skip=None):
        """
        Update the average assay values for this Compound and all of its DTX MMPs
        :param dtx_data: a dictionary of the form {dn_id: {assay_data}} with aggregated assay
            results from Dotmatics. If None, this function will make a call to Dotmatics.
        :param dns_to_skip: a list of DN IDs (strings) to skip. This is used for the newly assayed Compounds
            because brand-new assay data is not yet in the aggregate datasource, so we show the most
            recent assay value instead of aggregate data.
        """
        dns_to_skip = dns_to_skip or []
        dtx_mmps = self.mmps.exclude(dn_id="")
        if not dtx_data:
            dn_ids = list(dtx_mmps.values_list("dn_id", flat=True))
            if self.dn_id:
                dn_ids.append(self.dn_id)
            dtx_data = get_agg_ic50_data(dn_ids)
        if self.dn_id not in dns_to_skip:
            self.update_dtx_avg_assay_data(dtx_data)
        for mmp in dtx_mmps:
            if self.dn_id not in dns_to_skip:
                mmp.update_dtx_avg_assay_data(dtx_data)

    def update_dtx_avg_assay_data(self, dtx_data=None):
        """
        Update the average assay values for this Compound
        :param dtx_data: a dictionary of the form {dn_id: {assay_data}} with aggregated assay
            results from Dotmatics. If None, this function will make a call to Dotmatics.
        """
        if not self.dn_id:
            return  # Not in DTX, no assay data to pull
        if not dtx_data:
            dtx_data = get_agg_ic50_data([self.dn_id])
        assay_results = self.measured_data.get("assay_results", {})
        dtx_data = dtx_data.get(self.dn_id, {}).get("data", {})
        for assay_name, assay_data in dtx_data.items():
            if assay_name not in assay_results:
                assay_results[assay_name] = {}
            for analysis_name, analysis_data in assay_data.items():
                mean_result = analysis_data.get("RESULT_GEOM_MEAN")
                if mean_result:
                    assay_results[assay_name][analysis_name] = mean_result
        self.measured_data["assay_results"] = assay_results
        self.save()


class ParentsManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset().filter(parent_co__isnull=True)


class GemsManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset().filter(parent_co__isnull=False)


class CompoundOccurrence(models.Model):
    """
    Class containing information about the usage of a compound
    compound: FK to parent compound object
    owner: FK to user who created this occurrence
    generated: datetime when occurrence was generated
    molblock: molblock of specific 3D conformer (user uploaded)
    parent_co: CompoundOccurrence this CO is derived from
    gem_id: ID to use if this CO is saved as a Gem
    """

    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)
    owner = models.ForeignKey(BasechemUser, on_delete=models.PROTECT)
    generated = models.DateTimeField(auto_now_add=True)
    molblock = models.TextField(max_length=6000, blank=True)
    parent_co = models.ForeignKey(
        "main.CompoundOccurrence", on_delete=models.CASCADE, null=True, blank=True
    )
    saved_from = models.CharField(
        max_length=20, choices=SAVED_FROM_OPTIONS, default=USER_UPLOADED
    )
    gem_id = models.CharField(max_length=30, blank=True)

    objects = models.Manager()
    parents = ParentsManager()
    gems = GemsManager()

    def clean(self):
        """
        This is called as part of forms that create or edit a CompoundOccurrence instance.
        This function is called in the Django admin when editing CompoundOccurrences, which helps
        our maintainers avoid mistakes.
        """
        if self.parent_co:
            if self.compound != self.parent_co.compound:
                raise ValidationError(
                    {"compound": "must be the same as the parent_co's compound"}
                )
            if self.owner != self.parent_co.owner:
                raise ValidationError(
                    {"owner": "must be the same as the parent_co's owner"}
                )

    def display_gem_with_series(self):
        """
        If the object has a gem_id, return a user-readable name that includes the series
        :returns: a string (ex: "c10-co590-1 (series_name)" when gem_id is "c10-co590-s-2")
        """
        if not self.gem_id:
            return ""
        matches = re.match("(.*)-s-(\d+)", self.gem_id)
        if matches:
            conf_id = matches[1]
            series_name = Series.objects.get(pk=int(matches[2])).name
            return f"{conf_id} ({series_name})"
        return self.gem_id

    def display_gem_without_series(self):
        """
        If the object has a gem_id, return a user-readable name that does not include the series
        :returns: a string (ex: "c10-co590-1" when gem_id is "c10-co590-s-2")
        """
        if not self.gem_id:
            return ""
        return strip_series_from_conf_id(self.gem_id)

    def is_user_uploaded(self):
        return self.parent_co == None

    def is_2d(self):
        return self.molblock == ""

    def get_child_cos(self, collection=None):
        """
        :param collection: optional, a Collection object. If provided, only child
            objects that are in the provided collection will be returned
        :return: a queryset of CompoundOccurrence objects that were basechem-generated
            from this CompoundOccurrence
        """
        queryset = CompoundOccurrence.objects.filter(parent_co=self)
        if collection:
            queryset = queryset.filter(collection=collection)
        return queryset

    @property
    def mol(self):
        """
        Creates a mol object from the compound occurrence's saved `molblock`. If `molblock` is None,
        the compound occurrence has no particular conformation, so the compound's mol is returned
        (generated by smiles string)
        :returns: an RDKit mol object
        """
        if self.molblock:
            # If this occurrence has a 3D molblock, use it
            return Chem.MolFromMolBlock(self.molblock, removeHs=False)
        else:
            # No 3D molblock, use the 2D mol from the compound
            return self.compound.mol()

    def moltext(self, twoD=False):
        """
        Used in the torsion front end, this function returns moltext for the compound occurrence.
        If the object does not have a particular conformation (molblock != None), moltext is generated
        from the smiles string of the compound
        :param twoD: a boolean, if True returns 2D moltext even if the instance has a 3D conformation
        :returns: a string of moltext
        """
        if self.molblock and not twoD:
            return self.molblock
        else:
            return Chem.MolToMolBlock(self.compound.mol())

    def get_sdf_file(self):
        """
        Generates and saves an SDF file for this CompoundOccurrence. If using the Compound's
        sdf file, the file is not copied (it remains in the Compound's directory)
        :return: path to file (including filename), filename
        """
        if self.is_2d():
            return self.compound.get_sdf_file()
        else:
            filename = f"{self.pk}.sdf"
            localpath = compound_files_path(
                self.compound, filename, co=self, local=True
            )

            if not os.path.exists(localpath):
                writer = Chem.SDWriter(localpath)
                writer.write(self.mol)
                writer.close()

            return localpath, filename

    def get_3d_sdf_file(self, reference=None):
        """
        Returns the filepath and filename of an SDF file for this CompoundOccurrence with 3D.
        For 3D COs, this sdf file has the molblock stored in `self.molblock`. For 2D COs, this
        returns the path to the relevant flexibly-aligned SDF file.
        :param reference: for 2D COs, the parent object to superimpose to (Compound/Series/None).
            If None, the Compound's default series is used.
        :return: filepath (including filename), filename
        """
        if self.is_2d():
            lsalign_path = self.compound.get_lsalign_path(reference)
            _, filename = os.path.split(lsalign_path)
            return lsalign_path, filename
        else:
            return self.get_sdf_file()

    def get_propcalc_task_name(self, collection):
        """
        :param collection: the Collection that is generating this task
        :return: the name of the propcalc async task for this CompoundOccurrence
        """
        return f"propcalc_{self.pk}_{collection.pk}"

    #####################
    ### ALIGN METHODS ###
    #####################

    def get_align_task_name(self, ref_string):
        """
        :param ref_string: a string:
            - 'default': align all compounds to their assigned series
            - 's-(int)': align all compounds to the series w/ id (int)
        :return: the name of the async task for alignment with the ref_string
        """
        return f"{ALIGN}_{self.pk}_{ref_string}"

    def superimpose_to_ref(self, reference=None, torsion_alerts=True):
        """
        Returns a dictionary with the moltext of the conformers with the lowest
        rmsd and lowest energy aligned to the given reference compound.
        Default to return 8 (MAX_N_CONFS_DISPLAY) but if the lowest rmsd and lowest
        energy confs overlap, there will be less than 8
        :param reference: parent object to superimpose to; Compound, Series or
             None which defaults to this Compounds' Series
        :param torsion_alerts: a boolean, if True calculate torsion alerts for all returned conformers
        :return: a nested dictionary of the form
            {conf-id:  {"moltext": conf_moltext, "r_mmff_rel_energy": float, "r_bc_rmsd_to_lsalign":float}
        """
        _, ref_id = self.compound._pick_reference_file_and_id(reference)
        filepath = self.compound._get_superimposed_path(reference)

        if torsion_alerts:
            filepath = generate_torsion_alerts(filepath)
        suppl = Chem.SDMolSupplier(filepath, removeHs=False)

        confs = [mol for mol in suppl]
        for i, conf in enumerate(confs, start=1):
            conf.SetProp(
                "s_conf_number", f"c{self.compound.pk}-co{self.pk}-{i}-{ref_id}"
            )

        confs_dict = self.compound._select_best_confs(confs)

        return confs_dict

    ######################
    #### DOCK METHODS ####
    ######################

    def get_dock_task_name(self, ref_string):
        """
        :param ref_string: a string:
            - 'default': dock all compounds to the receptor of their assigned series
            - 's-(int)': dock all compounds to the receptor of the series w/ id (int)
        :return: the name of the async task for docking with the ref_string
        """
        return f"{DOCK}_{self.pk}_{ref_string}"

    def dock_to_receptor(self, reference=None, torsion_alerts=True, toklat=True):
        """
        Returns a dictionary with the moltext of the poses returned from running docking
        :param reference: a Series whose bound ligand and receptor should be used as a reference for docking this compound
        :param torsion_alerts: a boolean, if True calculate torsion alerts for all returned conformers
        :param toklat: a boolean, use the Toklat pose scoring model to score and annotate poses
        :return: a dictionary of the form {pose-id: {"moltext": conf_moltext, "dockingScore": score, "RMSDtoLSAligned": value}
        """
        if not reference:
            reference = self.compound.series

        # refresh from db incase the reference just had its rdock files generated
        reference.refresh_from_db()
        reference_prm_exists = os.path.exists(
            reference.path_to_local_file(reference.rdock_prm)
        )
        reference_mol2_exists = os.path.exists(
            reference.path_to_local_file(reference.receptor_file_mol2)
        )

        if reference_prm_exists and reference_mol2_exists:
            filepath = self._get_rdock_output_path(reference)

            if torsion_alerts and os.path.exists(filepath):
                torsion_filepath = generate_torsion_alerts(filepath)
                if torsion_filepath and os.stat(torsion_filepath).st_size != 0:
                    filepath = torsion_filepath
            if toklat and os.path.exists(filepath):
                filepath = RDKitWrappers.add_hydrogens_to_sdf(filepath)
                filepath = generate_toklat_scores(
                    filepath, reference.path_to_local_file(reference.receptor_file)
                )
            if os.path.exists(filepath):
                pose_dict = self._select_rdock_poses(filepath, reference, toklat)
                return pose_dict

        return {}

    def _select_rdock_poses(self, rdock_filepath, reference, toklat=True):
        """
        Helper function to select the 3 best scoring docking poses along with the 2 poses
        most similar to the flexibly aligned conformer for this compound.
        :param rdock_filepath: path to rdock output
        :param reference: the Series used for this docking run
        :param toklat: a boolean, use the Toklat pose scoring model to score and annotate poses
        :return: dictionary of poses of the form {pose-id:  {"moltext": conf_moltext, "dockingScore": score, "RMSDtoLSAligned": value}
        """
        _, ref_id = self.compound._pick_reference_file_and_id(reference)
        lsalign_path = self.compound.get_lsalign_path(reference)
        lsalign_mol = [mol for mol in Chem.SDMolSupplier(lsalign_path, removeHs=False)][
            0
        ]
        receptor_pdb_path = reference.path_to_local_file(reference.receptor_file)

        rdock_mols = [mol for mol in Chem.SDMolSupplier(rdock_filepath, removeHs=False)]
        for i, mol in enumerate(rdock_mols, start=1):
            conf_id = mol.GetProp("s_conf_id")
            pose_id = f"c{self.compound.pk}-co{self.pk}-{conf_id}-{i}-{ref_id}"
            mol.SetProp("s_pose_id", pose_id)
            try:
                rmsd_to_lsa = Chem.rdMolAlign.CalcRMS(mol, lsalign_mol)
            except:
                # If something is off about the conformer from ls_align, manually set rmsd
                rmsd_to_lsa = "100"
            mol.SetProp("r_bc_rdock_rmsd_to_lsalign", "%.3f" % float(rmsd_to_lsa))

        score_prop = "toklat_score" if toklat else "SCORE"
        sorted_by_toklat = sorted(
            rdock_mols, key=lambda x: float(x.GetProp("toklat_score")), reverse=False
        )[:3]
        sorted_by_rdock = sorted(
            rdock_mols, key=lambda x: float(x.GetProp("SCORE")), reverse=False
        )[:3]
        sorted_by_rmsd = sorted(
            rdock_mols,
            key=lambda x: float(x.GetProp("r_bc_rdock_rmsd_to_lsalign")),
            reverse=False,
        )[:2]
        poses = sorted_by_toklat + sorted_by_rdock + sorted_by_rmsd

        pose_dict = OrderedDict()
        for mol in poses:
            pose_id = mol.GetProp("s_pose_id")
            alerts, energy = mol_to_torsion_alerts(mol)
            toklat_annotations = mol_to_toklat_annotations(mol, receptor_pdb_path)
            pose_dict[pose_id] = {
                "moltext": Chem.MolToMolBlock(mol),
                "toklatScore": "%.3f" % float(mol.GetProp("toklat_score")),
                "rdockScore": "%.3f" % float(mol.GetProp("SCORE")),
                "RMSDtoLSAligned": "%.3f"
                % float(mol.GetProp("r_bc_rdock_rmsd_to_lsalign")),
                "torsionAlerts": alerts,
                "torsionAlertsTotalEnergy": energy,
                "toklatAnnotations": toklat_annotations,
            }

        return pose_dict

    def _get_rdock_output_path(self, reference):
        """
        Returns a filepath to the rdock docking output for the given compound and reference
        If the file does not yet exist, runs rdock
        :param reference: a Series whose bound ligand and receptor should be used as a reference for docking this compound
        :return: filepath to docked poses of this compound
        """
        rdock_output_filename = f"{self.id}_{reference.dn_id}_rdock_out.sd"
        localpath = compound_files_path(
            self.compound, rdock_output_filename, co=self, local=True
        )

        if not os.path.exists(localpath):
            os.makedirs(os.path.split(localpath)[0], exist_ok=True)
            self._run_rdock(reference, localpath)

        return localpath

    def _run_rdock(self, reference, output_filepath):
        """
        Sets up the information for a new rDock run for this compound and the given series
        :param reference: the series to dock to
        :param output_filepath: the path to the docking results
        """
        input_confs_path = self.compound._get_superimposed_path(reference)
        receptor_prm = reference.get_rdock_prm()
        ligand_prm = compound_files_path(
            self.compound, "ligand_rdock.prm", co=self, local=True
        )
        copy_rdock_ligand_prm(ligand_prm)

        run_rdock(input_confs_path, output_filepath, receptor_prm, ligand_prm)

    #####################
    #### ESP METHODS ####
    #####################
    def get_esp_task_name(self, ref_string):
        """
        :param ref_string: a string, the reference to use for flexibly aligning if this
            CompoundOccurrence is 2d (ex:"default" or "s-12")
        :returns: a string, the name for this ESP task
        """
        return f"{ESP}_{self.pk}_{ref_string}"

    def generate_esp_map(self, job_name=None, reference=None):
        """
        Returns a dictionary with the pqr of the esp maps generated by ESP-DNN
        :param reference: a Series whose bound ligand and receptor should be used as a reference for
            aligning this compound before generating the ESP map
        :param job_name: the task name for this task (if run asynchronously). This allows
            basechem to use a cached result if this exact task has been run before
        :return: a dictionary of the form {"pqr": conf_pqr_text, "related_series": "s-(pk of series)"}
        """
        # Check if this exact task has already been run successfully
        if job_name:
            cached_task = Task.objects.filter(success=True, name=job_name).first()
            if cached_task and cached_task.result and type(cached_task.result) == dict:
                return cached_task.result
        if not reference:
            reference = self.compound.series

        input_ligand_path, _ = self.get_3d_sdf_file(reference)
        espdnn_directory, expected_espdnn_file = self._get_espdnn_directory(
            input_ligand_path
        )
        if run_esp_predict(espdnn_directory):
            with open(expected_espdnn_file) as f:
                esp_text = f.read()
            dx_filepath = run_apbs(expected_espdnn_file)
            dx_text = ""
            if os.path.exists(dx_filepath):
                with open(dx_filepath) as f:
                    dx_text = f.read()
            if not dx_text:
                raise Exception("APBS failed to generate a dx file.")
            return {
                "pqr": esp_text,
                "related_series": f"s-{reference.pk}",
                "dx": dx_text,
            }

        return {}

    def _get_espdnn_directory(self, input_ligand_path):
        """
        Returns a directory path and expected outfile where ESP predicted maps for this structure can
        be found. ESP-DNN runs for all pdb structures in a directory so each flexibly aligned
        starting structure should be in its own directory
        :param input_ligand_path: path to input ligand (flexibly structure for 2D COs, sdf file for 3D COs)
        :return: path to directory, expected outfile
        """
        outdir = os.path.dirname(input_ligand_path)

        # input filename with .pdb extension -> ex. 75_lsaligned_to_s-2.pdb
        pdb_filename = f"{os.path.splitext(input_ligand_path)[0]}.pdb"
        if not os.path.exists(pdb_filename):
            # There should only be one mol in the input file
            mols = [
                mol for mol in Chem.SDMolSupplier(input_ligand_path, removeHs=False)
            ]
            Chem.MolToPDBFile(mols[0], pdb_filename)

        return outdir, f"{pdb_filename}.pqr"

    #########################
    #### TORSION METHODS ####
    #########################

    def get_torsion_task_name(self, pioneer_dihedral_smarts, dihedral_atoms):
        """
        :param pioneer_dihedral_smarts: a SMARTS string defining the atoms of interest in the pioneer
        :param dihedral_atoms: a comma separated string of 4 atom indices (ex: '3,5,8,9')
        :returns: a string, the name for this torsion task
        """
        return f"{TORSION}_{self.pk}_{pioneer_dihedral_smarts}_{dihedral_atoms}"

    def convert_atoms_to_smarts(self, selected_atoms):
        """
        Given 4 atom indices from this compound occurrence's molblock that define a rotatable bond,
        generate a SMARTS string that can be used to represent the atoms for a torsion scan
        :param selected_atoms: a comma separated string of 4 atom indices (ex: '3,5,8,9')
        :returns: a SMARTS string
        """
        mol = self.mol
        atom_indices = [int(index) for index in selected_atoms.split(",")]

        # Create a submolecule with only the selected atoms
        submol = Chem.RWMol(mol)
        for atom_index in range(mol.GetNumAtoms() - 1, -1, -1):
            if atom_index not in atom_indices:
                submol.RemoveAtom(atom_index)

        return Chem.MolToSmarts(submol)

    def convert_smarts_to_atoms(self, dihedral_smarts):
        """
        Given a SMARTS string that defines atoms to torsion scan, find the 4 atom indices that represent
        the same atoms for this CO
        :param dihedral_smarts: a SMARTS string
        :returns: a comma separated string of 4 atom indices (ex: '3,5,8,9')
        """
        mol = self.mol
        substruct_mol = Chem.MolFromSmarts(dihedral_smarts)
        substruct_match = mol.GetSubstructMatch(substruct_mol)

        if len(substruct_match) == 4:
            return ",".join(map(str, substruct_match))

    def pick_most_relevant_dihedral(self, pioneer_smarts, pioneer_atom_indices):
        """
        Given the smarts and atom indices of the chosen rotatable bond in a pioneer molecule, this function
        chooses the most relevant atoms in this compound occurrence. The most relevant atoms are
        those with the greatest overlap in atom indices with at most a 2-atom change from the pioneer
        smarts.
        :param pioneer_smarts: the smarts string that defines the chosen atoms in the pioneer
        :param pioneer_atom_indices: a comma-separated string of atom indices of interest in the pioneer
        :returns: a tuple (smarts, atom indices) where
            - [0] smarts: the smarts string that defines the chosen atoms (without wildcard atoms)
            - [1] atom indices: a comma-separated string of the atom indices of interest
        """
        mol = self.mol
        pioneer_atom_indices = [int(index) for index in pioneer_atom_indices.split(",")]
        all_matches = {}

        # Add any substructures that perfectly match the pioneer smarts
        substruct_mol = Chem.MolFromSmarts(pioneer_smarts)
        substruct_matches = mol.GetSubstructMatches(substruct_mol)
        all_matches[pioneer_smarts] = list(substruct_matches)

        # If a perfect match is found (smarts & atom indices), return it immediately
        for substruct_match in substruct_matches:
            if set(substruct_match) == set(pioneer_atom_indices):
                return pioneer_smarts, ",".join(str(index) for index in substruct_match)

        # Check every combination of 1 or 2 wildcard atoms and add their substructure matches
        for match1 in re.finditer("#\d+", pioneer_smarts):
            for match2 in re.finditer("#\d+", pioneer_smarts):
                if match1.start() != match2.start():
                    splits = sorted(
                        [match1.start(), match1.end(), match2.start(), match2.end()]
                    )
                    new_smarts = (
                        pioneer_smarts[: splits[0]]
                        + "*"
                        + pioneer_smarts[splits[1] : splits[2]]
                        + "*"
                        + pioneer_smarts[splits[3] :]
                    )
                else:
                    new_smarts = (
                        pioneer_smarts[: match1.start()]
                        + "*"
                        + pioneer_smarts[match1.end() :]
                    )

                if new_smarts not in all_matches:
                    substruct_mol = Chem.MolFromSmarts(new_smarts)
                    substruct_matches = mol.GetSubstructMatches(substruct_mol)
                    all_matches[new_smarts] = list(substruct_matches)

        best_overlap = 0
        best_atom_indices = None
        for smarts, matches in all_matches.items():
            for atom_indices in matches:
                # Calculate the overlap (number of atoms in common)
                atom_index_overlap = set(atom_indices).intersection(
                    set(pioneer_atom_indices)
                )
                pick_this_dihedral = False

                if len(atom_index_overlap) > best_overlap:
                    pick_this_dihedral = True
                elif (
                    len(atom_index_overlap) >= best_overlap and smarts == pioneer_smarts
                ):
                    pick_this_dihedral = True

                if pick_this_dihedral:
                    best_overlap = len(atom_index_overlap)
                    best_atom_indices = ",".join(str(index) for index in atom_indices)

        if best_atom_indices:
            return self.convert_atoms_to_smarts(best_atom_indices), best_atom_indices
        return None, None

    def generate_torsions(self, dihedral_atoms, dq_task_name):
        """
        Generate input and output file paths and submit those to the run_mc_torsion_scan analysis
        task. Returns a parsed dictionary of the analysis output_file of the format:
        {
            "delta_energy": closest_dihedral_to_input_e,
            "initial_dihedral": initial_dihedral,
            "torsions": {
                "co-pk-dihedral": {
                    "moltext": moltext,
                    "rel_energy": rel_e,
                    "torsion": torsion_dihedral }
                }
        }
        :param dihedral_atoms: a comma separated string of 4 atom indices (ex: '3,5,8,9')
        :param dq_task_name: django_q task name
        :return: dictionary of torsion scan results
        """
        # Check if this exact task has already been run successfully
        cached_task = Task.objects.filter(success=True, name=dq_task_name).first()
        if cached_task and cached_task.result and type(cached_task.result) == dict:
            return cached_task.result
        if not dihedral_atoms:
            # No torsion dihedral, return empty dict
            return {"torsions": {}, "delta_energy": ""}
        if self.molblock:
            mol = Chem.MolFromMolBlock(self.molblock, removeHs=False)
        else:
            lsalign_path = self.compound.get_lsalign_path()
            mol = Chem.SDMolSupplier(lsalign_path, removeHs=False)[0]

        # Add hydrogens if they aren't already there
        mol = Chem.rdmolops.AddHs(mol, addCoords=True)
        dihedral_smarts = self.convert_atoms_to_smarts(dihedral_atoms)
        filepath_atoms = dihedral_atoms.replace(",", "-")
        initial_dihedral = self._calculate_dihedral(mol, dihedral_atoms)

        output_filepath = compound_files_path(
            self.compound, f"torsions_{filepath_atoms}.sdf", co=self, local=True
        )
        input_filepath = compound_files_path(
            self.compound, f"torsions_input_{filepath_atoms}.sdf", co=self, local=True
        )

        writer = Chem.SDWriter(input_filepath)
        writer.write(mol)
        writer.close()

        output_filepath = run_mc_torsion_scan(
            input_filepath,
            dihedral_smarts,
            dihedral_atoms,
            output_filepath,
            dq_task_name,
        )

        if output_filepath and os.path.exists(output_filepath):
            torsion_output_dict = self.process_torsion_output(
                output_filepath, initial_dihedral
            )

            # Sort the output dihedrals in order from -180 to 180
            def get_sort_val(torsion_output):
                return int(torsion_output[0].split("-")[-1].replace("neg", "-"))

            torsion_output_dict["torsions"] = dict(
                sorted(torsion_output_dict["torsions"].items(), key=get_sort_val)
            )
            return torsion_output_dict

        return {}

    def process_torsion_output(self, output_file, initial_dihedral):
        """
        Convert the results of the torsion scan from the output sdf file to
        a dictionary format that will be used on the front end
        :param output_file: path to sdf file with the torsion scan output
        :param initial_dihedral: the dihedral of the input conformer
        :return: dictionary of results of the form
        {
            "delta_energy": closest_dihedral_to_input_e,
            "initial_dihedral": initial_dihedral,
            "torsions": {
                "co-pk-dihedral": {
                    "moltext": moltext,
                    "rel_energy": rel_e,
                    "torsion": torsion_dihedral }
                }
        }
        """
        results = {"torsions": {}, "initial_dihedral": initial_dihedral}
        suppl = Chem.SDMolSupplier(output_file, removeHs=False)
        torsion_dihedrals = []
        for mol in suppl:
            torsion = mol.GetProp("Torsion_Angle")
            conf_id = f"c{self.compound.pk}-co{self.pk}-{torsion}"
            if int(torsion) < 0:
                conf_id = f"c{self.compound.pk}-co{self.pk}-neg{abs(int(torsion))}"
            rel_e = mol.GetProp("Psi4_Relative_Energy (kcal/mol)")
            moltext = Chem.MolToMolBlock(mol)

            results["torsions"][conf_id] = {
                "moltext": moltext,
                "rel_energy": rel_e,
                "dihedral": torsion,
            }
            torsion_dihedrals.append(torsion)

        # Find energy of closest dihedral to the initial
        closest_dihedral = self._find_closest(torsion_dihedrals, initial_dihedral)
        conf_id = f"c{self.compound.pk}-co{self.pk}-{closest_dihedral}"
        if int(closest_dihedral) < 0:
            conf_id = f"c{self.compound.pk}-co{self.pk}-neg{abs(int(closest_dihedral))}"
        results["delta_energy"] = results["torsions"][conf_id]["rel_energy"]

        return results

    @staticmethod
    def _calculate_dihedral(mol, atom_indices):
        """
        Calculate the initial dihedral of the input mol
        :param mol: rdkit mol object
        :param atom_indices: atom indices of the form "1,2,3,4"
        :return: dihedral degree as a string
        """
        a, b, c, d = list(map(int, atom_indices.split(",")))
        conf = mol.GetConformer(0)
        degree = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
        return str(degree)

    @staticmethod
    def _find_closest(torsions, initial):
        """
        Find the closest dihedral in the list of torsions to the initial dihedral
        :param torsions: list of torsion dihedrals as strings
        :param initial: str(dihedral of input conformations)
        :return: closest dihedral to initial
        """
        closest = torsions[
            min(
                range(len(torsions)),
                key=lambda i: abs(int(torsions[i]) - float(initial)),
            )
        ]
        return str(closest)

    @staticmethod
    def clean_dihedrals(dihedral_smarts):
        """
        Check dihedral SMARTs pattern to make sure that atoms are not disconnected
        :param dihedral_smarts: smarts string
        :return: True if valid dihedral smarts
        """
        mol = Chem.MolFromSmarts(dihedral_smarts)

        # Check if atoms are disconnected
        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) > 1:
            return False

        return True

    #########################
    #### MMP METHODS ####
    #########################

    def get_mmp_task_name(self):
        """
        :return: the name of the MMP async task for this CompoundOccurrence
        """
        return f"{MMP}_{self.pk}"
