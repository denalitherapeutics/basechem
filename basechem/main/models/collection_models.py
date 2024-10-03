import copy
import logging
import os
import re
import sys
from io import StringIO

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.core.mail import mail_admins
from django.db import models
from django.db.models import Case, When
from django.urls import reverse
from django_q.tasks import async_task, fetch_group
from rdkit import Chem
from rdkit.rdBase import LogToPythonStderr

from basechem.common.analysis_utils import (
    collect_align_results,
    collect_dock_results,
    collect_esp_results,
    collect_torsion_results,
)
from basechem.common.constants import ADMIN_FAILURE
from basechem.common.dtx_utils import (
    check_dtx_for_inchi,
    get_agg_ic50_data,
    get_registered_structures,
)
from basechem.common.file_utils import get_tmp_filepath
from basechem.common.propcalc_utils import collect_propcalc_results
from basechem.common.rdkit_utils import RDKitWrappers
from basechem.main.constants import (
    ALIGN,
    ALOGD,
    AUTO,
    DOCK,
    DTX_MMP,
    ESP,
    HLM,
    MMP,
    PROPCALC,
    RLM,
    TORSION,
)
from basechem.main.models.compound_models import Compound, CompoundOccurrence, Series
from basechem.main.models.project_models import Project
from basechem.mni_common.storage import select_media_storage
from basechem.users.models import BasechemUser

logger = logging.getLogger("django")


def collection_files_path(instance, filename, local=False):
    """
    Directory for placing collection level files
    """
    dir_path = f"collections/{instance.owner.get_email_name()}/{instance.id}"
    if local:
        dir_path = f"{settings.MEDIA_ROOT}/{dir_path}"
        os.makedirs(dir_path, exist_ok=True)

    return f"{dir_path}/{filename}"


class Collection(models.Model):
    compound_occurrences = models.ManyToManyField(CompoundOccurrence)
    owner = models.ForeignKey(BasechemUser, on_delete=models.PROTECT)
    project = models.ForeignKey(Project, on_delete=models.PROTECT)
    created_on = models.DateTimeField(auto_now_add=True, verbose_name="Created Time")
    name = models.CharField(max_length=50, blank=True)
    sdf_file = models.FileField(
        upload_to=collection_files_path, storage=select_media_storage, blank=True
    )

    metadata = models.JSONField(null=True, blank=True, default=dict)

    def __str__(self):
        if self.name:
            return f"{self.pk}: {self.name}"
        else:
            return str(self.pk)

    def compounds(self):
        """
        Queryset of compounds for this collection
        """
        return (
            Compound.objects.filter(compoundoccurrence__collection__id=self.id)
            .distinct()
            .order_by("pk")
        )

    def get_co_order(self, cos=None):
        """
        Retrieves the value of self.metadata["co_order"], which is the desired viewing order
        for CompoundOccurrences in this collection. If `cos` contains CompoundOccurrences not
        in the current list, they are added to the list.
        :param cos: a queryset of CompoundOccurrence objects. If none, all COs in the collection are used
        :returns: a list of PKs, the desired order of CompoundOccurrence objects
        """
        if not cos:
            cos = self.compound_occurrences.all()
        co_order = self.metadata.get("co_order", [])

        for co in cos.select_related("parent_co").order_by("parent_co", "pk"):
            if co.pk not in co_order:
                if co.parent_co and co.parent_co.pk in co_order:
                    # Insert child COs directly after their parent CO when possible
                    insert_index = co_order.index(co.parent_co.pk) + 1
                    co_order = (
                        co_order[:insert_index] + [co.pk] + co_order[insert_index:]
                    )
                else:
                    co_order.append(co.pk)

        self.metadata["co_order"] = co_order
        self.save()
        return self.metadata["co_order"]

    def get_cos_in_order(self, cos=None):
        """
        Returns a queryset of CompoundOccurrence objects in the desired order (as saved in self.metadata["co_order"]).
        :param cos: a queryset of CompoundOccurrence objects. If none, all COs in the collection are used
        :returns: a CompoundOccurrence queryset, in the desired order
        """
        if not cos:
            cos = self.compound_occurrences.all()
        co_list = self.get_co_order(cos=cos)
        if co_list:
            preserved = Case(*[When(pk=pk, then=pos) for pos, pk in enumerate(co_list)])
            return cos.order_by(preserved)
        return cos

    def update_co_order(self, new_co_order):
        """
        Updates self.metadata["co_order"] so that the new_co_order is reflected in the object
        metadata. `new_co_order` may be a subset of CompoundOccurrences in this collection,
        in which case other CompoundOccurrence PKs may appear at the end of the list.
        :param new_co_order: a list of CompoundOccurrence PKs
        :returns: the new self.metadata["co_order"]
        """
        co_order = self.metadata.get("co_order", [])
        co_order = [pk for pk in co_order if pk not in new_co_order]
        self.metadata["co_order"] = new_co_order + co_order
        self.save()
        return self.metadata["co_order"]

    def get_cos_for_analysis(self, analysis):
        """
        Given an analysis, return CompoundOccurrences that should be shown in the table
        :param analysis: a string, the analysis being performed (align, dock, torsion, etc.)
        :returns: a queryset of CompoundOccurrences
        """
        if analysis in [TORSION, ESP]:
            return_pks = []
            # Return all basechem generated COs from other analyses, user uploaded 3D COs,
            # and user uploaded 2D COs with no child COs
            for co in self.compound_occurrences.exclude(saved_from=analysis):
                if (
                    not co.parent_co  # is a parent
                    and not co.molblock  # is 2D
                    # Has children that will be part of this analysis
                    and co.compoundoccurrence_set.filter(collection=self)
                    .exclude(saved_from=analysis)
                    .exists()
                ):
                    # If 2D and there are children in the collection, don't include
                    continue
                return_pks.append(co.pk)
            cos = CompoundOccurrence.objects.filter(pk__in=return_pks)
            return self.get_cos_in_order(cos=cos)
        elif analysis in [PROPCALC, ALIGN, DOCK]:
            co_pks = (
                self.compound_occurrences.filter(parent_co=None)
                .order_by("compound")
                .distinct("compound")
                .values_list("pk", flat=True)
            )
            cos = self.compound_occurrences.filter(pk__in=co_pks)
            return self.get_cos_in_order(cos=cos)
        elif analysis == MMP:
            # Currently the MMP analysis only supports one CompoundOccurrence at a time. Choose
            # a CompoundOccurrence that has data in the mmp_analysis metadata OR an arbitrary CompoundOccurrence
            # in the Collection if no CompoundOccurrence has data in the mmp_analysis metadata
            if len(self.metadata.get("mmp_analysis", {}).keys()):
                comp_pk = int(list(self.metadata["mmp_analysis"].keys())[0])
                cos = self.compound_occurrences.filter(
                    parent_co=None, compound__pk=comp_pk
                )
            else:
                cos = self.compound_occurrences.filter(parent_co=None)
            return self.compound_occurrences.filter(pk__in=cos.values_list("pk")[:1])

    def get_co_pks_used_in_group(self, group_name):
        """
        Given a group name, return a list of PKs of CompoundOccurrences that were included
        in the task group. This is used when people revisit analysis results to make sure that
        any post-analysis CompoundOccurrences don't show up in the table.
        :param group_name: the name of a djangoQ group of tasks
        :returns: a list of PKs
        """
        co_pks_used = []
        analysis = group_name.split("_")[0]
        for task in fetch_group(group_name, failures=True):
            if analysis in [ESP, TORSION]:
                co_pks_used.append(int(task.name.split("_")[1]))
        return co_pks_used

    @staticmethod
    def handle_sdf_upload(file):
        """
        Helper to handle sdf file uploads. LogToPythonStderr() along with the sio/stderr
        tracking is to capture RDKit warnings and errors based on this blog post:
        http://rdkit.blogspot.com/2016/03/capturing-error-information.html
        :param file: uploaded file
        :return: List of tuples of (RDKit Mol objects, 2D?)
        """
        tmpname = get_tmp_filepath(file, "sdf_upload.sdf")
        romols = []
        suppl = Chem.SDMolSupplier(tmpname, removeHs=False)

        LogToPythonStderr()
        sio = sys.stderr = StringIO()

        for mol in suppl:
            if mol:
                cleaned = RDKitWrappers.clean_mol_object(mol)
                romols.extend(cleaned)
                sio = sys.stderr = StringIO()

        sys.stderr = sys.__stderr__

        return romols

    def save_romols_to_sdf(self, romols):
        """
        Writes mol objects to an sdf file and saves it as the collection's `sdf_file`
        :param romols: list of tuples of (rdkit Mol objects, 2d?)
        """
        # Write mol objects to an SDF file to save original properties
        filename = f"{self.pk}_original_mols.sdf"
        localpath = collection_files_path(self, filename, local=True)
        writer = Chem.SDWriter(localpath)
        for mol, _ in romols:
            writer.write(mol)
        writer.close()
        with open(localpath, "rb") as data:
            self.sdf_file.save(filename, data)

    def new_collection_async(self, romols):
        """
        This wrapper is run asynchronously when a new collection is created. This prevents timeouts
        with large collections, because saving romols to the database can take several minutes due
        to Dotmatics calls and other parsing.
        :param romols: list of tuples of (rdkit Mol objects, bool(2d?))
        :returns: a dictionary of the form {"error": error}
        """
        error = None
        try:
            self.handle_romols(romols)
            # If this Collection was made from the MMPSubmitForm, add the PK of the Compound
            # added in handle_romols to the mmp_analysis metadata
            if "mmp_analysis" in self.metadata.keys():
                self.metadata["mmp_analysis"] = {
                    self.compounds().first().pk: self.metadata["mmp_analysis"]
                }
                self.save()
        except Exception as e:
            error = f"{type(e).__name__}: {e}"
        return {"error": error}

    def handle_romols(self, romols, test=False):
        """
        Helper to add CompoundOccurrences to the Collection from romols
        :param romols: list of tuples of (rdkit Mol objects, 2d?)
        :param test: True if this function is called for testing
        :returns: a list of CompoundOccurrences that were added to the collection
        """
        cos_added = []
        existing_cos = [co.pk for co in self.compound_occurrences.all()]
        for mol, twoD in romols:
            comp = self._create_compound(mol, test)
            if not comp:
                # There was an issue finding the related Compound object. Admins have been
                # emailed, skip this mol
                continue
            # Set `is_mmpdb_analogue` to False in case it was previously created by an MMP analysis
            # since now it has been uploaded by a user
            comp.is_mmpdb_analogue = False
            comp.save()
            cos = CompoundOccurrence.parents.filter(owner=self.owner, compound=comp)
            # If no CO exists, make one either 2D or 3D
            if cos.count() == 0:
                co = CompoundOccurrence.objects.create(compound=comp, owner=self.owner)
                if not twoD:
                    co.molblock = Chem.MolToMolBlock(mol)
                    co.save()

            elif twoD:
                # Find if any COs are 2D and reuse that if it exists
                cos = cos.filter(molblock="")
                if cos.exists():
                    co = cos.last()
                else:
                    # Otherwise make a new CO for the 2D comp
                    co = CompoundOccurrence.objects.create(
                        compound=comp, owner=self.owner
                    )
            # If 3D look for identical molblock, otherwise create a new CO
            else:
                molblock = Chem.MolToMolBlock(mol)
                cos = cos.filter(molblock=molblock)
                if cos.exists():
                    co = cos.last()
                else:
                    co = CompoundOccurrence.objects.create(
                        compound=comp, owner=self.owner, molblock=molblock
                    )
            self.compound_occurrences.add(co)
            if co.pk not in existing_cos and co not in cos_added:
                cos_added.append(co)
        return cos_added

    def _create_compound(self, mol, test=False):
        """
        Helper to get or create Compound objects for the given mol, only called in handle_romols()
        :param mol: rdkit Mol object
        :return: Compound object
        """
        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol, options="/suu")
        try:
            dn_id = mol.GetProp("dn_id")
        except:
            dn_id = ""

        matches = Compound.objects.filter(inchi=inchi)

        if len(matches) == 1:
            compound = matches[0]
            # Update DN when applicable from either prop value or dtx search
            if dn_id and not compound.dn_id:
                compound.update_dn_id(dn_id)
            if not compound.dn_id and not test:
                compound.update_dn_id(check_dtx_for_inchi(compound.inchi))
            # Update Project and Series if this compound was previously only Auto uploaded
            if compound.project.code == AUTO and not test:
                compound.project = self.project
                compound.series = compound.pick_series()
            # Try to set a series if the compound doesn't have one (probably b/c no series for Project when uploaded)
            if not compound.series and not test:
                compound.series = compound.pick_series()
            compound.save()

        elif len(matches) == 0:
            if not dn_id and not test:
                dn_id = check_dtx_for_inchi(inchi)
            compound = Compound.objects.create(
                inchi=inchi, smiles=smiles, project=self.project, dn_id=dn_id
            )
            # Generate the SDF file with the Dotmatics or user moltext so that it has preferable orientation
            if dn_id:
                dn_mol = Chem.MolFromMolBlock(get_registered_structures([dn_id])[dn_id])
                compound.get_sdf_file(dn_mol)
            else:
                compound.get_sdf_file(mol)
            compound.generate_svg_file()
        else:
            message = f"While creating Compounds for Collection {self.pk}, encountered an error: there are multiple Compounds ({', '.join([m.name for m in matches])}) with the inchi {inchi}"
            logger.error(message)
            mail_admins(ADMIN_FAILURE, message)
            return None

        return compound

    def get_sdf_file(self):
        """
        Gets this collection's sdf file. If this file does not exist, generates a file using
        the mol objects of the relevant compounds
        """
        filename = f"{self.pk}_original_mols.sdf"
        localpath = collection_files_path(self, filename, local=True)

        if not self.sdf_file:
            if not os.path.exists(localpath):
                writer = Chem.SDWriter(localpath)
                for c in self.compounds():
                    writer.write(c.mol())
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

    def most_relevant_series(self, reference, valid_series):
        """
        :param reference: the reference (ex: "s-10") used to generate results
        :param valid_series: a list of valid series identifiers to choose (ex: ["s-1", "s-10"])
        :returns: the most relevant series for this collection (ex: "s-10")
        """
        initial_series = None
        if reference and reference != "default":
            initial_series = reference
        else:
            # Count the number of times each series appears in this collection
            series_counts = list(
                self.compounds()
                .values("series__pk")
                .annotate(count=models.Count("series__pk"))
                .order_by("-count")
            )
            for series_count in series_counts:
                if f"s-{series_count['series__pk']}" in valid_series:
                    # Pick the series that appears the most times and is also a valid series
                    initial_series = f"s-{series_count['series__pk']}"
                    break

        # If no valid series is picked, use the first valid series
        if initial_series not in valid_series and len(valid_series) > 0:
            initial_series = valid_series[0]

        return initial_series

    def generate_file(
        self, current_view, test=False, group_name=None, selected_ids=None
    ):
        """
        Generates an sdf file based on the current view the user is on
        :param current_view: the current view the user is on (as url name)
        :param test: true if this method is called in a test, default False
        :param group_name: the name of the task group that is associated w/ the data to download
        :param selected_ids: a string with identifiers referring to which entities should be
            included in the download (used for downloading selected in ligand align)
        :return: a path to the generated file (including filename), the name of the file
        """
        filename = "tmp.sdf"
        filepath = collection_files_path(self, filename, local=True)
        if selected_ids:
            # convert the underscore-separated string of IDs ("co-142_co-143_s-11") to a list of IDs (["co-142", "co-143", "s-11"])
            selected_ids = selected_ids.split("_")

        if current_view == PROPCALC:
            filepath, filename = self._generate_allprops_sdf_file(group_name)

        elif current_view == ALIGN:
            filepath, filename = self._generate_align_sdf_file(group_name, selected_ids)

        elif current_view == DOCK:
            filepath, filename = self._generate_dock_sdf_file(group_name, selected_ids)

        elif current_view == ESP:
            filepath, filename = self._generate_esp_pqr_file(group_name, selected_ids)

        elif current_view == TORSION:
            filepath, filename = self._generate_torsion_sdf_file(
                group_name, selected_ids
            )

        elif current_view == DTX_MMP:
            filepath, filename = self._generate_assay_mmp_file()

        return filepath, filename

    def get_url(self, analysis_type, group_name=None):
        """
        Given an analysis type, returns the full url where the results of the analysis can be found
        :param analysis_type: the name of a url conf that includes a `collection_id`
        :param group_name: optional, the group name of the task that completed before sending this email
        """
        url = reverse(analysis_type, kwargs={"collection_id": self.id})
        if analysis_type != PROPCALC and group_name:
            url = reverse(
                analysis_type,
                kwargs={"collection_id": self.id, "group_name": group_name},
            )
        return settings.BASE_URL + url

    def run_analysis(self, current_view, **kwargs):
        """
        Wrapper to run a collection analysis asynchronously via DjangoQ task groups
        :param current_view: the current view the user is on which indicates analysis to perform
        :param kwargs: other kwargs needed for an analysis
        """
        if current_view == PROPCALC:
            self.propcalc_analysis()
        elif current_view == ALIGN:
            self.align_analysis(kwargs["ref_string"])
        elif current_view == DOCK:
            self.dock_analysis(kwargs["ref_string"])
        elif current_view == ESP:
            self.esp_analysis(kwargs["ref_string"])
        elif current_view == TORSION:
            self.torsion_analysis(
                kwargs["pioneer_pk"], kwargs["pioneer_dihedral_atoms"]
            )
        elif current_view == MMP:
            self.mmp_analysis()

    @staticmethod
    def _get_ref_obj(ref_string):
        """
        Returns the Series or Compound object from the given reference string
        """
        if ref_string == "default":
            return None

        ref, id = ref_string.split("-")
        if "s" in ref:
            return Series.objects.get(id=id)
        elif "c" in ref:
            return Compound.objects.get(id=id)
        else:
            return None

    ######################
    ##    MMP METHODS   ##
    ######################

    def get_mmp_group_name(self, comp_pk):
        """
        :param comp_pk: an integer, the PK of the Compound object whose mmps are being found in this task
        :return: a string, the name of the async group task for MMPs for this collection
        """
        return f"{MMP}_{self.pk}_{comp_pk}"

    def mmp_analysis(self):
        """
        Launch an MMP analysis for this Collection by starting an asynchronous tasks for the
        Compounds in this Collection
        """
        cos = self.get_cos_for_analysis(MMP)
        for co in cos:
            mmp_analysis_dict = self.metadata.get("mmp_analysis", {}).get(
                str(co.compound.pk), {}
            )
            async_task(
                co.compound.mmp_analysis,
                constant_smiles=mmp_analysis_dict.get("constant_smiles"),
                variable_smiles=mmp_analysis_dict.get("variable_smiles"),
                task_name=co.get_mmp_task_name(),
                group=self.get_mmp_group_name(co.compound.pk),
                hook="basechem.main.tasks.task_completion_hook",
            )

    def find_mmps(self):
        """
        Calls `find_mmps` for each Compound in this Collection
        """
        for comp in self.compounds():
            comp.find_mmps()

    def update_mmp_dtx_avg_assay_data(self, dns_to_skip=None):
        """
        Update average assay values for DTX Compounds in this Collection and all DTX MMPs
        of those Compounds, skipping any DNs that appear in `dns_to_skip`
        :param dns_to_skip: a list of DN IDs (strings) to skip. This is used for the newly assayed Compounds
            because brand-new assay data is not yet in the aggregate datasource, so we show the most
            recent assay value instead of aggregate data.
        """
        comps = self.compounds().prefetch_related("mmps")
        dn_ids = list(comps.exclude(dn_id="").values_list("dn_id", flat=True))
        for c in comps:
            dn_ids.extend(
                list(c.mmps.exclude(dn_id="").values_list("dn_id", flat=True))
            )
        ic50_data = get_agg_ic50_data(dn_ids)
        for c in comps:
            c.update_mmp_dtx_avg_assay_data(ic50_data, dns_to_skip)

    ######################
    ## PROPCALC METHODS ##
    ######################

    def get_propcalc_group_name(self):
        """
        :return: the name of the async group task for propcalc for this collection
        """
        return f"{PROPCALC}_{self.id}"

    def get_propcalc_column_headers(self):
        """
        The column headers in PROPCALC differ slightly from the name of the props (e.g. the RLM
        prop has two columns, "RLM Prediction" and "RLM Probabilities"). This function returns a list
        of column headers based on the props for this collection
        :returns: a list of column headers for propcalc
        """
        props = copy.deepcopy(self.metadata.get("props_to_show", []))
        # Some properties require multiple columns to display the data, adjust columns accordingly
        if ALOGD in props:
            props.append(f"{ALOGD} Prediction")
            props.remove(ALOGD)
        if RLM in props:
            props.extend([f"{RLM} Prediction", f"{RLM} Probabilities"])
            props.remove(RLM)
        if HLM in props:
            props.extend([f"{HLM} Prediction", f"{HLM} Probabilities"])
            props.remove(HLM)
        return props

    def inductive_in_props(self):
        """
        Are any of the props selected calculated using the InductiveBio models
        :returns: a boolean, do this collection's properties require InductiveBio
        """
        props = self.metadata.get("props_to_show", [])
        inductive = any([HLM in props, RLM in props, ALOGD in props])
        return inductive and settings.INDUCTIVE_BIO_ENABLED

    def propcalc_analysis(self):
        """
        Launch a property calculation analysis for this collection by starting an asynchronous
        task for each compound in the collection
        """
        props = self.metadata.get("props_to_show", [])
        if props:
            for co in self.get_cos_for_analysis(PROPCALC):
                async_task(
                    co.compound.run_propcalc,
                    inductive=self.inductive_in_props(),
                    task_name=co.get_propcalc_task_name(self),
                    group=self.get_propcalc_group_name(),
                    hook="basechem.main.tasks.task_completion_hook",
                    cluster="fast",
                )

    def _generate_allprops_sdf_file(self, group_name=None):
        """
        Generates an sdf file with all properties for all compounds in this collection
        :param group_name: optional, the name of the task group that calculated properties (so they don't need to be recalculated)
        :return: path to file and filename
        """
        filename = f"{self.id}_propcalc_results.sdf"
        filepath = collection_files_path(self, filename, local=True)

        writer = Chem.SDWriter(filepath)
        for c in self.compounds():
            if group_name:
                # Use the pre-computed properties from the task results
                propcalc_results = collect_propcalc_results(
                    fetch_group(group_name, failures=True)
                )
                mol = c.mol(hydrogens=True)
                for prop, value in propcalc_results.get(c.id, {}).items():
                    mol.SetProp(prop, str(value))

            else:
                inductive = self.inductive_in_props()
                mol = c.mol_w_properties(inductive=inductive, hydrogens=True)
            # Remove image paths from sdf download
            mol.ClearProp("rlm_probabilities")
            mol.ClearProp("hlm_probabilities")
            mol.ClearProp("rlm_interpretation")
            mol.ClearProp("hlm_interpretation")

            writer.write(mol)
        writer.close()

        return filepath, filename

    def get_refs_task_name(self, analysis):
        """
        :param analysis: a string, the analysis type (ex: ALIGN, DOCK, ESP)
        :returns: a string, the task_name for a "references" DjangoQ task for this collection (ex: "alignrefs_10")
        """
        return f"{analysis}refs_{self.pk}"

    #####################
    ### ALIGN METHODS ###
    #####################

    def get_align_group_name(self, ref_string):
        """
        :param ref_string: a string:
            - 'default': align all compounds to their assigned series
            - 's-(int)': align all compounds to the series w/ id (int)
        :return: the name of the async group for alignment with the ref_string
        """
        return f"{ALIGN}_{self.id}_{ref_string}"

    def align_analysis(self, ref_string):
        """
        Launch an async task for superimposing each conf in the collection and
        a task for formatting the reference compound
        :param ref_string: a string:
            - 'default': align all compounds to their assigned series
            - 's-(int)': align all compounds to the series w/ id (int)
        """
        ref_obj = self._get_ref_obj(ref_string)

        async_task(
            self._collect_align_references,
            task_name=self.get_refs_task_name(ALIGN),
            group=self.get_align_group_name(ref_string),
            hook="basechem.main.tasks.task_completion_hook",
        )

        for co in self.get_cos_for_analysis(ALIGN):
            async_task(
                co.superimpose_to_ref,
                reference=ref_obj,
                task_name=co.get_align_task_name(ref_string),
                group=self.get_align_group_name(ref_string),
                hook="basechem.main.tasks.task_completion_hook",
            )

    def _collect_align_references(self):
        """
        Helper to use when launching an align job to collect reference compounds
        """
        results = {"references": {}, "receptors": {}}
        series = Series.objects.filter(project=self.project, active=True)

        for s in series:
            results["references"][f"s-{s.pk}"] = Chem.MolToMolBlock(s.mol())
            if s.receptor_file:
                receptor = s.pdb_block()
                residues = []
                if s.receptor_residues:
                    residues = list(map(int, s.receptor_residues.split(",")))
                results["receptors"][f"s-{s.pk}"] = {
                    "moltext": receptor,
                    "residues": residues,
                }

        return results

    def _generate_align_sdf_file(self, group_name, selected_ids=None):
        """
        Generates an sdf file with all the aligned conformers of all the compounds in
        this collection
        :param group_name: a string, the name of the task group that generated this alignment
        :param selected_ids: a list of ids ["id-1", "id-2", "id-3"], where each id refers
            to an entity (reference or conformer) that should be included in the download
        :return: a tuple - (path to file, filename)
        """
        reference = group_name.split("_")[2]
        if selected_ids:
            filename = f"{self.id}_align_{reference}_selected.sdf"
        else:
            filename = f"{self.id}_align_{reference}_all.sdf"
        filepath = collection_files_path(self, filename, local=True)

        if os.path.exists(filepath):
            return filepath, filename

        align_results = collect_align_results(fetch_group(group_name, failures=True))

        writer = Chem.SDWriter(filepath)
        if align_results:
            for c_id, conformers in align_results["compounds"].items():
                for conf_id, data in conformers.items():
                    if not selected_ids or (selected_ids and conf_id in selected_ids):
                        mol = Chem.MolFromMolBlock(data["moltext"], removeHs=False)
                        mol = Chem.AddHs(mol, addCoords=True)
                        mol.SetProp("conf_id", conf_id)
                        mol.SetProp("r_mmff_rel_energy", data["r_mmff_rel_energy"])
                        mol.SetProp(
                            "r_bc_rmsd_to_lsalign", data["r_bc_rmsd_to_lsalign"]
                        )
                        writer.write(mol)

            for ref_id, moltext in align_results["references"].items():
                if not selected_ids or (selected_ids and ref_id in selected_ids):
                    mol = Chem.MolFromMolBlock(moltext, removeHs=False)
                    mol = Chem.AddHs(mol, addCoords=True)
                    writer.write(mol)
        writer.close()

        return filepath, filename

    #####################
    ### DOCK METHODS ####
    #####################

    def get_dock_group_name(self, ref_string):
        """
        :param ref_string: a string:
            - 'default': dock all compounds to their assigned series
            - 's-(int)': dock all compounds to the series w/ id (int)
        :return: the name of the async group for docking with the ref_string
        """
        return f"{DOCK}_{self.id}_{ref_string}"

    def dock_analysis(self, ref_string):
        """
        Launch an async task for docking each compound in the collection and
        a task for docking the reference compound
        :param ref_string: a string:
            - 'default': dock all compounds to their assigned series
            - 's-(int)': dock all compounds to the series w/ id (int)
        """
        ref_obj = self._get_ref_obj(ref_string)

        # Make sure all series available for this collection have the needed rDock setup files
        series = Series.objects.filter(project=self.project, active=True).exclude(
            receptor_file_mol2=""
        )
        for s in series:
            s.get_rdock_prm()

        async_task(
            self._collect_dock_references,
            task_name=self.get_refs_task_name(DOCK),
            group=self.get_dock_group_name(ref_string),
            hook="basechem.main.tasks.task_completion_hook",
        )

        for co in self.get_cos_for_analysis(DOCK):
            async_task(
                co.dock_to_receptor,
                reference=ref_obj,
                task_name=co.get_dock_task_name(ref_string),
                group=self.get_dock_group_name(ref_string),
                hook="basechem.main.tasks.task_completion_hook",
            )

    def _collect_dock_references(self):
        """
        Helper to use when launching a dock job to collect reference compounds
        :returns: dict with moltext for all Series' for this Collections' Project and
        receptors for the Series' that were assigned to the Collections' Compounds
        """
        results = {"references": {}, "receptors": {}}
        series = Series.objects.filter(project=self.project, active=True).exclude(
            receptor_file_mol2=""
        )

        for s in series:
            results["references"][f"s-{s.pk}"] = Chem.MolToMolBlock(s.mol())
            if s.receptor_file_mol2:
                receptor = s.pdb_block()
                residues = []
                if s.receptor_residues:
                    residues = list(map(int, s.receptor_residues.split(",")))
                results["receptors"][f"s-{s.pk}"] = {
                    "moltext": receptor,
                    "residues": residues,
                }

        return results

    def _generate_dock_sdf_file(self, group_name, selected_ids=None):
        """
        Generates an sdf file with the docked poses in this collection
        :param group_name: a string, the name of the task group for this analysis
        :param selected_ids: a list of ids ["id-1", "id-2", "id-3"], where each id refers to a pose
            that should be included in the download. If empty or None, includes all poses in the collection.
        :return: a tuple - (path to file, filename)
        """
        reference = group_name.split("_")[2]
        if selected_ids:
            filename = f"{self.id}_{DOCK}_{reference}_selected.sdf"
        else:
            filename = f"{self.id}_{DOCK}_{reference}_all.sdf"
        filepath = collection_files_path(self, filename, local=True)

        if os.path.exists(filepath):
            return filepath, filename

        dock_results = collect_dock_results(fetch_group(group_name, failures=True))

        writer = Chem.SDWriter(filepath)
        if dock_results:
            for c_id, poses in dock_results["compounds"].items():
                for pose_id, data in poses.items():
                    if not selected_ids or (selected_ids and pose_id in selected_ids):
                        mol = Chem.MolFromMolBlock(data["moltext"], removeHs=False)
                        mol = Chem.AddHs(mol, addCoords=True)
                        mol.SetProp("pose_id", pose_id)
                        mol.SetProp("toklatScore", data["toklatScore"])
                        mol.SetProp("rdockScore", data["rdockScore"])
                        mol.SetProp("RMSDtoLSAligned", data["RMSDtoLSAligned"])
                        writer.write(mol)

            for ref_id, moltext in dock_results["references"].items():
                if not selected_ids or (selected_ids and ref_id in selected_ids):
                    mol = Chem.MolFromMolBlock(moltext, removeHs=False)
                    mol = Chem.AddHs(mol, addCoords=True)
                    writer.write(mol)
        writer.close()

        return filepath, filename

    #####################
    ### ESP METHODS ####
    #####################

    def get_esp_group_name(self, ref_string):
        """
        :param ref_string: a string:
            - 'default': generate ESP maps using the series assigned to each compound all compounds to their assigned series
            - 's-(int)': generate ESP maps using the series with the given id (int)
        :return: the name of the async group for esp mapping with the ref_string
        """
        return f"{ESP}_{self.id}_{ref_string}"

    def esp_analysis(self, ref_string):
        """
        Launch an async task for generating ESP maps for each compound in the collection
        :param ref_string: a string:
            - 'default': generate ESP maps using the series assigned to each compound all compounds to their assigned series
            - 's-(int)': generate ESP maps using the series with the given id (int)
        """
        ref_obj = self._get_ref_obj(ref_string)

        async_task(
            self._collect_esp_references,
            task_name=self.get_refs_task_name(ESP),
            group=self.get_esp_group_name(ref_string),
            hook="basechem.main.tasks.task_completion_hook",
        )

        for co in self.get_cos_for_analysis(ESP):
            task_name = co.get_esp_task_name(ref_string)
            async_task(
                co.generate_esp_map,
                reference=ref_obj,
                job_name=task_name,
                task_name=task_name,
                group=self.get_esp_group_name(ref_string),
                hook="basechem.main.tasks.task_completion_hook",
            )

    def _collect_esp_references(self):
        """
        Helper to use when launching an ESP job to collect reference compounds
        :returns: dict with pqr text for all Series' for this Collections' Project
        """
        results = {"references": {}, "receptors": {}}
        series = Series.objects.filter(project=self.project, active=True)

        for s in series:
            esp_map_data = s.generate_ligand_esp_map()
            if not esp_map_data.get("dx"):
                logger.error(
                    f"Skipping series {series.name} in ESP map generation because APBS failed to generate a DX file."
                )
                continue
            results["references"][f"s-{s.pk}"] = esp_map_data
            if s.receptor_file:
                residues = []
                if s.receptor_residues:
                    residues = list(map(int, s.receptor_residues.split(",")))
                results["receptors"][f"s-{s.pk}"] = {
                    "pqr": s.generate_receptor_esp_map(),
                    "residues": residues,
                }

        return results

    def _generate_esp_pqr_file(self, group_name, selected_ids=None):
        """
        Generates a pqr file with the ESP maps in this collection
        :param group_name: a string, the name of the task group for this analysis
        :param selected_ids: a list of ids ["id-1", "id-2", "id-3"], where each id refers to a compound/reference
            whose ESP map should be included in the download. If empty or None, includes all ESP maps in the collection.
        :return: a tuple - (path to file, filename)
        """
        reference = group_name.split("_")[2]
        if selected_ids:
            filename = f"{self.id}_{ESP}_{reference}_selected.pqr"
        else:
            filename = f"{self.id}_{ESP}_{reference}_all.pqr"
        filepath = collection_files_path(self, filename, local=True)

        if os.path.exists(filepath):
            return filepath, filename

        esp_results = collect_esp_results(fetch_group(group_name, failures=True))

        if esp_results:
            with open(filepath, "w") as outfile:
                for c_id, c_data in esp_results["compounds"].items():
                    if not selected_ids or (selected_ids and c_id in selected_ids):
                        outfile.writelines(c_data["pqr"])

                for ref_id, pqr in esp_results["references"].items():
                    if not selected_ids or (selected_ids and ref_id in selected_ids):
                        outfile.writelines(pqr)
                outfile.close()

        return filepath, filename

    #######################
    ### TORSION METHODS ###
    #######################

    def get_torsion_group_name(self, pioneer_pk, pioneer_dihedral_atoms):
        """
        :param pioneer_pk: pk of the CompoundOccurrence object selected as the pioneer for torsion analysis
        :param pioneer_dihedral_atoms: a comma separated string of 4 atom indices (ex: '3,5,8,9')
        :return: the name of the async group for torsion scanning
        """
        url_atoms_str = pioneer_dihedral_atoms.replace(",", "-")
        return f"{TORSION}_{self.id}_{pioneer_pk}_{url_atoms_str}"

    def torsion_analysis(self, pioneer_pk, pioneer_dihedral_atoms):
        """
        Launch an async task for generating torsion scans for each CompoundOccurrence in the Collection
        :param pioneer_pk: pk of the CompoundOccurrence object selected as the pioneer for torsion analysis
        :param pioneer_dihedral_atoms: a comma separated string of 4 atom indices (ex: '3,5,8,9')
        """
        pioneer = CompoundOccurrence.objects.get(pk=pioneer_pk)
        pioneer_dihedral_smarts = pioneer.convert_atoms_to_smarts(
            pioneer_dihedral_atoms
        )

        for co in self.get_cos_for_analysis(TORSION):
            if co.pk == pioneer_pk:
                dihedral_atoms = pioneer_dihedral_atoms
            else:
                _, dihedral_atoms = co.pick_most_relevant_dihedral(
                    pioneer_dihedral_smarts, pioneer_dihedral_atoms
                )

            task_name = co.get_torsion_task_name(
                pioneer_dihedral_smarts, dihedral_atoms
            )
            async_task(
                co.generate_torsions,
                dihedral_atoms=dihedral_atoms,
                dq_task_name=task_name,
                group=self.get_torsion_group_name(pioneer_pk, pioneer_dihedral_atoms),
                hook="basechem.main.tasks.task_completion_hook",
                task_name=task_name,
                cluster="slow",
            )

    def _generate_torsion_sdf_file(self, group_name, selected_ids=None):
        """
        Generates an sdf file with the torsion scan results in this collection
        :param group_name: a string, the name of the task group for this analysis
        :param selected_ids: a list of ids ["co-pk1-dih1", "co-pk1-dih2", "co-pk2-dih2"], where each id
            refers to a particular torsion dihedral of a CompoundOccurrence that should be
            included in the download. If empty or None, includes all torsions
        :return: a tuple - (path to file, filename)
        """
        if selected_ids:
            filename = f"{self.id}_{TORSION}_selected.sdf"
        else:
            filename = f"{self.id}_{TORSION}_all.sdf"
            selected_ids = []
        filepath = collection_files_path(self, filename, local=True)

        if os.path.exists(filepath):
            return filepath, filename

        torsion_results = collect_torsion_results(
            fetch_group(group_name, failures=True)
        )

        writer = Chem.SDWriter(filepath)
        if torsion_results:
            for co_id in torsion_results.keys():
                torsion_dict = torsion_results[co_id]["torsions"]
                ids_regex = re.compile(r"(c\d+-co\d+-((neg)?\d+))")
                conf_ids_to_download = [
                    re.search(ids_regex, conf_id).groups()[0]
                    for conf_id in selected_ids
                    if co_id.replace("-", "") in conf_id
                ]
                # If there aren't any selected conformers, use all conformers
                if not selected_ids:
                    conf_ids_to_download = list(torsion_dict.keys())
                for conf_id in conf_ids_to_download:
                    mol = Chem.MolFromMolBlock(torsion_dict[conf_id]["moltext"])
                    mol.SetProp("dihedralAngle", torsion_dict[conf_id]["dihedral"])
                    mol.SetProp(
                        "Psi4_Relative_Energy (kcal/mol)",
                        torsion_dict[conf_id]["rel_energy"],
                    )
                    writer.write(mol)

        writer.close()

        return filepath, filename

    def _generate_assay_mmp_file(self):
        """
        Generates an sdf file with all MMP results for each Compound in this Collection
        :return: a tuple - (path to file, filename)
        """
        filename = f"{self.id}_{self.project.code}_{DTX_MMP}_all.sdf"
        filepath = collection_files_path(self, filename, local=True)

        if os.path.exists(filepath):
            return filepath, filename

        writer = Chem.SDWriter(filepath)
        compounds = self.compounds().order_by("dn_id")
        compound_pks = compounds.values_list("pk", flat=True)
        mmps_already_written = []

        for compound in compounds:
            c_mol = compound.mol()
            mmps = compound.mmps.exclude(dn_id="").order_by("dn_id")
            mmp_names = mmps.values_list("dn_id", flat=True)

            # Write compound and relevant info
            c_mol.SetProp("DN_MMPs", ", ".join(list(mmp_names)))
            c_mol.SetProp("Assay_Date", self.created_on.strftime("%m/%d/%Y"))

            assay_results = compound.measured_data.get("assay_results", {})
            for assay, assay_data in assay_results.items():
                for analysis, value in assay_data.items():
                    c_mol.SetProp(f"{assay}: {analysis}", str(value))
            writer.write(c_mol)

            # Write MMPs that haven't been written and relevant info
            for mmp in mmps:
                if mmp.pk not in compound_pks and mmp.pk not in mmps_already_written:
                    mmp_mol = mmp.mol()
                    assay_results = mmp.measured_data.get("assay_results", {})
                    for assay, assay_data in assay_results.items():
                        for analysis, value in assay_data.items():
                            mmp_mol.SetProp(f"{assay}: {analysis}", str(value))
                    writer.write(mmp_mol)
                    mmps_already_written.append(mmp.pk)

        writer.close()

        return filepath, filename
