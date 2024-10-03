from pathlib import Path

from django import forms
from django.contrib import admin
from django.core.exceptions import ObjectDoesNotExist
from django.core.files.base import File
from django.db import models
from django.db.models import Prefetch
from django.urls import reverse
from django.utils.safestring import mark_safe
from rdkit import Chem

from basechem.common.file_utils import get_tmp_filepath
from basechem.main.admin_filters import (
    CollectionFilter,
    CompoundDnIdFilter,
    CompoundOccurrenceFilter,
    CompoundPkFilter,
    ProjectFilter,
    UserFilter,
)
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, CompoundOccurrence, Series
from basechem.main.models.project_models import Project

#############################################
#                Admin Forms                #
#############################################


class DeleteOldFilesModelForm(forms.ModelForm):
    """
    This model form base class overrides the save method so that previous files are deleted
    if a file field is updated.
    """

    def save(self, commit=True):
        object_cls = type(self.instance)
        try:
            old_object = object_cls.objects.get(pk=self.instance.pk)
        except ObjectDoesNotExist:
            # First save, no previous object to compare too
            return super().save(commit=commit)
        for field in self.changed_data:
            if isinstance(old_object._meta.get_field(field), models.FileField):
                old_value = getattr(old_object, field)
                if old_value:
                    try:
                        old_value.delete()
                    except:
                        pass  # File doesn't exist
        return super().save(commit=commit)


class SeriesAdminForm(DeleteOldFilesModelForm):
    def clean_smiles(self):
        smiles = self.cleaned_data.get("smiles")
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise forms.ValidationError("Are you sure this is a valid SMILES string?")
        return smiles

    def clean_receptor_residues(self):
        residues_str = self.cleaned_data.get("receptor_residues")
        if residues_str:
            try:
                residues = residues_str.split(",")
                residues = [int(r.strip()) for r in residues]
            except:
                raise forms.ValidationError(
                    "Are you sure this is a comma-separated list of integers?"
                )
        return residues_str

    def clean(self, *args, **kwargs):
        """
        Help our maintainers avoid mistakes by checking that the required fields are filled appropriately.
        """
        cleaned_data = super().clean(*args, **kwargs)
        ground_state_file = cleaned_data.get("ground_state_file")
        bound_state_file = cleaned_data.get("bound_state_file")
        receptor_file = cleaned_data.get("receptor_file")
        receptor_residues = cleaned_data.get("receptor_residues")
        rdock_as = cleaned_data.get("rdock_as")
        rdock_grid = cleaned_data.get("rdock_grid")
        # Require either a bound or ground state SDF file
        if not ground_state_file and not bound_state_file:
            raise forms.ValidationError(
                {
                    "ground_state_file": "A series must have either a ground state or bound state file"
                }
            )

        # Confirm that the ground state SDF file has only one structure
        ground_state_mol = None
        if ground_state_file:
            tmp = get_tmp_filepath(
                ground_state_file, f"{self.instance.id}_ground_state.sdf"
            )
            suppl = Chem.SDMolSupplier(tmp, removeHs=False)
            ground_state_mol = suppl[0]
            if len(suppl) != 1:
                raise forms.ValidationError(
                    {
                        "ground_state_file": f"This file must contain exactly one structure, this one has {len(suppl)}!"
                    }
                )

        # Confirm that the bound state SDF file has only one structure
        bound_state_mol = None
        if bound_state_file:
            tmp = get_tmp_filepath(
                bound_state_file, f"{self.instance.id}_bound_state.sdf"
            )
            suppl = Chem.SDMolSupplier(tmp, removeHs=False)
            bound_state_mol = suppl[0]

            if len(suppl) != 1:
                raise forms.ValidationError(
                    {
                        "bound_state_file": f"This file must contain exactly one structure, this one has {len(suppl)}!"
                    }
                )

            # Remove any props associated with the mol - rDock fails on long prop values
            for prop in bound_state_mol.GetPropsAsDict().keys():
                bound_state_mol.ClearProp(prop)

            with Chem.SDWriter(tmp) as w:
                w.write(bound_state_mol)
            path = Path(tmp)
            original_name = bound_state_file.name
            f = path.open(mode="rb")
            bound_state_file = File(f, name=original_name)

        # If both ground and bound are specified, confirm that they are the same molecule
        if ground_state_mol and bound_state_mol:
            ground_state_inchi = Chem.MolToInchi(ground_state_mol)
            bound_state_inchi = Chem.MolToInchi(bound_state_mol)
            if ground_state_inchi != bound_state_inchi:
                raise forms.ValidationError(
                    {
                        "ground_state_file": "The ground state and bound state structures do not have the same Inchi. Double check the structures."
                    }
                )

        # Require binding pocket residues for a receptor file
        if receptor_file and not receptor_residues:
            raise forms.ValidationError(
                {
                    "receptor_residues": "A comma separated list of binding pocket residues (ex:'221,222,304') is required when a receptor file exists"
                }
            )

        # Require a receptor file for binding pocket residues
        if receptor_residues and not receptor_file:
            raise forms.ValidationError(
                {
                    "receptor_file": "A receptor file is required when receptor residues are specified"
                }
            )

        if (not rdock_grid and rdock_as) or (rdock_grid and not rdock_as):
            raise forms.ValidationError(
                {
                    "rdock_grid": "You must have both rdock_grid and rdock_as files",
                    "rdock_as": "You must have both rdock_grid and rdock_as files",
                }
            )

        return cleaned_data

    class Meta:
        model = Series
        exclude = []


class CompoundAdminForm(DeleteOldFilesModelForm):
    def save(self, commit=True):
        self.instance = super().save(commit=commit)
        # Call update_dn_id so the molecules in SDF files have their name updated to the current DN ID
        if "dn_id" in self.changed_data:
            self.instance.update_dn_id(self.instance.dn_id)
        return self.instance

    class Meta:
        model = Compound
        exclude = []


##############################################
#                Model Admins                #
##############################################


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ("code",)
    ordering = ["code"]


@admin.register(Series)
class SeriesAdmin(admin.ModelAdmin):
    list_display = ("dn_id", "name", "project", "active")
    form = SeriesAdminForm

    def get_queryset(self, request):
        return super().get_queryset(request).select_related("project")


@admin.register(Compound)
class CompoundAdmin(admin.ModelAdmin):
    list_display = ["virtual_id", "dn_id", "cos", "collections", "project"]
    list_filter = [
        ProjectFilter,
        CollectionFilter,
        CompoundPkFilter,
        CompoundDnIdFilter,
        CompoundOccurrenceFilter,
        UserFilter,
    ]
    raw_id_fields = ["series", "mmps"]
    form = CompoundAdminForm

    def get_queryset(self, request):
        parent_cos = (
            CompoundOccurrence.objects.filter(parent_co=None)
            .exclude(collection=None)
            .prefetch_related("collection_set")
            .order_by("pk")
        )
        compound_qs = (
            super()
            .get_queryset(request)
            .select_related("project")
            .prefetch_related(Prefetch("compoundoccurrence_set", queryset=parent_cos))
        )
        return compound_qs

    @admin.display(empty_value="", description="compound occurrences")
    def cos(self, comp):
        links = set()
        for co in comp.compoundoccurrence_set.all():
            link = '<a href="{}">{}</a>'.format(
                reverse("admin:main_compoundoccurrence_change", args=(co.pk,)), co.id
            )
            links.add(link)
        return mark_safe(", ".join(links))

    @admin.display(empty_value="", description="collections")
    def collections(self, comp):
        links = set()
        for co in comp.compoundoccurrence_set.all():
            for coll in co.collection_set.all():
                link = '<a href="{}">{}</a>'.format(
                    reverse("admin:main_collection_change", args=(coll.pk,)), coll.pk
                )
                links.add(link)
        return mark_safe(", ".join(links))


@admin.register(CompoundOccurrence)
class CompoundOccurrenceAdmin(admin.ModelAdmin):
    list_display = (
        "id",
        "owner",
        "compound_link",
        "parent_co_link",
        "gem_id",
        "dn_id",
        "collections",
    )
    list_filter = [
        ProjectFilter,
        CollectionFilter,
        CompoundPkFilter,
        CompoundDnIdFilter,
        CompoundOccurrenceFilter,
        UserFilter,
    ]
    raw_id_fields = ["compound", "parent_co"]

    def get_queryset(self, request):
        return (
            super()
            .get_queryset(request)
            .select_related("compound", "parent_co", "owner")
            .prefetch_related("collection_set")
        )

    @admin.display(empty_value="", description="compound")
    def compound_link(self, co):
        link = '<a href="{}">{}</a>'.format(
            reverse("admin:main_compound_change", args=(co.compound.pk,)),
            co.compound.id,
        )
        return mark_safe(link)

    @admin.display(empty_value="", description="parent CO")
    def parent_co_link(self, co):
        if co.parent_co:
            link = '<a href="{}">{}</a>'.format(
                reverse(
                    "admin:main_compoundoccurrence_change", args=(co.parent_co.pk,)
                ),
                co.parent_co.id,
            )
            return mark_safe(link)
        return ""

    @admin.display(empty_value="", description="DN ID")
    def dn_id(self, co):
        return co.compound.dn_id

    @admin.display(empty_value="", description="collections")
    def collections(self, co):
        colls = co.collection_set.all()
        links = []
        for coll in colls:
            link = '<a href="{}">{}</a>'.format(
                reverse("admin:main_collection_change", args=(coll.id,)), coll.id
            )
            links.append(link)
        return mark_safe(", ".join(links))


@admin.register(Collection)
class CollectionAdmin(admin.ModelAdmin):
    list_display = ("id", "owner", "project", "created_on", "compounds")
    list_filter = [ProjectFilter, CollectionFilter, UserFilter]
    raw_id_fields = ["compound_occurrences"]

    def get_queryset(self, request):
        co_qs = CompoundOccurrence.objects.all().select_related("compound")
        collection_qs = (
            super()
            .get_queryset(request)
            .select_related("project", "owner")
            .prefetch_related(Prefetch("compound_occurrences", queryset=co_qs))
        )
        return collection_qs

    @admin.display(empty_value="", description="compounds")
    def compounds(self, collection):
        links = set()
        for co in collection.compound_occurrences.all():
            link = '<a href="{}">{}</a>'.format(
                reverse("admin:main_compound_change", args=(co.compound.id,)),
                co.compound.id,
            )
            links.add(link)
        return mark_safe(", ".join(links))
