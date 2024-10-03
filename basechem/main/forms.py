import datetime
import json
import os
import re

from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from rdkit import Chem

from basechem.common.rdkit_utils import RDKitWrappers
from basechem.main.constants import *
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import CompoundOccurrence, Series
from basechem.main.models.project_models import Project
from basechem.main.utils import strip_series_from_conf_id
from basechem.mni_common.forms.fields import NoValidateMultipleChoiceField
from basechem.mni_common.forms.forms import CrispyForm, CrispyModelForm


def clean_moltext(moltext, invalid_message="Are you sure this is valid moltext?"):
    """
    Clean moltext fields by fragmenting them into the component molecules and constructing
    RDKit mol objects.
    :param moltext: a string of moltext
    :param invalid_message: a string, the error message that should be shown to the user if the moltext is invalid
    :returns: a list of tuples of the form (mol, is_2d) where mol is an rdkit mol object and is_2d is a boolean
    """
    if moltext:
        # Going to leave strictParsing as False to be more lenient for now
        m = Chem.MolFromMolBlock(moltext, strictParsing=False, removeHs=False)
        if m:
            return RDKitWrappers.clean_mol_object(m)
        else:
            raise ValidationError(invalid_message)
    return []


class MMPKetcherWidget(forms.widgets.Textarea):
    """
    Extended Ketcher Widget that includes buttons to set a variable region for MMP analysis
    """

    template_name = "main/components/widgets/mmp_ketcher_widget.html"


class MMPIntakeField(forms.CharField):
    widget = MMPKetcherWidget
    strip = False

    def _convert_indices(self, ketcher_indices, index_map):
        """
        Converts indices from the ketcher widget to corresponding moltext indices. This is
        needed because when fragments are deleted in the ketcher widget, the atom indices don't reset
        which leads to the indices not matching the moltext atom index (Ex: if a user draws
        benzene, deletes it, and redraws it, it will have atom indices 6-11 when it should have indices 0-5).
        :param ketcher_indices: a list of atom or bond indices from the ketcher widget
        :param index_map: a dictionary of the form {ketcher_index: moltext_index}. This dictionary
            is maintained in the frontend in mmp_ketcher.js
        :return: a list of moltext atom indices that correspond to ketcher_indices
        """
        return [index_map[str(idx)] for idx in ketcher_indices]

    def _clean_moltext(self, moltext):
        """
        Called in the clean method, this helper cleans the moltext provided from the front end,
        additionally confirming that only one molecule was drawn
        :param moltext: a string of moltext
        :return: an RDKit mol object
        """
        mols = clean_moltext(moltext, "Are you sure this is a valid molecule?")
        if not mols:
            raise forms.ValidationError(
                "You must sketch a compound and highlight a variable region"
            )
        elif len(mols) > 1:
            raise forms.ValidationError("You may only draw one molecule at a time")
        mol = mols[0][0]
        # Check that the mol has atoms (a blank sketcher returns valid moltext)
        if mol.GetNumAtoms() == 0:
            raise forms.ValidationError(
                "You must sketch a compound and highlight a variable region"
            )
        return mol

    def _clean_highlighted_region(self, widget_data, mol):
        """
        Parse the widget data into SMILES strings for the constant and variable regions and
        raise ValidationErrors if the highlighted region is invalid.
        :param widget_data: a dictionary of 'widgetData' from mmp_ketcher.js that takes the form
            {
                "moltext": moltext,
                "variable:{"atoms":[], "bonds":[]},
                "maps":{"atoms":{}, "bonds":{}}
            }
        :param mol: a mol object, the full molecule that was drawn
        :return: a tuple (variable_smiles, constant_smiles), the SMILES strings of the variable and constant regions
        """
        # Get the atom indices of the variable region
        variable_atom_indices = self._convert_indices(
            widget_data["variable"]["atoms"], widget_data["maps"]["atoms"]
        )
        if not variable_atom_indices:
            raise forms.ValidationError("You must highlight a variable region")

        # Check that no rings are split by the variable region
        for ring in mol.GetRingInfo().AtomRings():
            ring_atoms_selected = [
                atom_index in variable_atom_indices for atom_index in ring
            ]
            if any(ring_atoms_selected) and not all(ring_atoms_selected):
                raise forms.ValidationError(
                    "The variable region cannot include a partial ring"
                )

        # Find the bond that separates the variable and constant regions
        bond_index = None
        for atom_i in variable_atom_indices:
            atom = mol.GetAtomWithIdx(atom_i)
            for bond in atom.GetBonds():
                bonded_atom_i = bond.GetBeginAtomIdx()
                if bonded_atom_i == atom_i:
                    bonded_atom_i = bond.GetEndAtomIdx()
                if bonded_atom_i not in variable_atom_indices:
                    if bond_index != None:
                        raise forms.ValidationError(
                            f"The variable region cannot be disconnected. Please highlight a connected region."
                        )
                    bond_index = bond.GetIdx()

        if bond_index == None:
            raise forms.ValidationError(
                f"The variable region cannot be the entire molecule. Please highlight a substructure."
            )

        # Fragment the molecule on the bond separating the variable and constant regions
        frag_mol = Chem.FragmentOnBonds(mol, (bond_index,))
        frag_atom_indices = []
        frags = Chem.GetMolFrags(
            frag_mol, asMols=True, fragsMolAtomMapping=frag_atom_indices
        )
        if len(frags) > 2:
            raise forms.ValidationError(
                f"The variable region cannot be disconnected. Please highlight a connected region."
            )

        c_smiles, v_smiles = None, None
        for frag, atom_indices in zip(frags, frag_atom_indices):
            smiles = re.sub(r"\[\d+\*\]", "*", Chem.MolToSmiles(frag))
            if any(
                [atom_index in variable_atom_indices for atom_index in atom_indices]
            ):
                v_smiles = smiles
            else:
                c_smiles = smiles

        if not c_smiles and not v_smiles:
            raise forms.ValidationError(
                f"Could not parse the variable and constant regions. Try again or contact MnI."
            )
        return v_smiles, c_smiles

    def clean(self, widget_data):
        """
        Given a json object with widget data from the front end, return a dictionary with the
        RDKit mol object and the SMILES strings for the variable and constant regions.
        :param widget_data: json of 'widgetData' from mmp_ketcher.js that takes the form
            {
                "moltext": moltext,
                "variable:{"atoms":[], "bonds":[]},
                "maps":{"atoms":{}, "bonds":{}}
            }
        :return: a dictionary with keys 'mol', 'variable_smiles', and 'constant_smiles' where
            'mol' is an RDKit mol object and the variable and constant smiles are smiles
            strings defining the constant and variable regions for the MMP search.
        """
        if widget_data:
            widget_data = json.loads(widget_data)
        else:
            raise forms.ValidationError(
                "You must sketch a compound and highlight a variable region"
            )
        mol = self._clean_moltext(widget_data["moltext"])
        variable_smiles, constant_smiles = self._clean_highlighted_region(
            widget_data, mol
        )
        return {
            "mol": mol,
            "variable_smiles": variable_smiles,
            "constant_smiles": constant_smiles,
        }


class KetcherWidget(forms.widgets.Textarea):
    template_name = "main/components/widgets/ketcher_widget.html"


class KetcherField(forms.CharField):
    widget = KetcherWidget
    strip = False

    def clean(self, value):
        return clean_moltext(
            value, invalid_message="Are you sure this is a valid molecule?"
        )


class CompoundIntakeForm(CrispyModelForm):
    """
    Form to read in sdf file or mol text
    """

    compound_upload_fields = ["sketcher", "upload_file", "moltext"]
    collection = forms.ModelChoiceField(
        queryset=Collection.objects.all(),
        empty_label="Make New Collection",
        required=False,
    )
    name = forms.CharField(
        max_length=50,
        required=False,
        help_text="Optional: Give this collection a name to remember it by",
    )
    upload_file = forms.FileField(label="Choose an SDF file", required=False)
    moltext = forms.CharField(
        label="Paste MOL text", widget=forms.Textarea, required=False, strip=False
    )
    sketcher = KetcherField(required=False)
    project = forms.ModelChoiceField(queryset=Project.objects.all(), required=False)

    def __init__(self, current_user, collection_id=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        two_months_ago = datetime.datetime.now().date() - datetime.timedelta(days=90)
        self.available_collections = Collection.objects.filter(
            owner=current_user, created_on__gt=two_months_ago
        ).order_by("-pk")
        initial_collection = None
        if collection_id:
            try:
                initial_collection = self.available_collections.get(id=collection_id)
            except Exception:
                # If chosen collection is too old or not owned by the current user, add it to the queryset
                old_collection = Collection.objects.filter(id=collection_id)
                self.available_collections = self.available_collections | old_collection
                initial_collection = old_collection.first()

        self.fields["collection"].queryset = self.available_collections
        self.initial["collection"] = initial_collection

    class Meta:
        model = Collection
        fields = [
            "collection",
            "name",
            "project",
            "sketcher",
            "upload_file",
            "moltext",
        ]

    def clean(self):
        """
        Make sure either an sdf file or mol text are uploaded
        """
        cleaned_data = super().clean()
        file = self.files
        file_data = cleaned_data.get("upload_file")
        mol_data = cleaned_data.get("moltext")
        sketcher_data = cleaned_data.get("sketcher")
        collection = cleaned_data.get("collection")
        if not collection:
            if not cleaned_data.get("project"):
                self.add_error("project", "This field is required.")

            if not file and not mol_data and not sketcher_data:
                errors = any(
                    [
                        bool(self.errors.get(field))
                        for field in ["upload_file", "moltext", "sketcher"]
                    ]
                )
                if not errors:
                    self.add_error("upload_file", "")
                    self.add_error("moltext", "")
                    self.add_error("sketcher", "")
                    raise forms.ValidationError(
                        "You must upload an SDF file, paste in moltext, or sketch a compound"
                    )

            comps = []
            if type(file_data) == list:
                comps.extend(file_data)
            if type(mol_data) == list:
                comps.extend(mol_data)
            if type(sketcher_data) == list:
                comps.extend(sketcher_data)
            if len(comps) == 0:
                raise forms.ValidationError(
                    "Your uploaded data did not return any valid compounds"
                )

        return cleaned_data

    def clean_upload_file(self):
        """
        Check that the uploaded file has an sdf extension
        """
        data = self.cleaned_data["upload_file"]
        if data:
            filename, file_extension = os.path.splitext(data.name)
            if file_extension != ".sdf":
                self.add_error("upload_file", "Must be an sdf file.")
                return data
            data = Collection.handle_sdf_upload(data)
        return data

    def clean_moltext(self):
        """
        Check that the uploaded text is a mol object
        """
        data = self.cleaned_data["moltext"]
        return clean_moltext(data)


######################
### PROPCALC FORMS ###
######################


class PropCalcForm(CompoundIntakeForm):
    title = "Calculate Properties"
    tab_title = "PropCalc"
    directions = (
        "Please choose a previous collection, upload an sdf file, or paste in mol text. Then select the properties you would like displayed. "
        "If you receive an error, please contact the admins. "
    )

    COUNTS_PROPERTIES = [
        (ACCEPTORS, "Acceptors: Number of N and O atoms"),
        (AROMATIC_RINGS, "Aromatic Rings: Number of aromatic rings"),
        (CHARGE, "Charge: Formal charge at pH 7.4"),
        (
            DONORS,
            "Donors: Number of OH/NH groups in the neutral molecule, NH2 counts as two",
        ),
        (FRACTIONCSP3, "Fraction C SP3: Fraction of SP3 carbon atoms"),
        (RINGS, "Rings: Number of rings"),
        (ROTATABLEBONDS, "Rotatable Bonds: Number of rotatable bonds"),
        (HEAVYATOMS, "Heavy Atoms: Number of heavy atoms"),
    ]

    PHYSIOCHEMICAL_PROPERTIES = [
        (MW, "MW: Molecular Weight of the neutral molecule"),
        (SOLUBILITYINDEX, "Solubility Index: Aromatic Rings + MolLogP"),
        (TPSA, "TPSA: Topological polar surface area"),
        (CLOGP, "Calculated LogP: LogP is measured at pH where molecule is neutral"),
    ]
    initial_physiochemical = [MW, TPSA, HLM]
    if settings.INDUCTIVE_BIO_ENABLED:
        PHYSIOCHEMICAL_PROPERTIES.extend(
            [
                (ALOGD, "LogD: Predicted LodD using the InductiveBio ML model"),
                (RLM, "RLM: Predicted RLM stability using the InductiveBio ML model"),
                (HLM, "HLM: Predicted HLM stability using the InductiveBio ML model"),
            ]
        )
        initial_physiochemical.append(ALOGD)

    counts = forms.MultipleChoiceField(
        required=False, widget=forms.CheckboxSelectMultiple, choices=COUNTS_PROPERTIES
    )

    physiochemical = forms.MultipleChoiceField(
        required=False,
        widget=forms.CheckboxSelectMultiple,
        choices=PHYSIOCHEMICAL_PROPERTIES,
        initial=initial_physiochemical,
    )


######################
##### ALIGN FORMS ####
######################


class LigandAlignSubmitForm(CompoundIntakeForm):
    """
    Form that creates a Collection to be used for ligand alignment
    """

    title = "Align Ligands"
    tab_title = "Align"
    directions = (
        "Please choose a previous collection, upload an sdf file, or paste in mol text to upload your compounds for ligand alignment. "
        "If you receive an error, please contact the admins. "
    )


class ChooseReferenceForm(forms.ModelForm):
    """
    Form used to choose a reference for docking or align
    """

    reference = forms.ChoiceField(label="Reference Ligand")
    needs_structure = False

    def __init__(self, *args, **kwargs):
        self.needs_structure = kwargs.pop("needs_structure", False)
        super().__init__(*args, **kwargs)
        reference_choices = [("default", "Default - use assigned series")]
        if self.needs_structure:
            series = Series.objects.filter(
                project=self.instance.project, active=True
            ).exclude(receptor_file_mol2="")
        else:
            series = Series.objects.filter(project=self.instance.project, active=True)
        reference_choices.extend([(f"s-{s.id}", f"Series: {s.name}") for s in series])

        self.fields["reference"].choices = reference_choices

    def is_valid(self, *args, **kwargs):
        is_valid = super().is_valid()
        if (
            is_valid
            and self.needs_structure
            and self.cleaned_data["reference"] == "default"
        ):
            # Display an error if a CO in this collection has an assigned series without a receptor
            invalid_cos = self.instance.compound_occurrences.filter(
                parent_co=None, compound__series__receptor_file_mol2=""
            )
            if invalid_cos.exists():
                invalid_series = ", ".join(
                    invalid_cos.distinct("compound__series").values_list(
                        "compound__series__name", flat=True
                    )
                )
                message = f"At least one of the assigned Series don't have the necessary files to run Dock. Please select a reference from this dropdown. Offending Series: {invalid_series}"
                self.add_error(None, message)
                is_valid = False
        return is_valid

    class Meta:
        model = Collection
        fields = ("reference",)


class DockSubmitForm(CompoundIntakeForm):
    """
    Form that creates a Collection to be used for docking
    """

    title = "Dock Compounds"
    tab_title = "Dock"
    directions = (
        "Please choose a previous Collection, upload an sdf file, or paste in mol text to upload your compounds for docking. "
        "If you receive an error, please contact the admins. "
    )

    def __init__(self, current_user, collection_id, *args, **kwargs):
        super().__init__(current_user, collection_id, *args, **kwargs)
        project_pks_w_structure = [
            p.pk for p in Project.objects.all() if p.structure_available
        ]
        projects_w_structure = Project.objects.filter(pk__in=project_pks_w_structure)
        self.fields["project"].queryset = projects_w_structure
        self.available_collections = self.available_collections.filter(
            project__in=projects_w_structure
        )
        self.fields["collection"].queryset = self.available_collections
        try:
            self.initial["collection"] = self.available_collections.get(
                id=collection_id
            )
        except Exception:
            # Chosen collection is not structurally enabled
            self.initial["collection"] = None


class EspSubmitForm(CompoundIntakeForm):
    """
    Form that creates a Collection to be used for esp generation
    """

    title = "Generate ESP"
    tab_title = "ESP"
    directions = (
        "Please upload an sdf file or paste in mol text to upload your compounds to generate ESP maps. "
        "If you receive an error, please contact the admins. "
    )


class TorsionSubmitForm(CompoundIntakeForm):
    """
    Form that creates a Collection to be used for torsion scans
    """

    title = "Torsion Scan"
    tab_title = "Torsion"
    directions = (
        "Please upload an sdf file or paste in mol text to upload your compounds to run torsion scans maps. "
        "If you have a preferred pose you would like to use for torsion scanning, make sure you upload it as a 3D structure! "
        "If you upload a 2D structure, Basechem will generate a 3D conformation for you to use instead. "
        "If you receive an error, please try again later or contact MnI. "
    )


class MMPSubmitForm(CrispyModelForm):
    """
    Form that creates a Collection to be used for MMP Analysis
    """

    title = "MMP Analysis"
    tab_title = "MMP"
    directions = "Sketch a molecule and highlight a variable region to generate matched molecular pairs for this molecule."

    compound_upload_fields = ["sketcher"]
    name = forms.CharField(
        max_length=50,
        required=False,
        help_text="Optional: Give this collection a name to remember it by",
    )
    sketcher = MMPIntakeField(label="Sketch molecules", required=False)
    project = forms.ModelChoiceField(queryset=Project.objects.all(), required=True)

    def __init__(self, current_user, collection_id=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class Meta:
        model = Collection
        fields = [
            "name",
            "project",
            "sketcher",
        ]

    def clean(self):
        cleaned_data = super().clean()
        if "sketcher" in self.errors:
            # Errors are more visible at the top of the page than under the sketcher
            self.add_error(None, self.errors.pop("sketcher"))
        return cleaned_data


class HikeForm(CrispyModelForm):
    """
    This form is displayed in a modal to allow users to 'hike' from one analysis to another
    """

    directions = "Use this Collection for another analysis using the options below!"

    analysis = forms.ChoiceField(
        choices=[
            (PROPCALC, "Calculate Properties"),
            (ALIGN, "Align Ligands"),
            (DOCK, "Dock Compounds"),
            (ESP, "Generate ESP Maps"),
            (TORSION, "Explore Torsion Scans"),
        ]
    )

    def __init__(self, current_view=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        all_choices = self.fields["analysis"].choices
        disabled_choices = [current_view]
        if not self.instance.project.structure_available:
            disabled_choices.append(DOCK)
        analysis_choices = [c for c in all_choices if c[0] not in disabled_choices]
        self.fields["analysis"].choices = analysis_choices

    class Meta:
        model = Collection
        fields = ["analysis"]


class SaveGemsForm(CrispyModelForm):

    gems = NoValidateMultipleChoiceField(
        choices=[], widget=forms.CheckboxSelectMultiple, required=False, initial=[]
    )
    other_gems = NoValidateMultipleChoiceField(
        choices=[],
        widget=forms.CheckboxSelectMultiple,
        required=False,
        initial=[],
        help_text="If you remove a previously saved gem, you won't be able to use it unless you save it from the original analysis again.",
        label="Gems from other hikes",
    )

    def __init__(self, collection, confs, keys=None, *args, **kwargs):
        """
        :param collection: the Collection object which is having conformers saved to it
        :param confs: a dictionary of conformer data for this CompoundOccurrence (usually this
            is directly from the task's result dictionary)
        :param keys: a list of keys in the nested `confs` dict that should be used to display
            the conformers in the form. For example, if keys is ["scoreA", "scoreB"], then the
            conformers in the form will display as "conf_id: scoreA - scoreB", where scoreA
            and scoreB are the scores for that conf_id
        """
        super().__init__(*args, **kwargs)
        self.confs = confs
        saved_child_cos = self.instance.get_child_cos(collection=collection)
        # Construct choices for `conformers`
        choices = []
        for conf_id, conf in confs.items():
            conf_display = strip_series_from_conf_id(conf_id)
            if keys:
                conf_display += f': {", ".join(conf[key] for key in keys)}'
                self.fields[
                    "gems"
                ].help_text = f"IDs are formatted as --> conf_id: {', '.join(keys)}"
            choices.append((conf_id, conf_display))
        self.fields["gems"].choices = choices
        self.initial["gems"] = list(saved_child_cos.values_list("gem_id", flat=True))

        # Construct choices for `other_gems`
        other_gems_choices = []
        other_gems_initial = []
        for co in saved_child_cos:
            if co.gem_id not in self.confs:
                other_gems_choices.append(
                    (co.gem_id, f"{co.saved_from} {co.display_gem_with_series()}")
                )
                other_gems_initial.append(co.gem_id)
        self.fields["other_gems"].choices = other_gems_choices
        self.initial["other_gems"] = other_gems_initial
        if not other_gems_choices:
            help_text = "If there are no options here, there are no other gems saved for this compound."
            self.fields["other_gems"].help_text = help_text

        # Update ID so that each form has unique IDs
        self.fields["gems"].widget.attrs.update(
            {"id": f"save-gems-form-{self.instance.pk}_gems"}
        )
        self.fields["other_gems"].widget.attrs.update(
            {"id": f"save-gems-form-{self.instance.pk}_other_gems"}
        )

    def is_valid(self, *args, **kwargs):
        is_valid = super().is_valid(*args, **kwargs)
        if is_valid:
            if (
                len(self.cleaned_data["gems"]) + len(self.cleaned_data["other_gems"])
                > MAX_GEMS_PER_CO
            ):
                self.add_error(
                    None, f"You cannot save more than {MAX_GEMS_PER_CO} gems at a time."
                )
                is_valid = False
        return is_valid

    class Meta:
        model = CompoundOccurrence
        fields = ["gems", "other_gems"]


class SubstructureSearchForm(CrispyForm):

    search_type = forms.ChoiceField(
        choices=[("sss", "Substructure"), ("exact", "Exact")], required=True
    )
    sketcher = KetcherField(required=True)

    def clean_sketcher(self):
        data = self.cleaned_data["sketcher"]
        if not data:
            raise ValidationError("You must provide a valid structure")
        elif len(data) > 1:
            raise ValidationError("You may only search for one structure at a time")
        else:
            data = Chem.MolToSmiles(data[0][0])
            if "*" in data:
                raise ValidationError("This does not yet support wildcard atoms")

        return data


class AddCompoundForm(CrispyModelForm):

    cos_added = forms.CharField(required=False, widget=forms.HiddenInput)
    sketcher = KetcherField(label="Sketch molecule(s)", required=True)

    def __init__(self, analysis=None, moltext=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.analysis = analysis

    def is_valid(self):
        is_valid = super().is_valid()
        if is_valid:
            cos_added = self.instance.handle_romols(self.cleaned_data["sketcher"])
            self.cleaned_data["cos_added"] = cos_added
            if not cos_added:
                is_valid = False
                self.add_error(None, "This molecule is already in the collection")
        return is_valid

    class Meta:
        model = Collection
        fields = ["sketcher"]
