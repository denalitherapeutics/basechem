import logging

from django.contrib.auth.mixins import LoginRequiredMixin
from django.db.models import Case, IntegerField, Prefetch, Value, When
from django.http.response import HttpResponseRedirect, JsonResponse
from django.shortcuts import get_object_or_404
from django.urls import reverse
from django.views.generic import DetailView, ListView

from basechem.common.analytic import PAGE_VIEW, Analytic
from basechem.common.middleware import show_2d_warning
from basechem.main.constants import ALIGN, DOCK, ESP, MMP, PROPCALC, TORSION
from basechem.main.forms import ChooseReferenceForm
from basechem.main.mixins import AnalysisViewMixin
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, CompoundOccurrence
from basechem.main.tasks import get_analysis_results
from basechem.users.models import REQUIRED_EASTER_EGG_POINTS

logger = logging.getLogger("django")


class PropCalcView(LoginRequiredMixin, AnalysisViewMixin, ListView):
    """
    View to display property calculations
    """

    model = Compound
    context_object_name = "compounds"
    template_name = "main/propcalc/propcalc.html"

    def get_context_data(self):
        """
        Returns the context data which should contain the following keys:
            current_view = PROPCALC
            collection_id = id of the current collection
            properties = list of properties to display
            errors = list of errors that occurred during processing
            property_results = dict of {compound_id: list of properties}
        """
        context = super().get_context_data()
        collection = context["collection"]
        props = collection.get_propcalc_column_headers()
        if not props:
            logging.error(
                f"Collection {collection.id} is missing metadata for property calcs"
            )
        # Start or get results of group
        group_name = collection.get_propcalc_group_name()
        context["property_results"] = {}
        context["props"] = props
        context["group_name"] = group_name
        if props:
            get_analysis_results(collection, PROPCALC, group_name, {})
            # Log the page view to aws
            Analytic(
                PAGE_VIEW, collection.project.code, self.request.user, page=PROPCALC
            )
        context["loading_message"] = "Calculating properties..."
        return context

    def get_queryset(self):
        """
        Return a list of Compounds associated with the COs that are part of this collection
        :return: queryset of Compounds
        """
        collection = get_object_or_404(Collection, id=self.kwargs["collection_id"])
        return collection.compounds()


class LigandAlignView(LoginRequiredMixin, AnalysisViewMixin, DetailView):
    """
    View to display aligned ligands
    """

    model = Collection
    context_object_name = "collection"
    template_name = "main/align/align.html"

    slug_url_kwarg = "collection_id"
    slug_field = "id"

    def post(self, *args, **kwargs):
        collection = Collection.objects.get(id=kwargs.get("collection_id"))
        form = ChooseReferenceForm(self.request.POST, instance=collection)
        if form.is_valid():
            ref_string = form.cleaned_data["reference"]
            group_name = collection.get_align_group_name(ref_string)
            return JsonResponse(
                {
                    "errors": form.errors,
                    "redirect_url": reverse(ALIGN, args=[collection.id, group_name]),
                }
            )
        else:
            return JsonResponse({"errors": form.errors})

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        collection = context["collection"]
        context[
            "energytooltip"
        ] = "To show binding-relevant population ensemble, results are ordered by the four lowest RMSD confs, then any additional poses include the lowest energy conformers."
        group_name = self.kwargs.get("group_name")
        form = ChooseReferenceForm(instance=collection)
        if group_name:
            context["ref_string"] = group_name.split("_")[2]
            form.initial["reference"] = context["ref_string"]
            analysis_kwargs = {"ref_string": context["ref_string"]}
            context["failure"], _ = get_analysis_results(
                collection, ALIGN, group_name, analysis_kwargs
            )
            context["group_name"] = group_name
            context["loading_message"] = "Aligning compounds..."
            # Log the page view to aws
            Analytic(PAGE_VIEW, collection.project.code, self.request.user, page=ALIGN)

        context["form"] = form
        return context


class DockView(LoginRequiredMixin, AnalysisViewMixin, DetailView):
    """
    View to display docking results
    """

    model = Collection
    context_object_name = "collection"
    template_name = "main/dock/dock.html"

    slug_url_kwarg = "collection_id"
    slug_field = "id"

    def post(self, *args, **kwargs):
        collection = Collection.objects.get(id=kwargs.get("collection_id"))
        form = ChooseReferenceForm(
            self.request.POST, instance=collection, needs_structure=True
        )
        if form.is_valid():
            ref_string = form.cleaned_data["reference"]
            group_name = collection.get_dock_group_name(ref_string)
            return JsonResponse(
                {
                    "errors": form.errors,
                    "redirect_url": reverse(DOCK, args=[collection.id, group_name]),
                }
            )
        return JsonResponse({"errors": form.errors})

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        collection = context["collection"]
        group_name = self.kwargs.get("group_name")

        form = ChooseReferenceForm(instance=collection, needs_structure=True)
        if not group_name:
            # Start the default task if all compound occurrences in this collection ar assigned series
            # with mol2 files. Otherwise, don't start the default task.
            if not collection.compound_occurrences.filter(
                parent_co=None, compound__series__receptor_file_mol2=""
            ).exists():
                group_name = collection.get_dock_group_name("default")

        if group_name:
            context["ref_string"] = group_name.split("_")[2]
            form.initial["reference"] = context["ref_string"]
            analysis_kwargs = {"ref_string": context["ref_string"]}
            context["failure"], _ = get_analysis_results(
                collection, DOCK, group_name, analysis_kwargs
            )
            context["group_name"] = group_name
            context["loading_title"] = "Docking compounds..."
            # Log the page view to aws
            Analytic(PAGE_VIEW, collection.project.code, self.request.user, page=DOCK)
        context["form"] = form
        return context


class EspView(LoginRequiredMixin, AnalysisViewMixin, DetailView):
    """
    View to display ESP results
    """

    model = Collection
    context_object_name = "collection"
    template_name = "main/esp/esp.html"

    slug_url_kwarg = "collection_id"
    slug_field = "id"

    def post(self, *args, **kwargs):
        collection = Collection.objects.get(id=kwargs.get("collection_id"))
        form = ChooseReferenceForm(self.request.POST, instance=collection)
        if form.is_valid():
            ref_string = form.cleaned_data["reference"]
            group_name = collection.get_esp_group_name(ref_string)
            self.request.user.add_easter_egg_point(self.request, "ESP")
            return JsonResponse(
                {
                    "errors": form.errors,
                    "redirect_url": reverse(ESP, args=[collection.id, group_name]),
                }
            )
        else:
            return JsonResponse({"errors": form.errors})

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        collection = context["collection"]
        group_name = self.kwargs.get("group_name")

        form = ChooseReferenceForm(instance=collection)
        if not group_name:
            group_name = collection.get_esp_group_name("default")

        context["ref_string"] = group_name.split("_")[2]
        form.initial["reference"] = context["ref_string"]
        analysis_kwargs = {"ref_string": context["ref_string"]}
        context["failure"], _ = get_analysis_results(
            collection, ESP, group_name, analysis_kwargs
        )
        context["group_name"] = group_name
        context["loading_title"] = "Generating ESP maps..."
        # Log the page view to aws
        Analytic(PAGE_VIEW, collection.project.code, self.request.user, page=ESP)
        context["show_help_modal"] = int(
            self.request.user.easter_egg_points.get("ESP", 0)
            < REQUIRED_EASTER_EGG_POINTS["ESP"]
        )
        context["form"] = form

        return context


class TorsionView(LoginRequiredMixin, AnalysisViewMixin, DetailView):
    """
    View to display Torsion Scan results
    """

    model = Collection
    context_object_name = "collection"
    template_name = "main/torsion/torsion.html"

    slug_url_kwarg = "collection_id"
    slug_field = "id"

    def post(self, *args, **kwargs):
        collection = Collection.objects.get(id=kwargs.get("collection_id"))
        # Get data from form
        pioneer_pk = self.request.POST["pioneer"].split("-")[1]
        selected_atoms = self.request.POST["selected_atoms"]
        co = CompoundOccurrence.objects.get(pk=pioneer_pk)
        # Generate dihedral string and get group name
        dihedral_smarts = co.convert_atoms_to_smarts(selected_atoms)
        group_name = collection.get_torsion_group_name(pioneer_pk, selected_atoms)

        errors = {}
        if not co.clean_dihedrals(dihedral_smarts):
            errors["__all__"] = [
                "The 4 atoms you selected did not result in a valid dihedral to scan. Please try again."
            ]
            return JsonResponse(
                {
                    "errors": errors,
                    "redirect_url": reverse(TORSION, args=[collection.id]),
                }
            )

        return JsonResponse(
            {
                "errors": errors,
                "redirect_url": reverse(TORSION, args=[collection.id, group_name]),
            }
        )

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        collection = context["collection"]
        context[
            "energytooltip"
        ] = "Conformation with the lowest QM energy is 0. Relative energy of all other dihedrals are calculated in comparison to that lowest energy dihedral. The input conformation relative energy is derived from the QM energy of the conformation with the dihedral that is most similar to the input dihedral (e.g., input NCCO is 181, relative energy is approximated from conformation with NCCO of 180)."
        group_name = self.kwargs.get("group_name")
        context["pioneer_options"] = (
            context["compound_occurrences"]
            .order_by("compound", "pk")
            .distinct("compound")
        )

        if group_name:
            pioneer_pk = int(group_name.split("_")[2])
            context["pioneer"] = CompoundOccurrence.objects.get(pk=pioneer_pk)
            context["pioneer_dihedral_atoms"] = group_name.split("_")[3].replace(
                "-", ","
            )
            context["pioneer_dihedral_smarts"] = context[
                "pioneer"
            ].convert_atoms_to_smarts(context["pioneer_dihedral_atoms"])
            context["failure"], _ = get_analysis_results(
                collection,
                TORSION,
                group_name,
                {
                    "pioneer_pk": pioneer_pk,
                    "pioneer_dihedral_atoms": context["pioneer_dihedral_atoms"],
                },
            )
            context["group_name"] = group_name
            context["loading_title"] = "Scanning torsion dihedrals..."
            context[
                "loading_message"
            ] = "This will probably take a long time...come back later! You will get an email when the QM job completes."
            # Add 1 to the pioneer dihedral atom indices to get the atom numbers
            # because atom indices are 0 indexed and atom numbers are 1 indexed
            context["pioneer_dihedral_atom_nums"] = ", ".join(
                map(lambda x: str(int(x) + 1), group_name.split("_")[3].split("-"))
            )
            # Log the page view to aws
            Analytic(
                PAGE_VIEW, collection.project.code, self.request.user, page=TORSION
            )
            self.request.user.add_easter_egg_point(self.request, "Torsion")
            context["show_help_modal"] = (
                self.request.user.easter_egg_points["Torsion"]
                < REQUIRED_EASTER_EGG_POINTS["Torsion"]
            )
        else:
            # Show users a warning if compounds are in 2D
            show_2d_warning(self.request, collection, TORSION)

        return context


class DTXMMPView(LoginRequiredMixin, AnalysisViewMixin, DetailView):
    """
    View to display MMPs found in Dotmatics for each Compound in a Collection
    """

    model = Collection
    context_object_name = "collection"
    template_name = "main/mmps/dtx_mmps.html"

    slug_url_kwarg = "collection_id"
    slug_field = "id"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        collection = context["collection"]
        compounds = collection.compounds()
        compound_pks = compounds.values_list("pk", flat=True)
        dtx_cpds = Compound.objects.exclude(dn_id="").annotate(
            assayed_comp=Case(
                When(pk__in=compound_pks, then=Value(1)),
                default=Value(0),
                output_field=IntegerField(),
            )
        )
        compound_qs = compounds.prefetch_related(
            Prefetch("mmps", queryset=dtx_cpds, to_attr="dtx_mmps")
        )

        compounds_list = []
        for c in compound_qs:
            c.assayed_comp = True
            s_mmps = sorted(
                list(c.dtx_mmps),
                key=lambda x: (
                    x.assayed_comp,
                    float(x.metadata.get("similarity", {}).get(str(c.pk), 0)),
                ),
                reverse=True,
            )
            c_tup = (c, s_mmps[:20])
            compounds_list.append(c_tup)
        context["compounds"] = compounds_list

        self.request.user.add_easter_egg_point(self.request, "Assay MMPs")
        return context


class MMPSearchView(LoginRequiredMixin, DetailView):
    """
    View to display MMPs found from MMP Analysis Form
    """

    model = Collection
    context_object_name = "collection"
    template_name = "main/mmps/mmps.html"

    slug_url_kwarg = "collection_id"
    slug_field = "id"

    def get(self, request, *args, **kwargs):
        collection_id = self.kwargs.get("collection_id")
        group_name = self.kwargs.get("group_name")
        if not group_name:
            collection = Collection.objects.get(pk=collection_id)
            co = collection.get_cos_for_analysis(MMP)[0]
            group_name = collection.get_mmp_group_name(co.compound.pk)
            redirect = reverse(
                MMP, kwargs={"collection_id": collection_id, "group_name": group_name}
            )
            return HttpResponseRedirect(redirect)
        return super().get(request, *args, **kwargs)

    def get_context_data(self, *args, **kwargs):
        context = super().get_context_data(*args, **kwargs)
        collection = context["collection"]
        group_name = self.kwargs["group_name"]
        comp_pk = group_name.split("_")[-1]
        context["co"] = collection.compound_occurrences.filter(
            compound__pk=comp_pk
        ).first()
        _, _ = get_analysis_results(collection, MMP, group_name, {})
        context["group_name"] = group_name
        context["loading_message"] = "Generating MMPs..."
        # Log the page view to aws
        Analytic(PAGE_VIEW, collection.project.code, self.request.user, page=MMP)
        return context
