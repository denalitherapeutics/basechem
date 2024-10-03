import logging
import os
from wsgiref.util import FileWrapper

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core import signing
from django.core.mail import mail_admins
from django.db.models import Prefetch
from django.http import HttpResponse, HttpResponseRedirect
from django.http.response import JsonResponse
from django.shortcuts import render
from django.urls import reverse, reverse_lazy
from django.views.decorators.clickjacking import xframe_options_sameorigin
from django.views.generic import FormView, ListView, TemplateView, UpdateView, View
from django_q.models import Task
from django_q.tasks import async_task, fetch_group, result
from rdkit import Chem

from basechem.common.analytic import PAGE_VIEW, Analytic
from basechem.common.constants import ADMIN_FAILURE
from basechem.main.constants import ALIGN, AUTO, DOCK, ESP, MMP, PROPCALC, TORSION
from basechem.main.forms import (
    AddCompoundForm,
    CompoundIntakeForm,
    DockSubmitForm,
    EspSubmitForm,
    HikeForm,
    LigandAlignSubmitForm,
    MMPSubmitForm,
    PropCalcForm,
    SaveGemsForm,
    SubstructureSearchForm,
    TorsionSubmitForm,
)
from basechem.main.mixins import SaveGemsMixin
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, CompoundOccurrence
from basechem.main.models.project_models import Project
from basechem.mni_common.rdkit_utils import moltext_to_svg

logger = logging.getLogger("django")


class HomePageView(LoginRequiredMixin, ListView):
    """
    Home page to choose what actions to perform
    """

    model = Project
    context_object_name = "projects"
    template_name = "main/homepage.html"

    def get_queryset(self):
        """
        Only show filtered Projects
        """
        Analytic(PAGE_VIEW, None, self.request.user, page="home")
        qs = super().get_queryset()
        return qs.exclude(code=AUTO)


class SubmitCompoundsView(LoginRequiredMixin, FormView):
    """
    View to collect compounds from sdf or pasted mol
    """

    model = Collection
    template_name = "main/submitcompounds.html"
    success_url = reverse_lazy("homepage")

    def get(self, request, *args, **kwargs):
        """
        If the new_collection_async task has already completed without errors, redirect to the next page
        """
        nextview = self.kwargs["nextview"]
        task_id = self.kwargs.get("task_id")
        if task_id:
            # Check task results
            task_result = result(task_id)
            if task_result != None and not task_result.get("error"):
                # New collection task is complete and successful
                if nextview == "homepage":
                    redirect = reverse("homepage")
                else:
                    redirect = reverse(
                        self.kwargs["nextview"],
                        kwargs={"collection_id": self.kwargs["collection_id"]},
                    )
                return HttpResponseRedirect(redirect)
        return super().get(request, *args, **kwargs)

    def get_form_class(self):
        """
        Dynamically sets the form class based off the url arg
        """
        nextview = self.kwargs["nextview"]
        if nextview == PROPCALC:
            return PropCalcForm
        elif nextview == ALIGN:
            return LigandAlignSubmitForm
        elif nextview == DOCK:
            return DockSubmitForm
        elif nextview == ESP:
            return EspSubmitForm
        elif nextview == TORSION:
            return TorsionSubmitForm
        elif nextview == MMP:
            return MMPSubmitForm
        else:
            return CompoundIntakeForm

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["current_user"] = self.request.user
        kwargs["collection_id"] = self.kwargs.get("collection_id")
        return kwargs

    def get_context_data(self, **kwargs):
        """
        Sets additional context data from the form
        """
        context = super(SubmitCompoundsView, self).get_context_data(**kwargs)
        context["loading_title"] = "Processing compounds..."
        nextview = self.kwargs["nextview"]
        task_id = self.kwargs.get("task_id")

        # Here we will swap out subclassed forms based on what analysis the user
        # wants to get to next to display additional metadata fields when needed
        if nextview != "homepage":
            context["nextview"] = nextview
            context["title"] = context["form"].title
            context["tab_title"] = context["form"].tab_title
            context["directions"] = context["form"].directions

        if task_id:
            task_result = result(task_id)
            if task_result == None:
                # New collection task is still going
                context["task_id"] = task_id
            elif task_result.get("error"):
                # New collection task is complete, but there's an error, send admins an email
                message = f'New collection task {task_id} for collection {self.kwargs.get("collection_id")} failed with "{task_result["error"]}"'
                mail_admins(ADMIN_FAILURE, message)
                # And show users a message

                if "DtxException" in task_result["error"]:
                    messages.error(
                        self.request,
                        f"There was a Dotmatics error while processing this collection: '{task_result['error'].split(': ')[1]}'",
                    )
                else:
                    messages.error(
                        self.request,
                        "There was an error while processing this collection. Please double check your uploaded structures and reach out to MnI.",
                    )
        return context

    def update_metadata(self, collection, form):
        """
        Updates relevant metadata based on the submitted form and deletes any existing
        propcalc task if the next view is propcalc
        :param collection: the collection being processed
        :param form: the CompoundIntake form being processed
        """
        metadata = {}
        if self.get_form_class() == PropCalcForm:
            metadata["props_to_show"] = form.cleaned_data.get(
                "counts", []
            ) + form.cleaned_data.get("physiochemical", [])

        if self.get_form_class() == MMPSubmitForm:
            metadata["mmp_analysis"] = {
                "constant_smiles": form.cleaned_data["sketcher"]["constant_smiles"],
                "variable_smiles": form.cleaned_data["sketcher"]["variable_smiles"],
            }
        collection.metadata = metadata
        collection.save()

        if self.kwargs["nextview"] == PROPCALC:
            # Delete existing propcalc tasks if they exist
            tasks = fetch_group(collection.get_propcalc_group_name())
            if tasks:
                tasks.delete()

    def process_new_collection(self, form):
        """
        Creates a new collection when an existing collection was not selected from
        the dropdown menu. Starts an async task to process compounds and compound occurrences
        :param form: the CompoundIntakeForm being processed
        :returns: the url path that this collection should redirect to
        """
        # Extract romols from form data
        project = form.cleaned_data["project"]
        created_by = self.request.user
        extracted_romol = []
        for mol_field in form.compound_upload_fields:
            data = form.cleaned_data[mol_field]
            if data:
                if isinstance(data, dict):
                    extracted_romol.append((data["mol"], True))
                else:
                    extracted_romol.extend(data)

        # Create a new collection for this submission
        collection = form.save(commit=False)
        collection.owner = created_by
        collection.project = project
        collection.save()
        self.update_metadata(collection, form)
        # Save romols to sdf to maintain uploaded properties
        collection.save_romols_to_sdf(extracted_romol)
        task_id = async_task(
            collection.new_collection_async,
            romols=extracted_romol,
            task_name=f"new-collection_{collection.pk}_{self.kwargs['nextview']}",
            group="new-collection",
            cluster="fast",
        )

        return reverse("submit", args=[self.kwargs["nextview"], collection.id, task_id])

    def form_valid(self, form):
        response = super().form_valid(form)
        collection = form.cleaned_data.get("collection")
        if collection:
            # Collection already made, go straight to next view
            self.update_metadata(collection, form)
            if self.kwargs["nextview"] == "homepage":
                redirect_url = reverse("homepage")
            else:
                redirect_url = reverse(
                    self.kwargs["nextview"], kwargs={"collection_id": collection.pk}
                )
        else:
            redirect_url = self.process_new_collection(form)

        return HttpResponseRedirect(redirect_url)

    @staticmethod
    def handle_exception(exception):
        """
        Returns an error message for commonly seen exceptions
        :param exception: a python exception that was raised during the view and caught by the middleware
        :return: a string that is a specific error message to be displayed
        """
        if type(exception) == KeyError:
            return (
                "We encountered an error creating this collection. "
                "There might be an issue with your compound occurrence."
            )
        else:
            return "We encountered an error creating this collection."


class DownloadResultFileView(LoginRequiredMixin, View):
    """
    View that helps download a file
    """

    @staticmethod
    def get(request, **kwargs):
        collection = Collection.objects.get(id=kwargs["collection_id"])
        filepath, filename = collection.generate_file(
            kwargs.get("current_view"),
            test=settings.Q_CLUSTER.get("sync", False),
            group_name=kwargs.get("group_name"),
            selected_ids=kwargs.get("selected_ids"),
        )

        wrapper = FileWrapper(open(filepath, "rb"))
        response = HttpResponse(wrapper, content_type="chemical/x-mdl-sdf")
        response["Content-Disposition"] = f'attachment; filename="{filename}"'

        # Delete file once it's attached to the response
        os.remove(filepath)

        return response


class ProjectView(LoginRequiredMixin, ListView):
    model = Compound
    context_object_name = "compounds"
    template_name = "main/project_home.html"

    def get_queryset(self):
        """
        Show all compounds for the relevant Project
        """
        compound_qs = super().get_queryset()
        try:
            self.project = Project.objects.get(code=self.kwargs["project_code"])
            # Retrieve a queryset of all compounds for this project, prefetching `compoundoccurrence_set
            # and the CompoundOccurrences' `collection_set`. These sets are prefetched to avoid
            # making a separate DB call for every compound in the collection (when the template loads)
            co_w_owner = (
                CompoundOccurrence.objects.select_related("owner")
                .order_by("pk")
                .prefetch_related("collection_set")
            )
            compound_qs = (
                compound_qs.filter(project=self.project, is_mmpdb_analogue=False)
                .select_related("series")
                .prefetch_related(
                    Prefetch("compoundoccurrence_set", queryset=co_w_owner)
                )
            )
            # If a structure search is specified, filter queryset
            search_type = self.kwargs.get("search_type")
            encrypted_smiles = self.kwargs.get("encrypted_smiles")
            if search_type and encrypted_smiles:
                smiles = signing.loads(encrypted_smiles)
                search_qs = []
                if search_type == "sss":
                    search_qs = Compound.objects.has_substruct(smiles)
                elif search_type == "exact":
                    search_qs = Compound.objects.is_equal(smiles)
                search_qs_pks = [c.pk for c in search_qs]
                compound_qs = compound_qs.filter(pk__in=search_qs_pks)

        except Project.DoesNotExist:
            return None

        # Generate SVG files for any compounds that are missing them
        for comp in compound_qs.filter(svg_file=""):
            comp.generate_svg_file()
        return compound_qs

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["project"] = self.project
        context["sss_form"] = SubstructureSearchForm()
        context["search_type"] = self.kwargs.get("search_type")
        if self.kwargs.get("encrypted_smiles"):
            context["smiles"] = signing.loads(self.kwargs.get("encrypted_smiles"))
            moltext = Chem.MolToMolBlock(Chem.MolFromSmiles(context["smiles"]))
            context["smiles_svg"] = moltext_to_svg(moltext, size_x=200, size_y=100)

        Analytic(
            PAGE_VIEW,
            self.kwargs["project_code"],
            self.request.user,
            page="project_view",
        )

        return context


class SubstructureSearchView(LoginRequiredMixin, FormView):
    form_class = SubstructureSearchForm
    http_method_names = ["post"]

    def form_valid(self, form, *args, **kwargs):
        project_code = self.kwargs["project_code"]
        search_type = form.cleaned_data["search_type"]
        smiles = form.cleaned_data["sketcher"]
        encrypted_smiles = signing.dumps(smiles)
        redirect_url = reverse(
            "project", args=[project_code, search_type, encrypted_smiles]
        )
        return JsonResponse({"errors": form.errors, "redirect_url": redirect_url})

    def form_invalid(self, form, *args, **kwargs):
        return JsonResponse({"errors": form.errors})


class HikeView(LoginRequiredMixin, UpdateView):
    model = Collection
    form_class = HikeForm

    slug_field = "id"
    slug_url_kwarg = "collection_id"
    http_method_names = ["post"]

    def form_valid(self, form, *args, **kwargs):
        analysis = form.cleaned_data["analysis"]
        if analysis == PROPCALC:
            redirect_url = reverse(
                "submit",
                kwargs={"nextview": analysis, "collection_id": form.instance.id},
            )
        else:
            redirect_url = reverse(analysis, kwargs={"collection_id": form.instance.id})

        self.request.user.add_easter_egg_point(self.request, "Hiking")
        return HttpResponseRedirect(redirect_url)


class SaveGemsView(LoginRequiredMixin, SaveGemsMixin, UpdateView):
    """
    View to change which gems are saved on a CompoundOccurrence model
    """

    form_class = SaveGemsForm
    model = CompoundOccurrence
    context_object_name = "parent_co"
    template_name = "main/homepage.html"

    slug_field = "pk"
    slug_url_kwarg = "compound_occurrence_id"
    http_method_names = ["post"]

    def get_form_kwargs(self, *args, **kwargs):
        analysis = self.kwargs["group_name"].split("_")[0]
        kwargs = super().get_form_kwargs(*args, **kwargs)
        if "confs" not in kwargs:
            task = (
                Task.objects.filter(
                    name=self.kwargs["task_name"], group=self.kwargs["group_name"]
                )
                .order_by("started")
                .last()
            )
            kwargs["confs"] = self.get_confs_dict_from_task(analysis, task)
        if "collection" not in kwargs:
            kwargs["collection"] = Collection.objects.get(
                pk=self.kwargs["collection_id"]
            )
        return kwargs

    def form_valid(self, form):
        """
        Processes a `SaveGemsForm`. Gems currently in the collection that are not in `gems`
        or `other_gems` are removed from the collection. New gems are added to the collection,
        creating a new CompoundOccurrence object if one does not yet exist.
        :param form: an instance of `SaveGemsForm`
        :returns: a JsonResponse to be parsed in the `postModalForm` function
        """
        collection = Collection.objects.get(pk=self.kwargs["collection_id"])
        analysis = self.kwargs["group_name"].split("_")[0]
        if form.is_valid():
            parent_co = self.object
            # Remove any CompoundOccurrences that were previously saved, but were unchecked in the form
            gems_to_keep = form.cleaned_data["other_gems"] + form.cleaned_data["gems"]
            unsaved_gem_ids = self.remove_unsaved_cos(
                gems_to_keep, collection, parent_co
            )
            # For all conformers in the form
            for conf_id in list(form.confs.keys()):
                mol = Chem.MolFromMolBlock(form.confs[conf_id]["moltext"])
                conf_molblock = Chem.MolToMolBlock(mol)
                # Get the existing CO object if it exists
                co = self.get_existing_co(parent_co, conf_id, conf_molblock)
                # If the conformer was selected, add it to the Collection
                if conf_id in form.cleaned_data["gems"]:
                    if not co:
                        co = CompoundOccurrence.objects.create(
                            parent_co=parent_co,
                            compound=parent_co.compound,
                            owner=parent_co.owner,
                            molblock=conf_molblock,
                            gem_id=conf_id,
                            saved_from=analysis,
                        )
                    # Add new compound occurrence to collection (if not there already)
                    collection.compound_occurrences.add(co)
                # If the gem was not selected, remove it from the Collection
                else:
                    if co:
                        collection.compound_occurrences.remove(co)
        return JsonResponse(
            {
                "errors": form.errors,
                "close_modal": True,
                "success_method": f"updateSavedGems({self.object.compound.pk}, {self.object.pk}, {form.cleaned_data['gems']}, {unsaved_gem_ids}, '{self.get_select_class(analysis)}')",
            }
        )

    def form_invalid(self, form):
        form.is_valid()
        return JsonResponse({"errors": form.errors})


class AddCompoundView(LoginRequiredMixin, UpdateView):
    model = Collection
    form_class = AddCompoundForm

    slug_field = "id"
    slug_url_kwarg = "collection_id"
    http_method_names = ["post"]

    def form_valid(self, form, *args, **kwargs):
        self.analysis = self.kwargs["current_view"]
        self.group_name = self.kwargs.get("group_name")
        self.collection = self.object
        cos_added = form.cleaned_data["cos_added"]
        table_rows = [self.get_table_row(co) for co in cos_added]
        grid_items = []
        if self.analysis == PROPCALC:
            grid_items = [self.get_grid_item(co) for co in cos_added]
        if self.group_name:
            for co in cos_added:
                self.start_analysis(co)

        return JsonResponse(
            {
                "errors": form.errors,
                "close_modal": True,
                "success_method": f"addComps('{self.analysis}', {table_rows}, {grid_items})",
            }
        )

    def form_invalid(self, form):
        return JsonResponse({"errors": form.errors})

    def get_table_row(self, co):
        """
        Generate a <tr> HTML string to add to the table w/ the newly added CompoundOccurrence
        :param co: a CompoundOccurrence object to add to the table
        :returns: a string of HTML
        """
        context = {
            "collection": self.collection,
            "co": co,
            "group_name": self.group_name,
            "new_row": "true",
        }
        if self.analysis == PROPCALC:
            context["props"] = self.collection.get_propcalc_column_headers()
        elif self.analysis in [ALIGN, DOCK, ESP]:
            context["ref_string"] = (
                self.group_name.split("_")[-1] if self.group_name else ""
            )
        elif self.analysis == TORSION and self.group_name:
            pioneer_pk = int(self.group_name.split("_")[2])
            self.pioneer = CompoundOccurrence.objects.get(pk=pioneer_pk)
            self.pioneer_dihedral_atoms = self.group_name.split("_")[3].replace(
                "-", ","
            )
            self.pioneer_dihedral_smarts = self.pioneer.convert_atoms_to_smarts(
                self.pioneer_dihedral_atoms
            )
            context["pioneer_dihedral_atoms"] = self.pioneer_dihedral_atoms
            context["pioneer_dihedral_smarts"] = self.pioneer_dihedral_smarts

        table_row = render(
            self.request, f"main/{self.analysis}/table_row.html", context
        ).content.decode()
        return table_row

    def get_grid_item(self, co):
        """
        Generate a grid-item HTML string to add to the propcalc grid w/ the newly added CompoundOccurrence
        :param co: a CompoundOccurrence object to add to the table
        :returns: a string of HTML
        """
        context = {"co": co, "property_results": {}}
        grid_item = render(
            self.request, f"main/{self.analysis}/grid_item.html", context
        ).content.decode()
        return grid_item

    def start_analysis(self, co):
        """
        Starts a DjangoQ task for the newly added CompoundOccurrence
        :param co: a CompoundOccurrence object
        """
        if self.analysis == PROPCALC:
            async_task(
                co.compound.run_propcalc,
                inductive=self.collection.inductive_in_props(),
                task_name=co.get_propcalc_task_name(self.collection),
                group=self.collection.get_propcalc_group_name(),
                cluster="fast",
            )
        elif self.analysis == ALIGN:
            ref_string = self.group_name.split("_")[-1]
            async_task(
                co.superimpose_to_ref,
                reference=self.collection._get_ref_obj(ref_string),
                task_name=co.get_align_task_name(ref_string),
                group=self.collection.get_align_group_name(ref_string),
            )
        elif self.analysis == DOCK:
            ref_string = self.group_name.split("_")[-1]
            async_task(
                co.dock_to_receptor,
                reference=self.collection._get_ref_obj(ref_string),
                task_name=co.get_dock_task_name(ref_string),
                group=self.collection.get_dock_group_name(ref_string),
            )
        elif self.analysis == ESP:
            ref_string = self.group_name.split("_")[-1]
            task_name = co.get_esp_task_name(ref_string)
            async_task(
                co.generate_esp_map,
                reference=self.collection._get_ref_obj(ref_string),
                job_name=task_name,
                task_name=task_name,
                group=self.collection.get_esp_group_name(ref_string),
            )
        elif self.analysis == TORSION:
            _, dihedral_atoms = co.pick_most_relevant_dihedral(
                self.pioneer_dihedral_smarts, self.pioneer_dihedral_atoms
            )
            task_name = co.get_torsion_task_name(
                self.pioneer_dihedral_smarts, dihedral_atoms
            )
            async_task(
                co.generate_torsions,
                dihedral_atoms=dihedral_atoms,
                dq_task_name=task_name,
                group=self.collection.get_torsion_group_name(
                    self.pioneer.pk, self.pioneer_dihedral_atoms
                ),
                task_name=task_name,
                cluster="slow",
            )


class KetcherView(TemplateView):
    template_name = "main/components/widgets/ketcher_template.html"

    @xframe_options_sameorigin
    def get(self, request, *args, **kwargs):
        return super().get(request, *args, **kwargs)
