import json
import re
from collections import defaultdict

from django.contrib.auth.mixins import LoginRequiredMixin
from django.http.response import JsonResponse
from django.shortcuts import render
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import View
from django_q.models import OrmQ, Task
from django_q.tasks import result
from rdkit import Chem

from basechem.main.constants import *
from basechem.main.forms import SaveGemsForm
from basechem.main.mixins import SaveGemsMixin
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import CompoundOccurrence
from basechem.main.tasks import get_group_results


class CheckTaskView(LoginRequiredMixin, View):
    """
    View to get the result of an asynchronous task or task group
    """

    def get_group_response(self, group_name):
        _, _, group_result = get_group_results(group_name)
        if group_result is not None:
            return JsonResponse({"result": group_result})
        return JsonResponse({})

    def get_task_response(self, task_id):
        task_result = result(task_id)
        if task_result != None:
            return JsonResponse({"result": task_result})
        return JsonResponse({})

    def get(self, request, *args, **kwargs):
        if request.GET.get("task_type") == "group":
            group_name = request.GET.get("identifier")
            return self.get_group_response(group_name)

        elif request.GET.get("task_type") == "task":
            task_id = request.GET.get("identifier")
            return self.get_task_response(task_id)


@method_decorator(csrf_exempt, name="dispatch")
class SaveViewerGemsView(LoginRequiredMixin, SaveGemsMixin, View):
    """
    This view handles ajax POST requests when a user clicks the "save all conformers in the viewer" button
    """

    def post(self, *args, **kwargs):
        request_dict = json.loads(self.request.body)
        _, _, task_results = get_group_results(request_dict["groupName"])
        analysis = request_dict["groupName"].split("_")[0]
        collection = Collection.objects.get(pk=kwargs["collection_id"])
        parent_cos = collection.compound_occurrences.filter(parent_co=None)
        errors = []
        # Construct dictionaries of previously saved and to-save conformers
        prev_saved = defaultdict(list)
        for co in collection.compound_occurrences.exclude(
            parent_co=None
        ).select_related("parent_co", "compound"):
            prev_saved[co.parent_co.pk].append(co.gem_id)
        to_save = defaultdict(list)
        for conf_id in request_dict["viewerConfIds"]:
            parent_co_pk = int(re.match("c\d+-co(\d+)", conf_id)[1])
            to_save[parent_co_pk].append(conf_id)

        cos_to_save = []
        for parent_co in parent_cos:
            confs_to_save = to_save[parent_co.pk]
            already_saved = prev_saved[parent_co.pk]
            total_gems = len(set(confs_to_save + already_saved))
            if total_gems > MAX_GEMS_PER_CO:
                # Confirm that at most 3 gems will be saved to each CO
                errors.append(
                    f"Attempting to save {total_gems} conformers on {parent_co.compound.name}."
                )
            else:
                # Save new conformer IDs
                confs = self.get_confs_dict_from_group(
                    analysis, task_results, parent_co.pk
                )
                for conf_id in confs_to_save:
                    mol = Chem.MolFromMolBlock(confs[conf_id]["moltext"])
                    conf_molblock = Chem.MolToMolBlock(mol)
                    # Get the existing CO object if it exists
                    co = self.get_existing_co(parent_co, conf_id, conf_molblock)
                    if not co:
                        co = CompoundOccurrence.objects.create(
                            parent_co=parent_co,
                            compound=parent_co.compound,
                            owner=parent_co.owner,
                            molblock=conf_molblock,
                            gem_id=conf_id,
                            saved_from=analysis,
                        )
                    cos_to_save.append(co)
        # Add all new compound_occurrences at once to reduce db hits
        collection.compound_occurrences.add(*cos_to_save)
        success_methods = []
        for parent_co in parent_cos:
            confs_to_save = to_save[parent_co.pk] + prev_saved[parent_co.pk]
            if confs_to_save:
                success_methods.append(
                    f"updateSavedGems({parent_co.compound.pk}, {parent_co.pk}, {confs_to_save}, [], '{self.get_select_class(analysis)}')"
                )

        if errors:
            errors.append(f"Maximum is {MAX_GEMS_PER_CO}.")
        return JsonResponse({"errors": errors, "success_methods": success_methods})


@method_decorator(csrf_exempt, name="dispatch")
class UpdateCollectionView(LoginRequiredMixin, View):
    def post(self, *args, **kwargs):
        request_dict = json.loads(self.request.body)
        collection = Collection.objects.get(pk=kwargs["collection_id"])
        update_field = request_dict.get("update_field")
        if update_field == "co_order":
            collection.update_co_order(request_dict["co_order"])
        return JsonResponse({})


class AjaxCollectTaskView(LoginRequiredMixin, View):
    def get(self, *args, **kwargs):
        """
        Expects a request with `collection_pk`, `group_name`, and `task_name`. Returns a JsonResponse
        that includes task status and results.
        :returns: a JsonResponse w/ variable keys depending on the given task_name. The response
            always includes a key "taskStatus".
        """
        self.collection = Collection.objects.get(
            pk=self.request.GET.get("collection_pk")
        )
        self.group_name = self.request.GET.get("group_name")
        self.task_name = self.request.GET.get("task_name")
        self.task = (
            Task.objects.filter(name=self.task_name, group=self.group_name)
            .order_by("started")
            .last()
        )
        # Stop here if the task doesn't have results
        self.update_task_status()
        if self.task_status != COMPLETE:
            return JsonResponse({"taskStatus": self.task_status})
        if PROPCALC in self.group_name:
            return self.collect_propcalc_task()
        elif ALIGN in self.group_name:
            if "alignrefs" in self.task_name:
                return self.collect_alignrefs_task()
            else:
                return self.collect_align_task()
        elif DOCK in self.group_name:
            if "dockrefs" in self.task_name:
                return self.collect_dockrefs_task()
            else:
                return self.collect_dock_task()
        elif ESP in self.group_name:
            if "esprefs" in self.task_name:
                return self.collect_esprefs_task()
            else:
                return self.collect_esp_task()
        elif TORSION in self.group_name:
            return self.collect_torsion_task()
        elif MMP in self.group_name:
            return self.collect_mmp_task()

    def update_task_status(self):
        """
        Updates `self.task_status` based on the status of `self.task_name`
        """
        self.task_status = COMPLETE
        if not self.task:
            queued = [
                oq for oq in OrmQ.objects.all() if oq.task.get("name") == self.task_name
            ]
            if queued:
                self.task_status = IN_PROGRESS
            else:
                self.task_status = DROPPED
        elif not self.task.success:
            self.task_status = ERROR

    ######################################
    #              PROPCALC              #
    ######################################

    def collect_propcalc_task(self):
        """
        Called when `self.task_name` is a propcalc task, this function collects the results of the
        task and generates HTML for the row in table-view and grid item in grid-view so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of conformers (if the task completed successfully),
            tableRow: a string of HTML containing a <tr> element with data to populate the table view,
            gridItem: a string of HTML containing a <div> element with data to populate the grid view,
        }
        """
        self.co = CompoundOccurrence.objects.get(pk=self.request.GET.get("co_pk"))
        context = {
            "co": self.co,
            "collection": self.collection.pk,
            "group_name": self.group_name,
            "property_results": {self.co.pk: self.task.result},
            "props": self.collection.get_propcalc_column_headers(),
        }
        table_row = render(
            self.request, f"main/propcalc/table_row.html", context
        ).content.decode()
        grid_item = render(
            self.request, f"main/propcalc/grid_item.html", context
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "tableRow": table_row,
                "gridItem": grid_item,
            }
        )

    ###################################
    #              ALIGN              #
    ###################################

    def collect_alignrefs_task(self):
        """
        Called when `self.task_name` is an alignrefs task, this function collects the results of the
        task and generates HTML for the reference dropdown so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the alignrefs task,
            taskResult: a dictionary of reference and receptor data (if the task completed successfully),
            refOptions: a string of HTML with <option> tags to populate the display-reference dropdown
        }
        """
        ref_string = self.task_name.split("_")[-1]
        initial_series = self.collection.most_relevant_series(
            ref_string, list(self.task.result["references"].keys())
        )
        reference_options = render(
            self.request,
            "main/components/ajax/reference_options.html",
            {
                "references": self.task.result["references"],
                "initial_series": initial_series,
            },
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "refOptions": reference_options,
            }
        )

    def collect_align_task(self):
        """
        Called when `self.task_name` is an align task, this function collects the results of the
        task and generates HTML for the conformer dropdown and save_gems_modal so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of conformers (if the task completed successfully),
            saveGemsModal: a string of HTML with the save-gems modal for this CompoundOccurrence,
            confOptions: a string of HTML with <option> tags to populate the select-conformers dropdown
        }
        """
        self.co = CompoundOccurrence.objects.get(pk=self.request.GET.get("co_pk"))
        form = SaveGemsForm(
            instance=self.co,
            collection=self.collection,
            confs=self.task.result,
            keys=["r_bc_rmsd_to_lsalign", "r_mmff_rel_energy"],
        )
        save_gems_modal = render(
            self.request,
            "main/components/modals/save_gems.html",
            {
                "collection": self.collection,
                "co": self.co,
                "group_name": self.group_name,
                "task_name": self.task_name,
                "form": form,
            },
        ).content.decode()

        select_conf_options = render(
            self.request,
            "main/align/conf_options.html",
            {
                "confs": self.task.result,
                "co": self.co,
                "saved_gem_ids": list(
                    self.collection.compound_occurrences.filter(
                        parent_co=self.co
                    ).values_list("gem_id", flat=True)
                ),
            },
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "saveGemsModal": save_gems_modal,
                "confOptions": select_conf_options,
            }
        )

    ###################################
    #              DOCK               #
    ###################################

    def collect_dockrefs_task(self):
        """
        Called when `self.task_name` is a dockrefs task, this function collects the results of the
        task and generates HTML for the display reference and display receptor dropdowns so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of reference and receptor data (if the task completed successfully),
            refOptions: a string of HTML with <option> tags to populate the display-reference dropdown,
            recOptions: a string of HTML with <option> tags to populate the display-receptor dropdown
        }
        """
        ref_string = self.task_name.split("_")[-1]
        initial_series = self.collection.most_relevant_series(
            ref_string, list(self.task.result["references"].keys())
        )
        reference_options = render(
            self.request,
            "main/components/ajax/reference_options.html",
            {
                "references": self.task.result["references"],
                "initial_series": initial_series,
            },
        ).content.decode()
        receptor_options = render(
            self.request,
            "main/components/ajax/reference_options.html",
            {
                "references": self.task.result["references"],
                "initial_series": initial_series,
            },
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "refOptions": reference_options,
                "recOptions": receptor_options,
            }
        )

    def collect_dock_task(self):
        """
        Called when `self.task_name` is a dock task, this function collects the results of the
        task and generates HTML for the pose dropdown and save_gems_modal so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of poses(if the task completed successfully),
            saveGemsModal: a string of HTML with the save-gems modal for this CompoundOccurrence,
            confOptions: a string of HTML with <option> tags to populate the select-poses dropdown
        }
        """
        self.co = CompoundOccurrence.objects.get(pk=self.request.GET.get("co_pk"))
        form = SaveGemsForm(
            instance=self.co,
            collection=self.collection,
            confs=self.task.result,
            keys=["toklatScore", "rdockScore", "RMSDtoLSAligned"],
        )
        save_gems_modal = render(
            self.request,
            "main/components/modals/save_gems.html",
            {
                "collection": self.collection,
                "co": self.co,
                "group_name": self.group_name,
                "task_name": self.task_name,
                "form": form,
            },
        ).content.decode()
        select_pose_options = render(
            self.request,
            "main/dock/pose_options.html",
            {
                "poses": self.task.result,
                "co": self.co,
                "saved_gem_ids": list(
                    self.collection.compound_occurrences.filter(
                        parent_co=self.co
                    ).values_list("gem_id", flat=True)
                ),
            },
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "saveGemsModal": save_gems_modal,
                "confOptions": select_pose_options,
            }
        )

    ###################################
    #               ESP               #
    ###################################

    def collect_esprefs_task(self):
        """
        Called when `self.task_name` is an esprefs task, this function collects the results of the
        task and generates HTML for the reference dropdown so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of reference and receptor data (if the task completed successfully),
            refOptions: a string of HTML with <option> tags to populate the display-reference dropdown
        }
        """
        ref_string = self.task_name.split("_")[-1]
        initial_series = self.collection.most_relevant_series(
            ref_string, list(self.task.result["references"].keys())
        )
        reference_options = render(
            self.request,
            "main/components/ajax/reference_options.html",
            {
                "references": self.task.result["references"],
                "initial_series": initial_series,
            },
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "refOptions": reference_options,
            }
        )

    def collect_esp_task(self):
        """
        Called when `self.task_name` is an esp task, this function collects the results of the
        task so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of ESP map data (if the task completed successfully),
        }
        """
        self.co = CompoundOccurrence.objects.get(pk=self.request.GET.get("co_pk"))
        return JsonResponse(
            {"taskStatus": self.task_status, "taskResult": self.task.result}
        )

    #####################################
    #              TORSION              #
    #####################################

    def collect_torsion_task(self):
        """
        Called when `self.task_name` is a torsion task, this function collects the results of the
        task and generates HTML for the save_gems_modal so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the  task,
            taskResult: a dictionary of torsion data (if the task completed successfully),
            saveGemsModal: a string of HTML with the save-gems modal for this CompoundOccurrence,
        }
        """
        self.co = CompoundOccurrence.objects.get(pk=self.request.GET.get("co_pk"))
        form = SaveGemsForm(
            instance=self.co,
            collection=self.collection,
            confs=self.task.result["torsions"],
            keys=["rel_energy"],
        )
        save_gems_modal = render(
            self.request,
            "main/components/modals/save_gems.html",
            {
                "collection": self.collection,
                "co": self.co,
                "group_name": self.group_name,
                "task_name": self.task_name,
                "form": form,
            },
        ).content.decode()
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "saveGemsModal": save_gems_modal,
            }
        )

    ######################################
    #              MMP              #
    ######################################

    def collect_mmp_task(self):
        """
        Called when `self.task_name` is an MMP task, this function collects the results of the
        task and generates HTML for the grid item in grid-view so the DOM can be updated appropriately.
        :returns: a JsonResponse of the form: {
            taskStatus: a string, the status of the task,
            taskResult: a dictionary of conformers (if the task completed successfully),
            dnGrid: a string of HTML containing a <div> element with a grid view of DN MMPs,
            ideaGrid: a string of HTML containing a <div> element with a grid view of non-DN MMPS,
        }
        """
        self.co = CompoundOccurrence.objects.get(pk=self.request.GET.get("co_pk"))
        # Get queryset of MMPs filtered by constant region
        mmps = self.co.compound.mmps.filter(pk__in=self.task.result)
        dn_mmps = mmps.exclude(dn_id="")
        vsm_mmps = mmps.filter(dn_id="", is_mmpdb_analogue=False)

        dn_grid = render(self.request, f"main/mmps/mmps_grid.html", {"mmps": dn_mmps})
        idea_grid = render(
            self.request, f"main/mmps/mmps_grid.html", {"mmps": vsm_mmps}
        )
        return JsonResponse(
            {
                "taskStatus": self.task_status,
                "taskResult": self.task.result,
                "dnGrid": dn_grid.content.decode(),
                "ideaGrid": idea_grid.content.decode(),
            }
        )
