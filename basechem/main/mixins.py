from django.core.mail import mail_admins
from django.views.generic.base import ContextMixin

from basechem.common.constants import ADMIN_FAILURE, LOADING_MESSAGE
from basechem.main.constants import ALIGN, DOCK, TORSION
from basechem.main.forms import AddCompoundForm, HikeForm
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import CompoundOccurrence


class AnalysisViewMixin(ContextMixin):
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Get current view and collection from url
        [current_view, collection_id] = self.request.path.split("/")[1:3]
        collection = Collection.objects.get(id=collection_id)
        context["current_view"] = current_view
        context["collection"] = collection
        context["hike_form"] = HikeForm(instance=collection, current_view=current_view)
        context["add_comp_form"] = AddCompoundForm(analysis=current_view)
        context["compound_occurrences"] = collection.get_cos_for_analysis(current_view)
        context["saved_gem_ids"] = list(
            collection.compound_occurrences.all().values_list("gem_id", flat=True)
        )
        context["loading_message"] = LOADING_MESSAGE
        return context


class SaveGemsMixin(ContextMixin):
    def remove_unsaved_cos(self, gems_to_keep, collection, parent_co):
        """
        Remove CompoundOccurrence objects from the collection if they were unsaved
        :param gems_to_keep: a list of gem_ids that should not be removed
        :param collection: the collection being edited
        :param parent_co: the parent CompoundOccurrence object
        :param return: a list of gem_ids that were unsaved from this collection
        """
        cos_to_remove = CompoundOccurrence.objects.filter(
            collection=collection, owner=parent_co.owner, parent_co=parent_co
        ).exclude(gem_id__in=gems_to_keep)
        collection.compound_occurrences.remove(*cos_to_remove)
        return [co.gem_id for co in cos_to_remove]

    def get_existing_co(self, parent_co, conf_id, conf_molblock):
        """
        Given a conformer ID and molblock, check if there are any existing CompoundOccurrences
        that match. If so, return the matching CompoundOccurrence.
        :param parent_co: the parent CompoundOccurrence object
        :param conf_id: a string, the conf_id of the gem being processed
        :param conf_molblock: a string, the gem's molblock
        :returns: an existing CompoundOccurrence object for this gem if one exists
        """
        matching_cos = CompoundOccurrence.gems.filter(
            parent_co=parent_co, gem_id=conf_id, molblock=conf_molblock
        )
        co = None
        if matching_cos:
            co = matching_cos.first()
        if matching_cos.count() > 1:
            admin_message = f"Multiple CompoundOccurrences found with the following attributes:\n- parent_co: {parent_co.pk}\n- gem_id: {conf_id}\n- conf_molblock:\n{conf_molblock}"
            mail_admins(ADMIN_FAILURE, admin_message)
        return co

    def get_confs_dict_from_group(self, analysis, task_results, co_pk):
        """
        Given the task results for an analysis, return a dictionary of conformers for the given `co_pk`
        :param analysis: a string, a basechem analysis (ALIGN, DOCK, etc.)
        :param task_results: a dictionary of results for the task group
        :param co_pk: an integer, the pk of the compound occurrence whose conformers are being retrieved
        :returns: a dictionary of conformers where the keys are conf_ids
        """
        if not task_results:
            return {}
        elif analysis in [ALIGN, DOCK]:
            return task_results["compounds"][f"co-{co_pk}"]
        elif analysis == TORSION:
            return task_results[f"co-{co_pk}"]["torsions"]
        return {}

    def get_confs_dict_from_task(self, analysis, task):
        """
        Given a DjangoQ task for an analysis, return a dictionary of conformers
        :param analysis: a string, a basechem analysis (ALIGN, DOCK, etc.)
        :param task: a DjangoQ Task object
        :returns: a dictionary of conformers where the keys are conf_ids
        """
        if not task:
            return {}
        elif analysis in [ALIGN, DOCK]:
            return task.result
        elif analysis == TORSION:
            return task.result["torsions"]
        return {}

    def get_select_class(self, analysis):
        """
        Get html class that should be used as a parameter to `updateSavedGems` for this analysis
        :param analysis: a string, a basechem analysis (ALIGN, DOCK, etc.)
        :returns: a string, the html class of the select elements where users can choose
        conformers to display in the viewer
        """
        if analysis == ALIGN:
            return "select-conformer"
        elif analysis == DOCK:
            return "select-pose"
        return ""
