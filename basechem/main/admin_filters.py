from django.contrib import admin

from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, CompoundOccurrence


class FreeTextFilter(admin.SimpleListFilter):
    """
    Adds a filter to the django admin where the user can search using a free text input.
    Filters should inherit this class, stipulate a `parameter_name`, `title`, and an
    override for `queryset`. Add filters to an admin page using the `list_filter` property
    in the `ModelAdmin` object.
    """

    template = "main/admin/text_filter.html"

    def lookups(self, request, model_admin):
        # Dummy, required to show the filter.
        return ((),)

    def choices(self, changelist):
        choice = next(super().choices(changelist))
        yield choice


class ProjectFilter(FreeTextFilter):
    parameter_name = "project"
    title = "Project (code)"

    def queryset(self, request, queryset):
        code = self.value()
        if not code:
            return
        try:
            if queryset.model is Compound:
                return queryset.filter(project__code=code)
            elif queryset.model is CompoundOccurrence:
                return queryset.filter(compound__project__code=code)
            elif queryset.model is Collection:
                return queryset.filter(project__code=code)
        except Exception:
            return queryset.none()


class CompoundOccurrenceFilter(FreeTextFilter):
    parameter_name = "compound_occurrence"
    title = "Compound occurrence (PK)"

    def queryset(self, request, queryset):
        co_id = self.value()
        if not co_id:
            return
        try:
            if queryset.model is Compound:
                comp = CompoundOccurrence.objects.get(id=co_id).compound
                return queryset.filter(id=comp.id)
            elif queryset.model is CompoundOccurrence:
                return queryset.filter(id=co_id)
        except Exception:
            return queryset.none()


class CollectionFilter(FreeTextFilter):
    parameter_name = "collection"
    title = "Collection (PK)"

    def queryset(self, request, queryset):
        coll_id = self.value()
        if not coll_id:
            return
        try:
            if queryset.model is Compound:
                return Collection.objects.get(id=coll_id).compounds()
            elif queryset.model is CompoundOccurrence:
                coll = Collection.objects.get(id=coll_id)
                queryset &= coll.compound_occurrences.all()
                return queryset
            elif queryset.model is Collection:
                return queryset.filter(id=coll_id)
        except Exception:
            return queryset.none()


class CompoundPkFilter(FreeTextFilter):
    parameter_name = "compound"
    title = "Compound (PK)"

    def queryset(self, request, queryset):
        comp_id = self.value()
        if not comp_id:
            return
        try:
            if queryset.model is Compound:
                return queryset.filter(id=comp_id)
            elif queryset.model is CompoundOccurrence:
                return queryset.filter(compound__id=comp_id)
        except Exception:
            return queryset.none()


class CompoundDnIdFilter(FreeTextFilter):
    parameter_name = "compound_dn"
    title = "Compound (DN ID)"

    def queryset(self, request, queryset):
        dn_id = self.value()
        if not dn_id:
            return
        try:
            if queryset.model is Compound:
                return queryset.filter(dn_id__contains=dn_id)
            elif queryset.model is CompoundOccurrence:
                return queryset.filter(compound__dn_id__contains=dn_id)
        except Exception:
            return queryset.none()


class UserFilter(FreeTextFilter):
    parameter_name = "owner"
    title = "Owner (last name)"

    def queryset(self, request, queryset):
        last_name = self.value()
        if not last_name:
            return
        try:
            if queryset.model is Compound:
                cos = CompoundOccurrence.objects.filter(
                    owner__last_name__iexact=last_name
                )
                comp_ids = [co.compound_id for co in cos]
                return queryset.filter(id__in=comp_ids)
            elif queryset.model is CompoundOccurrence:
                return queryset.filter(owner__last_name__iexact=last_name)
            elif queryset.model is Collection:
                return queryset.filter(owner__last_name__iexact=last_name)
        except Exception:
            return queryset.none()
