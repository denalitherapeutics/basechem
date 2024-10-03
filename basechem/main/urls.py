from django.urls import path
from django.views.generic import TemplateView

from basechem.main.constants import ALIGN, DOCK, DTX_MMP, ESP, MMP, PROPCALC, TORSION
from basechem.main.views.ajax_views import (
    AjaxCollectTaskView,
    CheckTaskView,
    SaveViewerGemsView,
    UpdateCollectionView,
)
from basechem.main.views.analysis_views import (
    DockView,
    DTXMMPView,
    EspView,
    LigandAlignView,
    MMPSearchView,
    PropCalcView,
    TorsionView,
)
from basechem.main.views.base_views import (
    AddCompoundView,
    DownloadResultFileView,
    HikeView,
    HomePageView,
    ProjectView,
    SaveGemsView,
    SubmitCompoundsView,
    SubstructureSearchView,
)

urlpatterns = [
    path(r"homepage/", HomePageView.as_view(), name="homepage"),
    path(r"faq/", TemplateView.as_view(template_name="main/faq.html"), name="faq"),
    path(r"project/<str:project_code>/", ProjectView.as_view(), name="project"),
    path(
        r"project/<str:project_code>/<str:search_type>/<str:encrypted_smiles>",
        ProjectView.as_view(),
        name="project",
    ),
    path(r"sss/<str:project_code>/", SubstructureSearchView.as_view(), name="sss"),
    path(r"<str:nextview>/submit/", SubmitCompoundsView.as_view(), name="submit"),
    path(
        r"<str:nextview>/submit/<int:collection_id>",
        SubmitCompoundsView.as_view(),
        name="submit",
    ),
    path(
        r"<str:nextview>/submit/<int:collection_id>/<str:task_id>",
        SubmitCompoundsView.as_view(),
        name="submit",
    ),
    path(r"propcalc/<int:collection_id>", PropCalcView.as_view(), name=PROPCALC),
    path(r"align/<int:collection_id>", LigandAlignView.as_view(), name=ALIGN),
    path(
        r"align/<int:collection_id>/<str:group_name>",
        LigandAlignView.as_view(),
        name=ALIGN,
    ),
    path(r"dock/<int:collection_id>", DockView.as_view(), name=DOCK),
    path(r"dock/<int:collection_id>/<str:group_name>", DockView.as_view(), name=DOCK),
    path(r"esp/<int:collection_id>", EspView.as_view(), name=ESP),
    path(r"esp/<int:collection_id>/<str:group_name>", EspView.as_view(), name=ESP),
    path(r"torsion/<int:collection_id>", TorsionView.as_view(), name=TORSION),
    path(
        r"torsion/<int:collection_id>/<str:group_name>",
        TorsionView.as_view(),
        name=TORSION,
    ),
    path(r"dtx-mmp/<int:collection_id>", DTXMMPView.as_view(), name=DTX_MMP),
    path(r"mmp/<int:collection_id>", MMPSearchView.as_view(), name=MMP),
    path(
        r"mmp/<int:collection_id>/<str:group_name>", MMPSearchView.as_view(), name=MMP
    ),
    ############################### Modal/Ajax URLs ###############################
    path(r"ajax/check-group/", CheckTaskView.as_view(), name="ajax_check_task_group"),
    path(
        r"ajax/update-collection/<int:collection_id>/",
        UpdateCollectionView.as_view(),
        name="ajax_update_collection",
    ),
    path(
        r"ajax/save-viewer-gems/<int:collection_id>/",
        SaveViewerGemsView.as_view(),
        name="ajax_save_viewer_gems",
    ),
    path(r"hike/<int:collection_id>", HikeView.as_view(), name="hike"),
    path(
        r"save-gems/<int:collection_id>/<int:compound_occurrence_id>/<str:group_name>/<str:task_name>",
        SaveGemsView.as_view(),
        name="save_gems",
    ),
    path(
        r"<str:current_view>/add-comp/<int:collection_id>",
        AddCompoundView.as_view(),
        name="add_comp",
    ),
    path(
        r"<str:current_view>/add-comp/<int:collection_id>/<str:group_name>",
        AddCompoundView.as_view(),
        name="add_comp",
    ),
    path(
        r"ajax/collect-task/", AjaxCollectTaskView.as_view(), name="ajax_collect_task"
    ),
    ############################### Download URLs ###############################
    path(
        r"<str:current_view>/<int:collection_id>/download/<str:group_name>/<str:selected_ids>",
        DownloadResultFileView.as_view(),
        name="downloadresultfile",
    ),
    path(
        r"<str:current_view>/<int:collection_id>/download/<str:group_name>/",
        DownloadResultFileView.as_view(),
        name="downloadresultfile",
    ),
    path(
        r"<str:current_view>/<int:collection_id>/download",
        DownloadResultFileView.as_view(),
        name="downloadresultfile",
    ),
]
