from django.conf import settings
from django.conf.urls import include
from django.conf.urls.static import static
from django.contrib import admin
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import LogoutView
from django.urls import path
from django.views.generic import TemplateView

from basechem.main.views.base_views import KetcherView

urlpatterns = [
    path(
        r"home/",
        login_required(TemplateView.as_view(template_name="base.html")),
        name="home",
    ),
    path(r"ketcher/", KetcherView.as_view(), name="ketcher"),
    path(r"admin/", admin.site.urls),
    path(r"logout/", LogoutView.as_view(), name="logout"),
    path(r"", include("basechem.sso.urls")),
    path(r"", include("basechem.main.urls")),
    path(r"", include("basechem.mni_common.urls")),
]

if settings.DEBUG:
    import debug_toolbar

    urlpatterns += [
        path(r"__debug__/", include(debug_toolbar.urls)),
        path("sentry-debug/", lambda _: 1 / 0),
    ]
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
