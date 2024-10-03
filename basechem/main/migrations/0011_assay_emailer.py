# Generated by Django 3.2.7 on 2022-07-14 08:55

import django.core.serializers.json
from django.conf import settings
from django.db import migrations, models

import basechem.main.models.project_models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ("main", "0010_collection_name"),
    ]

    operations = [
        migrations.AddField(
            model_name="project",
            name="assays",
            field=models.JSONField(
                blank=True,
                decoder=basechem.main.models.project_models.AssaysDecoder,
                default=basechem.main.models.project_models.assays_default_value,
                encoder=django.core.serializers.json.DjangoJSONEncoder,
                null=True,
            ),
        ),
        migrations.AddField(
            model_name="project",
            name="leads",
            field=models.ManyToManyField(
                related_name="leads", to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name="project",
            name="subscribers",
            field=models.ManyToManyField(
                related_name="subscribers", to=settings.AUTH_USER_MODEL
            ),
        ),
        migrations.AddField(
            model_name="project",
            name="target",
            field=models.CharField(default="", max_length=20),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name="series",
            name="assay_data",
            field=models.JSONField(blank=True, default=dict, null=True),
        ),
    ]
