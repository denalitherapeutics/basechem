# Generated by Django 3.2.7 on 2021-12-06 17:25

from django.db import migrations, models

import basechem.main.models.compound_models
import basechem.mni_common.storage


class Migration(migrations.Migration):

    dependencies = [
        ("main", "0004_auto_20211103_1907"),
    ]

    operations = [
        migrations.AddField(
            model_name="compound",
            name="confgenx_file",
            field=models.FileField(
                blank=True,
                null=True,
                storage=basechem.mni_common.storage.select_media_storage,
                upload_to=basechem.main.models.compound_models.compound_files_path,
            ),
        ),
        migrations.AlterField(
            model_name="series",
            name="bound_state_file",
            field=models.FileField(
                blank=True,
                null=True,
                storage=basechem.mni_common.storage.select_media_storage,
                upload_to=basechem.main.models.compound_models.series_files_path,
            ),
        ),
        migrations.AlterField(
            model_name="series",
            name="ground_state_file",
            field=models.FileField(
                blank=True,
                null=True,
                storage=basechem.mni_common.storage.select_media_storage,
                upload_to=basechem.main.models.compound_models.series_files_path,
            ),
        ),
    ]
