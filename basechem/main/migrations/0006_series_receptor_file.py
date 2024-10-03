# Generated by Django 3.2.7 on 2022-02-25 19:27

from django.db import migrations, models

import basechem.main.models.compound_models
import basechem.mni_common.storage


class Migration(migrations.Migration):

    dependencies = [
        ("main", "0005_auto_20211206_1225"),
    ]

    operations = [
        migrations.AddField(
            model_name="series",
            name="receptor_file",
            field=models.FileField(
                blank=True,
                null=True,
                storage=basechem.mni_common.storage.select_media_storage,
                upload_to=basechem.main.models.compound_models.series_files_path,
            ),
        ),
    ]
