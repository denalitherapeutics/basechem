# Generated by Django 3.2.7 on 2023-02-13 23:01

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("main", "0021_remove_nulls_from_char_text_file"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="compound",
            name="prop_alogd",
        ),
    ]
