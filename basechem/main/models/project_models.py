import datetime
import json

from django.core.serializers.json import DjangoJSONEncoder
from django.db import models
from django.db.models.functions import Upper

from basechem.main.models.compound_models import Series
from basechem.users.models import BasechemUser


def assays_default_value():
    """
    Default dict to use for assays JSON field
    The assays field is JSON of the form
    {
        data_exp_id: int(exp_id),   <- highest experiment ID included in a data email
        last_ping_sent: datetime, <- when the last data email was sent
        control: dn_id, <- TODO: will be used to track assay shift
        assay_exp_info: {
            assay1: [int(exp_id), exp_date:datetime] <- highest experiment ID included in a ping email, datetime of experiment
            assay2: [int(exp_id), exp_date:datetime] <- highest experiment ID included in a ping email, datetime of experiment
        }
    }
    """
    return {
        "data_exp_id": 0,
        "last_ping_sent": datetime.datetime.now(),
        "control": None,  # TODO: Will be used for tracking assay shift
        "assay_exp_info": {},
    }


class AssaysDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, assays):
        if assays.get("last_ping_sent"):
            assays["last_ping_sent"] = datetime.datetime.strptime(
                assays["last_ping_sent"], "%Y-%m-%dT%H:%M:%S.%f"
            )
        for _, data_tup in assays.get("assay_exp_info", {}).items():
            data_tup[1] = datetime.datetime.strptime(data_tup[1], "%Y-%m-%dT%H:%M:%S")
        return assays


class Project(models.Model):
    """
    Class containing project information
    """

    code = models.CharField(max_length=10)
    target = models.CharField(max_length=20)
    leads = models.ManyToManyField(BasechemUser, related_name="leads")
    subscribers = models.ManyToManyField(
        BasechemUser, related_name="subscribers", blank=True
    )
    assays = models.JSONField(
        blank=True,
        null=True,
        encoder=DjangoJSONEncoder,
        decoder=AssaysDecoder,
        default=assays_default_value,
    )
    maestro_prj = models.URLField(max_length=300, blank=True, null=True)

    class Meta:
        ordering = (Upper("code"),)

    def __str__(self):
        return self.code

    @property
    def structure_available(self):
        return (
            Series.objects.filter(project=self, active=True)
            .exclude(receptor_file_mol2="")
            .exists()
        )
