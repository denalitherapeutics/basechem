import datetime

from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand
from django_q.models import Schedule

from basechem.common.mmpdb_utils import initialize_mmpdb
from basechem.main.constants import AUTO, TEST
from basechem.main.models.project_models import Project


class Command(BaseCommand):
    help = "Add base admin, TEST project, and output dirs upon start up if they don't exist"

    def handle(self, *args, **options):
        Users = get_user_model()
        if not Users.objects.filter(last_name="ADMIN").exists():
            Users.objects.create_superuser(
                first_name="ADMIN",
                last_name="ADMIN",
                username=settings.ADMIN_USER,
                password=settings.ADMIN_PASSWORD,
                email=settings.ADMIN_EMAIL,
            )

        # Set up project codes
        if settings.ENVIRONMENT == "local":
            if not Project.objects.filter(code=TEST).exists():
                Project.objects.create(code=TEST, target=TEST)

        if not Project.objects.filter(code=AUTO).exists():
            Project.objects.create(code=AUTO, target=AUTO)

        if not settings.ENVIRONMENT or settings.ENVIRONMENT != "prod":
            initialize_mmpdb("test")

        #######################
        ### SCHEDULED TASKS ###
        #######################

        # Slurm JWT tokens are generated weekly since they expire every 999999 seconds
        name = "Generate Slurm JWT Token"
        if settings.SLURM_REST_API_HOST:
            if not Schedule.objects.filter(name=name).exists():
                Schedule.objects.create(
                    name=name,
                    func="basechem.common.slurm_utils.generate_slurm_api_token",
                    schedule_type=Schedule.WEEKLY,
                    kwargs={"q_options": {"broker_name": "default"}},
                )
        else:
            Schedule.objects.filter(name=name).delete()

        if settings.ENVIRONMENT != "local" and settings.DTX_HOST:
            name = "DTX Propcalc"
            if not Schedule.objects.filter(name=name).exists():
                Schedule.objects.create(
                    name=name,
                    func="basechem.main.tasks.dtx_propcalc",
                    schedule_type=Schedule.DAILY,
                    kwargs={"q_options": {"broker_name": "default"}},
                )

            name = "Update aLogD Model Data"
            if not Schedule.objects.filter(name=name).exists():
                Schedule.objects.create(
                    name=name,
                    func="basechem.main.tasks.update_logd_model_data",
                    schedule_type=Schedule.WEEKLY,
                    kwargs={"q_options": {"broker_name": "default"}},
                )

            name = "Assay Emailer"
            if not Schedule.objects.filter(name=name).exists():
                Schedule.objects.create(
                    name=name,
                    func="basechem.main.tasks.check_for_new_assay_data",
                    schedule_type=Schedule.CRON,
                    # run every 15min, Wednesday - Friday
                    cron="*/15 * * * wed-fri",
                    kwargs={"q_options": {"broker_name": "default"}},
                )

            name = "Update MMPDB with DTX Data"
            if not Schedule.objects.filter(name=name).exists():
                Schedule.objects.create(
                    name=name,
                    func="basechem.main.tasks.update_mmpdb",
                    schedule_type=Schedule.MONTHLY,
                    kwargs={"q_options": {"broker_name": "default"}},
                )

            name = "Monitor Toklat scoring"
            today = datetime.date.today()
            if not Schedule.objects.filter(name=name).exists():
                Schedule.objects.create(
                    name=name,
                    func="basechem.main.tasks.monitor_toklat_scoring",
                    schedule_type=Schedule.WEEKLY,
                    kwargs={"q_options": {"broker_name": "slow"}},
                    # First run should occur at 11:00 PM next Friday
                    next_run=datetime.datetime.combine(
                        today + datetime.timedelta(days=(4 - today.weekday())),
                        datetime.time(hour=23, minute=0),
                    ),
                )
