import datetime
import logging
import os
import re
import subprocess
from zipfile import ZipFile

import pandas as pd
from django.conf import settings
from django.core.mail import EmailMultiAlternatives, mail_admins, send_mail
from django.template import loader
from django.urls import reverse
from django_q.models import OrmQ
from django_q.tasks import async_task, fetch_group
from rdkit import Chem

from basechem.common.analysis_utils import (
    collect_align_results,
    collect_dock_results,
    collect_esp_results,
    collect_mmp_results,
    collect_torsion_results,
)
from basechem.common.analytic import TASK, Analytic
from basechem.common.constants import ADMIN_FAILURE, ADMIN_NOTIFICATION
from basechem.common.dtx_utils import (
    ASSAY_DATA_COLUMNS,
    DTX_LM_STABILITY_EXP_ID,
    DTX_LM_STABILITY_SCRIPT_ID,
    DTX_PROPCALC_EXP_ID,
    DTX_PROPCALC_SCRIPT_ID,
    get_all_dn_after,
    get_all_dtx_properties,
    get_ic50_data,
    get_logd_agg_data,
    get_registered_structures,
    upload_file_to_dtx,
)
from basechem.common.inductive_utils import update_inductive_logd_data
from basechem.common.mmpdb_utils import _fragment_mmpdb, _index_mmpdb, loadprops_mmpdb
from basechem.common.propcalc_utils import (
    collect_propcalc_results,
    generate_dtx_lm_stability_csv,
    generate_dtx_propcalc_csv,
)
from basechem.main.assay_emailer_utils import depreciate_old_assays, process_new_data
from basechem.main.constants import ALIGN, AUTO, DOCK, ESP, MMP, PROPCALC, TORSION
from basechem.main.models.collection_models import Collection, collection_files_path
from basechem.main.models.compound_models import CompoundOccurrence
from basechem.main.models.project_models import Project
from basechem.users.models import BasechemUser

logger = logging.getLogger("django")


##########################
# Scheduled Task Methods #
##########################


def check_for_new_assay_data():
    """
    Start async_tasks for each project to poll new assay data, if there is no current task running
    """
    projects = Project.objects.exclude(assays={})
    for project in projects:
        task_name = f"check_for_{project.code}_assay_data"
        queued = [oq for oq in OrmQ.objects.all() if oq.task.get("name") == task_name]
        if not queued:
            async_task(
                check_for_new_assay_data_by_project,
                project=project,
                task_name=task_name,
            )


def check_for_new_assay_data_by_project(project):
    """
    Poll Dotmatics for new assay data between now and the last time an assay email was sent for the given project
    :param project: a Project object, what project should be queried
    """
    depreciate_old_assays(project)

    assay_data_in_range, analyses = get_ic50_data(
        project.assays["data_exp_id"], project.target
    )
    if assay_data_in_range:
        assay_data_as_series = [
            (k, *tup) for k, val in assay_data_in_range.items() for tup in val
        ]
        all_data_df = pd.DataFrame(assay_data_as_series, columns=ASSAY_DATA_COLUMNS)
        process_new_data(project=project, assay_data_df=all_data_df, analyses=analyses)

        project.save()


def dtx_propcalc(date=None):
    """
    Gets all DNs registered on a date and performs the following analyses on all compounds:
    1. Runs property calculator
    2. Requests predictions from the InductiveBio RLM/HLM models
    For each of these analyses, a csv file is generated with results and is posted to Dotmatics
    using the appropriate experiment and script IDs.
    :param date: (optional) a datetime object, use all compounds registered on the given date,
        by default, the date is yesterday
    :return: a boolean, was this task successful
    """
    if not date:
        date = datetime.datetime.today() - datetime.timedelta(days=1)
    dn_ids = get_all_dn_after(date=date)
    if not dn_ids:  # No dns registered on date, stop here
        return True
    mol_blocks = get_registered_structures(dn_ids)
    mols = []
    for dn_id, mol_block in mol_blocks.items():
        try:
            mol = Chem.MolFromMolBlock(mol_block)
            mol.SetProp("dn_id", dn_id)
            mols.append((mol, True))  # Always 2D coming from DTX
        except:
            logging.error(f"Properties could not be set for: {dn_id}")

    # Set up collection
    collection = Collection()
    collection.project = Project.objects.get(code=AUTO)
    collection.owner = BasechemUser.objects.get(first_name="ADMIN")
    collection.save()
    collection.handle_romols(mols)
    ############################################
    # Make and upload csv to propcalc protocol #
    ############################################
    file_suffix = f"{date.strftime('%m-%d-%Y')}_{dn_ids[-1]}-{dn_ids[0]}.csv"
    props_filename = f"basechem_propcalc_{file_suffix}"
    props_filepath = collection_files_path(collection, props_filename, local=True)
    generate_dtx_propcalc_csv(collection, props_filepath)

    props_file = open(props_filepath, "rb")
    upload_file_to_dtx(
        props_file, props_filename, DTX_PROPCALC_EXP_ID, DTX_PROPCALC_SCRIPT_ID
    )

    ################################################
    # Make and upload csv to LM stability protocol #
    ################################################
    lm_filename = f"basechem_lm_stability_{file_suffix}"
    lm_filepath = collection_files_path(collection, lm_filename, local=True)
    try:
        generate_dtx_lm_stability_csv(collection, lm_filepath)
    except:
        message = f"Failed to generate LM stability CSV: {lm_filename}"
        mail_admins(ADMIN_FAILURE, message)
    lm_file = open(lm_filepath, "rb")
    upload_file_to_dtx(
        lm_file, lm_filename, DTX_LM_STABILITY_EXP_ID, DTX_LM_STABILITY_SCRIPT_ID
    )

    collection.delete()


def update_logd_model_data():
    """
    Polls Dotmatics for new logD data from the past week, sends to Inductive
    and checks current prediction accuracy
    """
    last_week = datetime.datetime.today() - datetime.timedelta(days=8)
    last_week_fmt = last_week.strftime("%Y%m%d")
    logd_data_mols = get_logd_agg_data(last_week_fmt)

    if logd_data_mols:
        tmp_data_file = f"/tmp/{last_week_fmt}_logd_dtx.sdf"
        writer = Chem.SDWriter(tmp_data_file)
        for mol in logd_data_mols:
            writer.write(mol)
        writer.close()

        # Extra check to not send InductiveBio DTX TEST data
        upload_response = ""
        if "test" not in settings.DTX_HOST:
            upload_response = update_inductive_logd_data(tmp_data_file)

        return upload_response


def update_mmpdb(dn_id=None):  # pragma: no cover
    """
    Pull compounds from DTX registered after the given DN ID and create a new MMPDB
    that includes these compounds
    :param dn_id: the dn_id after which to add new compounds
    """
    if not dn_id:
        cached_smi = f"{settings.MMPDB_DIR}/mmpdb_data/cache.smi"
        # Get last line in file with tail so whole file isn't read in
        last_line = subprocess.check_output(["tail", "-1", cached_smi])
        dn_id = re.search("DN\d{7}", str(last_line)).group(0)

    if dn_id:
        new_dn_ids = get_all_dn_after(dn_id=dn_id)
        if not new_dn_ids:
            mail_admins(
                ADMIN_NOTIFICATION,
                f"MMPDB was not updated due to no DNs found in DTX since {dn_id}.",
            )
            return

        mol_blocks = get_registered_structures(new_dn_ids)
        mols = []
        for dn, m in mol_blocks.items():
            mol = Chem.MolFromMolBlock(m)
            mol.SetProp("_Name", dn)
            mols.append(mol)

        smiles = [(Chem.MolToSmiles(mol), mol.GetProp("_Name")) for mol in mols]
        smiles.sort(key=lambda x: x[1])

        # Append to the existing smiles file on EFS so we don't need to pull all data
        with open(cached_smi, "a") as fp:
            for tup in smiles:
                fp.write(f"{tup[0]}   {tup[1]}\n")

        # Generate fragments and index the db for the latest data
        frag_db = _fragment_mmpdb(cached_smi)
        if not frag_db:
            return
        indexed = _index_mmpdb(frag_db)
        if not indexed:
            return

        mail_admins(
            ADMIN_NOTIFICATION,
            f"MMPDB successfully updated with {len(new_dn_ids)} new compounds: {smiles[0][1]} to {smiles[-1][1]}.",
        )

    # Whenever new compounds are added, also recompute property statistics
    update_mmpdb_props()


def update_mmpdb_props():  # pragma: no cover
    """
    Pull the DENALI_COMPOUND_ASSAY_PROPS datasource from DTX and update MMPDB
    with those property values
    """
    propfile_path = f"{settings.MMPDB_DIR}/mmpdb_data/cache.props"
    all_cpds = get_all_dtx_properties()

    all_df = pd.DataFrame.from_dict(all_cpds, orient="index")
    data_df = all_df.loc[all_df["data"] != {}]
    data = data_df["data"].to_list()
    results_df = pd.DataFrame(data)

    first_column = results_df.pop("COMPOUND_ID")
    results_df.insert(0, "ID", first_column)
    results_df = results_df[["ID", "LOGD_AVG", "PIC50"]]
    # mmpdb needs empty values to be *
    results_df.fillna("*", inplace=True)

    results_df.to_csv(propfile_path, sep="\t", index=False)

    loaded = loadprops_mmpdb(propfile_path)
    if loaded:
        mail_admins(
            ADMIN_NOTIFICATION,
            f"MMPDB loadprops succeeded. {len(results_df)} were added to mmpdb with the columns: {list(results_df.columns)}",
        )


def monitor_toklat_scoring():
    """
    This method runs rDock for all Compounds that were uploaded by users in the past 7 days
    and determines if Toklat picked different poses than rDock would have. Sends a report
    to the Django admins with the results.
    """
    co_dict = {"0": [], "1": [], "2": [], "3": [], "skipped": []}
    file_dict = {"0": [], "1": [], "2": [], "3": []}
    errors = []
    # Get compound occurrences uploaded in the last week
    today = datetime.datetime.today()
    collections = Collection.objects.filter(
        created_on__gte=today - datetime.timedelta(days=7)
    )
    cos = (
        CompoundOccurrence.objects.filter(collection__in=collections)
        .order_by("compound__pk")
        .distinct()
    )
    if not cos:
        return
    for co in cos:
        try:
            co.admin_url = settings.BASE_URL + reverse(
                "admin:main_compoundoccurrence_change", args=(co.pk,)
            )
            # Perform docking analysis
            reference = co.compound.series
            _, ref_id = co.compound._pick_reference_file_and_id(reference)
            co.dock_to_receptor(reference=reference, torsion_alerts=True, toklat=True)
            rdock_filepath = co._get_rdock_output_path(reference)
            # If mayachem fails for some reason, the toklat path name will be different
            toklat_filepath = (
                f"{os.path.splitext(rdock_filepath)[0]}_alerts_all_scored.sdf"
            )
            if not os.path.exists(toklat_filepath):
                toklat_filepath = f"{os.path.splitext(rdock_filepath)[0]}_scored.sdf"
            mols = [mol for mol in Chem.SDMolSupplier(toklat_filepath, removeHs=False)]
            for i, mol in enumerate(mols, start=1):
                conf_id = mol.GetProp("s_conf_id")
                pose_id = f"c{co.compound.pk}-co{co.pk}-{conf_id}-{i}-{ref_id}"
                mol.SetProp("s_pose_id", pose_id)

            sorted_by_toklat = sorted(
                mols, key=lambda x: float(x.GetProp("toklat_score")), reverse=False
            )[:3]
            sorted_by_rdock = sorted(
                mols, key=lambda x: float(x.GetProp("SCORE")), reverse=False
            )[:3]

            toklat_poses = set([mol.GetProp("s_pose_id") for mol in sorted_by_toklat])
            rdock_poses = set([mol.GetProp("s_pose_id") for mol in sorted_by_rdock])

            num_diff = len(toklat_poses.difference(rdock_poses))
            co_dict[str(num_diff)].append(co)
            file_dict[str(num_diff)].append(toklat_filepath)
        except Exception as e:
            co_dict["skipped"].append(co)
            errors.append(f"Error running dock for {co.compound.name}: {e}")

    recipients = [email for _, email in settings.ADMINS]
    subject = f"Toklat vs rDock Scoring Comparison ({today.strftime('%m-%d-%Y')})"
    msg = EmailMultiAlternatives(
        subject, subject, settings.DEFAULT_FROM_EMAIL, recipients
    )
    # Zip files and attach them to email
    for num_diff, filepaths in file_dict.items():
        try:
            if filepaths:
                zip_filename = f"toklat_compare_{num_diff}_diff.zip"
                with ZipFile(zip_filename, "w") as zfp:
                    for filepath in filepaths:
                        _, filename = os.path.split(filepath)
                        zfp.write(filepath, arcname=filename)
                msg.attach_file(zip_filename)
        except Exception as e:
            errors.append(f"Error zipping files for {num_diff} differences: {e}")
    # Add HTML content
    html_content = loader.render_to_string(
        "main/components/toklat_comparison_email.html",
        {
            "compound_occurrences": co_dict,
            "errors": errors,
            "total_num_cos": cos.count(),
        },
    ).strip()
    msg.attach_alternative(html_content, "text/html")
    msg.send()


#########################
# DjangoQ Queue Methods #
#########################


def get_group_results(group_name):
    """
    Aggregates results of the tasks in a group. This is also called in the
    ajax_view to check when results are ready before reloading the analysis page.
    :param group_name: the name of the group to collect results from
    :return: a tuple - (failed, completed, results)
        - [0] failed: did all of the tasks in the group fail
        - [1] completed: True if all tasks are complete
        - [2] result: dict result of the async task group (None if a task in the group is still running)
    """
    tasks = fetch_group(group_name, failures=True)
    failed = not any(task.success for task in tasks) if tasks else None
    queued = [oq for oq in OrmQ.objects.all() if oq.task.get("group") == group_name]

    # If a job is restarted due to failure, it will be in the queue as well as failed jobs.
    # So even if the number of tasks in the group matches the number of expected tasks,
    # we assume that the group has not completed all tasks while at least one of the tasks
    # is in the queue. The failed job is removed when the restarted job succeeds.
    # If all retries fail, the failed job stays such that `count_group` returns the
    # expected number of tasks.
    if tasks and not queued:
        completed = all([task.stopped for task in tasks])
        if completed:
            if PROPCALC in group_name:
                results = collect_propcalc_results(tasks)
            elif ALIGN in group_name:
                results = collect_align_results(tasks)
            elif DOCK in group_name:
                results = collect_dock_results(tasks)
            elif ESP in group_name:
                results = collect_esp_results(tasks)
            elif TORSION in group_name:
                results = collect_torsion_results(tasks)
            elif MMP in group_name:
                results = collect_mmp_results(tasks)

            return failed, True, results

    return failed, False, None


def task_completion_hook(task):
    """
    A wrapper that triggers actions when a django q task completes
    :param task: a djangoQ task, this is generally called as a hook
    """
    if not task.group:
        # This hook is being run because the task group was deleted: don't send any emails
        return
    failed, _, _ = get_group_results(task.group)

    # While this hook is running, `task` is still on the queue even though it is complete
    # (because it's running this hook). `task`'s group can be considered complete if `task`
    # is the only task from this group on the queue
    queued = [oq for oq in OrmQ.objects.all() if oq.task.get("group") == task.group]
    completed = False

    if not queued or len(queued) == 1 and queued[0].task_id() == task.id:
        completed = True

    # Determine how long tasks took to run
    tasks = fetch_group(task.group, failures=True)
    duration = max([task.time_taken() for task in tasks])

    # Send notifications to requestors/error notifications to admins
    send_task_notifications(task, failed, completed, duration)
    if completed:
        # Log the completed task
        log_taskgroup_analytics(tasks, duration)


def simple_email_task_hook(task):
    """
    A wrapper that triggers a completion or failure admin email notification when a django q task completes
    :param task: djangoQ task, this is called as a hook
    """
    time_taken = round(task.time_taken())

    if task.success and len(task.result) == 2:
        outfile, errors = task.result
        if os.path.exists(outfile):
            task_url = reverse("admin:django_q_success_change", args=(task.id,))
            message = f"The task {task.name} ({settings.BASE_URL}{task_url}) successfully completed in {time_taken} seconds with {len(errors)} errors. The results file can be found at {outfile}"
            mail_admins(ADMIN_NOTIFICATION, message)
    else:
        task_url = reverse("admin:django_q_failure_change", args=(task.id,))
        message = f"The task {task.name} ({settings.BASE_URL}{task_url}) failed in {time_taken} seconds. These errors exist: {task.result}"
        mail_admins(ADMIN_FAILURE, message)


def log_taskgroup_analytics(tasks, duration):
    """
    Logs the duration of a group of tasks to AWS for analytics. One statement is logged
    per group by only logging when the last task of the group completes
    :param tasks: a queryset of all tasks in the the task group
    :param duration: the number of seconds this task group took to complete
    """
    num_failed = sum([not task.success for task in tasks])
    collection = Collection.objects.get(id=tasks[0].group.split("_")[1])
    num_comps = collection.compounds().count()
    task_type = tasks[0].group.split("_")[0]
    Analytic(
        TASK,
        collection.project.code,
        collection.owner,
        task_type=task_type,
        duration=duration,
        num_comps=num_comps,
        num_failed=num_failed,
    )


def send_task_notifications(task, failed, completed, duration, duration_threshold=60):
    """
    Sends email notifications related to tasks after a task finishes. The emails sent are:
        - if the task failed and the taskgroup is not complete, sends an email to django admins
        - if the taskgroup is complete with no failures, sends an email to the collection owner
        - if the taskgroup is complete with failures, sends an email to the collection owner and django admins
    :param task: a djangoQ task, this is generally called as a hook
    :param failed: a boolean, did all tasks in this taskgroup fail
    :param completed: a boolean, are all tasks in this taskgroup finished
    :param duration: the number of seconds this task group took to complete
    :param duration_threshold: the minimum task duration that should trigger an email notification (in seconds)
    """
    # Not all tasks have completed, still send an email to django admins if there's a failure
    collection = Collection.objects.get(id=task.group.split("_")[1])
    owner = collection.owner
    if not task.success:
        # Send admins an email immediately if this task but not all tasks failed
        time_taken = round(task.time_taken())
        collection_url = reverse("admin:main_collection_change", args=(collection.id,))
        task_url = reverse("admin:django_q_failure_change", args=(task.id,))
        message = f"The task {task.name} ({settings.BASE_URL}{task_url}) failed in {time_taken} seconds. The Collection ({settings.BASE_URL}{collection_url}) belongs to {owner}."
        mail_admins(ADMIN_FAILURE, message)

    # Send email to the collection owner if all tasks are complete and they took longer than the threshold
    if completed and duration > duration_threshold:
        analysis_type = task.group.split("_")[0]

        analysis_dict = {
            PROPCALC: "property calculation",
            ALIGN: "ligand alignment",
            DOCK: "compound docking",
            ESP: "ESP map generation",
            TORSION: "torsion scan",
            MMP: "MMP analysis",
        }
        recipients = [collection.owner.email]
        url = collection.get_url(analysis_type, task.group)
        if not failed:
            subject = "Basechem Task Complete"
            message = f"Your {analysis_dict.get(analysis_type)} for collection {collection.id} is now complete."
            plain_text_message = (
                f"{message}\nClick this link to see your results: {url}"
            )
        else:
            # Users/Admins will only get this email if all tasks in the group have failed
            subject = "Basechem Task Failed"
            message = f"Your {analysis_dict.get(analysis_type)} for collection {collection.id} failed. Admins have been notified of this failure."
            plain_text_message = f"{message}\nClick this link to go to the align page for this collection: {url}"
            recipients.extend([email for admin, email in settings.ADMINS])

        html_message = loader.render_to_string(
            "main/components/analysis_email.html",
            {"subject": subject, "message": message, "url": url, "success": not failed},
        )
        send_mail(
            subject,
            plain_text_message,
            settings.DEFAULT_FROM_EMAIL,
            recipients,
            fail_silently=False,
            html_message=html_message,
        )


def hide_outdated_results(collection, group_name, results):
    """
    If the currently saved task results for `group_name` don't include the necessary
    CompoundOccurrences, hide the existing results so that new tasks are run.
    :param collection: the Collection running this torsion analysis
    :param group_name: a string, the name of the async group for the analysis
    :param results: a dictionary, the currently saved results dictionary for this group
    :returns: a boolean, was the task group hidden?
    """
    analysis = group_name.split("_")[0]
    cos = collection.get_cos_for_analysis(analysis)
    # All COs we want included in the analysis
    desired_set = set([f"co-{co.pk}" for co in cos])
    # All COs that are in the cached set
    if TORSION in group_name:
        results_set = set(results.keys())
    elif ESP in group_name:
        results_set = set(results["compounds"].keys())
    if not desired_set.issubset(results_set):
        # If the keys we expect are not all present in the current results, hide the tasks
        # so the analysis will re-run
        hide_task_group(group_name)
        return True
    return False


def get_analysis_results(collection, current_view, group_name, run_analysis_kwargs):
    """
    Either creates a new async task group or gets the result of an existing async
    group task for the given analysis, this is called in get_context_data in each
    analysis view
    :param collection: the collection being analyzed
    :param current_view: the analysis view/type of analysis
    :param group_name: the name of the async group for the analysis
    :param run_analysis_kwargs: a dictionary of kwargs for the `run_analysis` function
    :return: a tuple - (failed, result)
        - [0] failed: did any of the tasks in the group fail
        - [1] result: dict result of the async task group (None if task is still running)
    """
    tasks = fetch_group(group_name, failures=True)
    # fetch_group doesn't checked queued tasks, only successful/failed
    queued = [oq.task.get("group") == group_name for oq in OrmQ.objects.all()]
    if not tasks and not any(queued):
        collection.run_analysis(current_view, **run_analysis_kwargs)

    failed, _, result = get_group_results(group_name)

    if not failed and result and current_view in [ESP, TORSION]:
        # Ensure that the task results include the expected CompoundOccurrences
        rerun = hide_outdated_results(collection, group_name, result)
        if rerun:
            return get_analysis_results(
                collection, current_view, group_name, run_analysis_kwargs
            )

    return failed, result


def hide_task_group(group_name):
    """
    Given a group name, edit all tasks in the group to have no group name. This is used when
    the task results do not include the CompoundOccurrences that are expected (because users
    changed their saved conformers). "Hiding" the tasks instead of deleting them allows the
    saved results to still be used when possible.
    :param group_name: the name of the async group to be "hidden"
    """
    tasks = fetch_group(group_name, failures=True)
    if tasks:
        for task in tasks:
            task.group = None
            task.save()
