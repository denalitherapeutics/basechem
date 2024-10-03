import base64
import datetime
import io
import os

import numpy as np
from django.conf import settings
from django.core.mail import EmailMultiAlternatives, mail_admins, send_mail
from django.template import loader
from django.template.loader import render_to_string
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
from weasyprint import CSS, HTML

from basechem.common.constants import ADMIN_FAILURE
from basechem.common.dtx_utils import get_registered_structures
from basechem.main.constants import MAX_HOURS_TO_WAIT
from basechem.main.models.collection_models import Collection
from basechem.main.models.compound_models import Compound, Series
from basechem.mni_common.storage import save_media_file
from basechem.users.models import BasechemUser


def process_new_data(project, assay_data_df, analyses):
    """
    Process new assay data to send ping and/or data emails when appropriate
    :param project: the Project being processed
    :param assay_data_df: dataframe with new assay data, each line is sample in an assay run
    :param analyses: a dictionary of the form {assay_name: [analysis_names]}, the analyses
        whose data are in assay_data_df
    """
    if not assay_data_df.empty:
        # Include new series assay data for comparison in the data email
        project_series = [
            s.dn_id for s in Series.objects.filter(project=project, active=True)
        ]
        series_df = assay_data_df[assay_data_df["DN_ID"].isin(project_series)]
        if not series_df.empty:
            update_series_data(series_df)

        # Send ping emails
        send_new_data_ping(project, assay_data_df)

        # If data email is ready to send, create a collection and find DTX MMPs
        if _all_assays_in(project) or _impatient_assays(project):
            mol_blocks = get_registered_structures(
                list(set(assay_data_df["DN_ID"].to_list()))
            )
            assay_data_df["Moltext"] = assay_data_df["DN_ID"].map(mol_blocks)
            assay_data_col = create_assay_emailer_collection(
                project, assay_data_df, analyses
            )
            try:
                assayed_dns = assay_data_df["DN_ID"].unique().tolist()
                assay_data_col.find_mmps()
                assay_data_col.update_mmp_dtx_avg_assay_data(assayed_dns)
            except:
                # Don't let an MMP error stop the email from being sent
                pass

            update_assay_results_from_new_data_df(assay_data_df)
            send_new_data_assay_email(project, assay_data_df, assay_data_col.pk)


def create_assay_emailer_collection(project, assay_data, analyses):
    """
    Creates a collection from assay_data to be used for MMP.
    :param project: the Project being processed
    :param assay_data: dataframe with assay data for the project
    :param analyses: a dictionary of the form {assay_name: [analysis_names]}, the analyses
        whose data are in assay_data_df
    :return: the assay emailer collection
    """
    collection = Collection(
        project=project,
        owner=BasechemUser.objects.get(first_name="ADMIN"),
        metadata={"assays": analyses},
    )
    collection.save()
    mols = []
    for index in assay_data.index:
        mol = Chem.MolFromMolBlock(assay_data.loc[index, "Moltext"])
        if not mol:
            admin_email_message = f"The moltext from DTX API was not able to be parsed into a mol. DN: {assay_data.loc[index, 'DN_ID']} with moltext: \n {assay_data.loc[index, 'Moltext']}"
            mail_admins(ADMIN_FAILURE, admin_email_message)
            continue
        mol.SetProp("dn_id", assay_data.loc[index, "DN_ID"])
        mols.append((mol, False))
    collection.handle_romols(mols)
    return collection


def update_assay_results_from_new_data_df(assay_data_df):
    """
    Update measured_data.assay_results for Compounds in the given assay data df. This is
    called when preparing an assay emailer collection so that the newly assayed Compounds show
    the most recent assay value instead of the average values that may exist in Dotmatics.
    :param assay_data_df: dataframe with new assay data
    """
    dn_ids = assay_data_df["DN_ID"].unique().tolist()
    comps = Compound.objects.filter(dn_id__in=dn_ids)
    for comp in comps:
        assay_results = comp.measured_data.get("assay_results", {})
        for i, row in assay_data_df.loc[
            assay_data_df["DN_ID"] == comp.dn_id
        ].iterrows():
            if row["Assay_Name"] not in assay_results:
                assay_results[row["Assay_Name"]] = {}
            assay_results[row["Assay_Name"]][row["Analysis_Name"]] = row["IC50_Value"]
        comp.measured_data["assay_results"] = assay_results
        comp.save()


######################
#### PING METHODS ####
######################


def send_new_data_ping(project, assay_data_df):
    """
    Pings the leads of the project that there is new data available
    :param project: the Project being processed
    :param assay_data: dataframe with assay data for the project
    """
    need_to_send = _assays_needing_pings(project, assay_data_df)
    if need_to_send:
        assays = ", ".join(sorted(need_to_send))
        subject = f"Dotmatics got {project.code} data!"
        message = f"New data for {project.code}: {assays} has been added to Dotmatics. Check it out!"
        recipients = [lead.email for lead in project.leads.all()]
        html_message = loader.render_to_string(
            "main/components/newdata_email.html",
            {"subject": subject, "message": message, "ping": True},
        )
        send_mail(
            subject,
            message,
            settings.DEFAULT_FROM_EMAIL,
            recipients,
            fail_silently=False,
            html_message=html_message,
        )
        project.assays["last_ping_sent"] = datetime.datetime.now()


def _assays_needing_pings(project, assay_data_df):
    """
    Updates the given projects assay json field with new assays
    :param project: the Project being processed
    :param assay_data: dataframe with assay data for the project
    :return: a set of unique assays that need to send pings
    """
    # Return a list of tuples of the unique Assay_Name/Experiment_ID/Experiment_Date triples
    assays_w_exp_ids = list(
        dict(
            assay_data_df.groupby(
                ["Assay_Name", "Experiment_ID", "Experiment_Date"]
            ).size()
        ).keys()
    )
    need_to_send = set()
    # Check if each assay has already been sent recently
    for assay, new_exp_id, new_exp_date in assays_w_exp_ids:
        if assay not in project.assays["assay_exp_info"].keys() or (
            project.assays["assay_exp_info"][assay][0] < int(new_exp_id)
        ):
            project.assays["assay_exp_info"][assay] = (int(new_exp_id), new_exp_date)
            need_to_send.add(assay)

    return need_to_send


######################
#### DATA METHODS ####
######################


def send_new_data_assay_email(project, assay_data, collection_pk):
    """
    Creates and sends an email with the given assay data to project subscribers
    :param project: the Project that has new data
    :param assay_data: dataframe of the given project's new data
    :param collection_pk: the PK of the Collection created for the assayed compounds
    """
    today = datetime.datetime.now().strftime("%m_%d_%Y")
    pdf_path = f"{settings.MEDIA_ROOT}/project_{project.code}/{today}_{project.code}_assayreport.pdf"
    html_table = _create_assay_html_table(assay_data, project)
    HTML(string=html_table).write_pdf(
        pdf_path, stylesheets=[CSS("basechem/static/css/pdf.css")]
    )

    recipients = [sub.email for sub in project.subscribers.all()]
    subject = f"You've got {project.code} assay data!"
    html_content = render_to_string(
        "main/components/newdata_email.html",
        {
            "subject": subject,
            "project": project,
            "data": True,
            "collection_pk": collection_pk,
        },
    ).strip()

    msg = EmailMultiAlternatives(
        subject, subject, settings.DEFAULT_FROM_EMAIL, recipients
    )
    msg.attach_alternative(html_content, "text/html")
    msg.attach_file(pdf_path)
    msg.send()

    # Update project data to reflect email send
    project.assays["data_exp_id"] = max(
        [data_tup[0] for _, data_tup in project.assays["assay_exp_info"].items()]
    )
    project.save()

    # Delete files now that PDF has been sent
    _delete_files(project, pdf_path)


def _create_assay_html_table(assay_data, project):
    """
    Build html table from the `assay_data` dataframe
    :param assay_data: dataframe with assay data for the project
    :param project: the Project being processed
    """
    assay_data.sort_values("DN_ID", inplace=True)
    assay_data["IC50_Value"] = (
        assay_data["IC50_Value"].replace("", np.nan).astype(float)
    )

    assay_data["Structure_File"] = assay_data.apply(
        lambda row: _create_mol_image_file(row["DN_ID"], row["Moltext"], project.code),
        axis=1,
    )
    assay_data["Curve_File"] = assay_data.apply(
        lambda row: _create_curve_file(
            row["DN_ID"], row["IC50_Curve"], row["Assay_Name"], row.name, project.code
        ),
        axis=1,
    )

    assay_data = _clean_assay_names(assay_data)
    html_table = assay_data.to_html(
        columns=[
            "Structure_File",
            "DN_ID",
            "Stereo",
            "Assay_Name",
            "Classification",
            "IC50_Value",
            "Curve_File",
        ],
        index=False,
        float_format="%.3g",
        justify="center",
        classes="assay_table",
        na_rep="",
        escape=False,
    )
    return html_table


def _create_mol_image_file(dn_id, moltext, project):
    """
    Helper to create base64 image of structure
    :param dn_id: the DN ID of the compound whose mol image file is being generated
    :param mol: moltext from Dotmatics structure
    :param project: the Project being processed
    :return: an html string with an image tag pointing to the generated file
    """
    filepath = f"project_{project}/images/{dn_id}_structure.png"
    localpath = f"{settings.MEDIA_ROOT}/{filepath}"

    try:
        mol = Chem.MolFromMolBlock(moltext)
        drawmol = Draw.PrepareMolForDrawing(mol)
        os.makedirs(os.path.dirname(localpath), exist_ok=True)
        Draw.MolToFile(drawmol, localpath, size=(200, 100))
        with open(localpath, "rb") as fp:
            media_url = save_media_file(filepath, fp)

        img_tag = f"<img src={media_url} width='200px' height='131px'>"
    except:
        # We mostly see this when DTX API returns a moltext that is not a valid mol
        # Admin email sent in `create_assay_emailer_collection`
        img_tag = f"<p>Structure unable to be parsed.</p>"

    return img_tag


def _create_curve_file(dn_id, img_data, assay_name, df_index, project):
    """
    Helper to add html tag to base64 images
    :param dn_id: the DN ID of the compound whose image is being generated
    :param img_data: a base64 image of an IC50 curve for this DN ID
    :param assay_name: a string, the name of the assay
    :param df_index: an integer, the index of the row in the assay data df that is calling
        this function. This number is used to distinguish the which image belongs to which row
        when assays are run in replicate (same dn_id and assay_name, but different rows in the df).
    :param project: the Project being processed
    :return: an html string with an image tag pointing to the generated file
    """
    img_data = img_data.encode("utf_8")
    safe_assay_name = assay_name.replace(" ", "_").replace("/", "-")
    filepath = (
        f"project_{project}/images/{dn_id}_{safe_assay_name}_curve_{df_index}.png"
    )
    localpath = f"{settings.MEDIA_ROOT}/{filepath}"
    os.makedirs(os.path.dirname(localpath), exist_ok=True)

    # Save image locally
    img = Image.open(io.BytesIO(base64.decodebytes(img_data)))
    new_img = img.resize((186, 120))  # original dimensions: 320x210
    new_img.save(localpath, format="PNG")

    # Save image to S3
    with open(localpath, "rb") as fp:
        media_url = save_media_file(filepath, fp)

    img_tag = f"<img src={media_url} width='200px' height='131px'>"
    return img_tag


def _delete_files(project, pdf_path):
    """
    Helper to delete all the mol/ curve image files and pdf file
    after the pdf is sent
    :param project: the Project being processed
    :param pdf_path: path to generated pdf file
    """
    images_dir = f"{settings.MEDIA_ROOT}/project_{project}/images"
    for file in os.listdir(images_dir):
        os.remove(os.path.join(images_dir, file))

    if os.path.exists(pdf_path):
        os.remove(pdf_path)


########################
#### HELPER METHODS ####
########################


def update_series_data(series_df):
    """
    Update the activity value(s) for the series in the df
    :param series_df: a dataframe containing just series assay data
    """
    for _, row in series_df.iterrows():
        ser = Series.objects.get(dn_id=row["DN_ID"])
        if ser:
            ser.assay_data[row["Assay_Name"]] = row["IC50_Value"]
            ser.save()


def _all_assays_in(project):
    """
    Returns True if all assays have new data that has not been sent yet
    :param project: the Project object to check assays for
    """
    return not any(
        [
            data_tup[0] <= int(project.assays["data_exp_id"])
            for _, data_tup in project.assays["assay_exp_info"].items()
        ]
    )


def _any_assays_in(project):
    """
    Return True if any of the assays have new data that has not been sent yet
    :param project: the Project object to check assays for
    """
    return any(
        [
            data_tup[0] > project.assays["data_exp_id"]
            for _, data_tup in project.assays["assay_exp_info"].items()
        ]
    )


def _impatient_assays(project):
    """
    Returns True if a data email should be sent even if all assays are not in
    :param project: the Project object to check assays for
    :return: True if not all assays are in but we should send an email anyways (False otherwise)
    """
    missing_assays = not _all_assays_in(project)
    hours_since_last_ping = (
        datetime.datetime.now() - project.assays["last_ping_sent"]
    ).total_seconds() // 3600
    return (
        missing_assays
        and _any_assays_in(project)
        and hours_since_last_ping >= MAX_HOURS_TO_WAIT
    )


def depreciate_old_assays(project):
    """
    If any assays in the Project have not had new data for 90 days,
    remove from the assay field
    :param project: the Project to check assays for
    """
    assays = project.assays["assay_exp_info"]
    for name in list(assays.keys()):
        if (assays[name][1] + datetime.timedelta(days=90)) < project.assays[
            "last_ping_sent"
        ]:
            del project.assays["assay_exp_info"][name]

    project.save()


def _clean_assay_names(assay_data_df):
    """
    Shortens assay names to 3 words for the html. This function creates a new dataframe,
    because the rest of the assay emailer uses the full assay names, so we don't want
    to edit the original dataframe
    :assay_data_df: a dataframe of assay data
    :returns: a copy of assay_data_df with shortened assay names
    """
    new_assay_data_df = assay_data_df.copy()
    new_assay_data_df["Assay_Name"] = assay_data_df["Assay_Name"].apply(
        lambda x: " ".join(x.split(" ")[:3])
    )

    return new_assay_data_df
