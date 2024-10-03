import datetime

import pandas as pd
from django.conf import settings
from django_q.tasks import fetch_group

from basechem.common.inductive_utils import run_inductive_lm_predict
from basechem.main.constants import ALL_PROPS, ALOGD, HLM, RLM, ROTATABLEBONDS


def get_dtx_prop_name(prop):
    """
    Given a basechem property name, return the property name that Dotmatics expects
    :param prop: a string, the name of a property
    """
    if prop == ROTATABLEBONDS:
        return "RotBonds"
    elif prop == ALOGD:
        return "AlogD"
    elif prop == "dn_id":
        return "Compound ID"
    else:
        return prop


def generate_dtx_lm_stability_csv(collection, lm_filepath):
    """
    Generates a csv file of lm stability predictions for the compounds in the collection
    that matches what the DTX "Activity Predictions Upload" script expects
    :param collection: the collection to generate properties for
    :param lm_filepath: path to the csv file to populate
    """
    sdf_filepath, _ = collection.get_sdf_file()
    rlm_output = run_inductive_lm_predict(sdf_filepath, "R", False)
    hlm_output = run_inductive_lm_predict(sdf_filepath, "H", False)
    rlm_df = pd.DataFrame.from_dict(rlm_output)
    hlm_df = pd.DataFrame.from_dict(hlm_output)

    # model_version and latest_data_date will be the same for all entries so can just take the first
    model_version = list(set(hlm_df["model_version"]))[0]
    lm_data_date = list(set(hlm_df["latest_data_date"]))[0]

    joint_df = rlm_df.merge(hlm_df, on="name", suffixes=["_rlm", "_hlm"])
    # Add constant columns
    joint_df["assay"] = "Inductive Bio GCNN"
    joint_df["prediction_date"] = datetime.datetime.today().strftime("%m/%d/%Y")
    joint_df["model_version"] = f"{model_version}: {lm_data_date}"
    # blank column so the CSV matches the Dotmatics processing script
    joint_df["skip"] = ""

    # Add calculated columns
    joint_df["out_of_domain_flag"] = joint_df["out_of_domain_rlm"].apply(
        lambda x: "out-of-domain" if x == "True" else ""
    )
    joint_df["pStable"] = joint_df["probs_list_hlm"].apply(lambda x: f"{x[0]:.3f}")
    joint_df.to_csv(
        lm_filepath,
        columns=[
            "name",
            "assay",
            "skip",
            "skip",
            "skip",
            "out_of_domain_flag",
            "skip",
            "pStable",
            "skip",
            "skip",
            "prediction_date",
            "model_version",
            "prediction_hlm",
            "skip",
            "prediction_rlm",
            "skip",
        ],
        index=False,
    )


def generate_dtx_propcalc_csv(collection, props_filepath):
    """
    Generates a csv file of property values for the compounds in the collection
    :param collection: the collection to generate properties for
    :param props_filepath: path to the csv file to populate
    """
    props = ["dn_id"]
    props.extend(ALL_PROPS)
    props.remove(RLM)
    props.remove(HLM)
    if not settings.INDUCTIVE_BIO_ENABLED:
        props.remove(ALOGD)
    collection.metadata["props_to_show"] = props
    collection.save()
    collection.propcalc_analysis()
    group_name = collection.get_propcalc_group_name()

    # Wait indefinitely until group returns
    tasks = fetch_group(group_name, failures=True, count=collection.compounds().count())
    prop_dict = collect_propcalc_results(tasks)

    df = pd.DataFrame.from_dict(prop_dict, orient="index")
    if settings.INDUCTIVE_BIO_ENABLED:
        df[ALOGD.lower()] = df["logd_prediction"]

    # Drop columns that aren't in the properties
    cleaned_props = [prop.lower().replace(" ", "_") for prop in props]

    drop_cols = [col for col in df.columns if col not in cleaned_props]
    df = df.drop(columns=drop_cols)
    # Rename the columns to match Dotmatics names and put them in the correct order
    new_col_names = {}
    new_col_order = []
    for prop in props:
        dtx_prop_name = get_dtx_prop_name(prop)
        new_col_names[prop.lower().replace(" ", "_")] = dtx_prop_name
        new_col_order.append(dtx_prop_name)

    df.rename(columns=new_col_names, inplace=True)
    df = df.reindex(columns=new_col_order)
    df.sort_values(by=["Compound ID"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.to_csv(props_filepath, index=False)


def collect_propcalc_results(tasks):
    """
    Helper for constructing propcalc results dictionary
    :return: dictionary of properties {compound_id: dict-of-properties}
    """
    results = {}
    for task in tasks:
        co_id = int(task.name.split("_")[1])
        results[co_id] = task.result

    return results
