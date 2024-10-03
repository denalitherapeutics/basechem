import requests
from django.conf import settings
from django.core.mail import mail_admins
from rdkit import Chem

from basechem.common.constants import ADMIN_FAILURE

##################
#   LM METHODS   #
##################


def run_inductive_lm_predict(
    input_file, species, image=True, version=settings.INDUCTIVE_VERSION
):
    """
    Send the given sdf file to the InductiveBio predict endpoint
    to retrieve LM stability predictions for Human or Rat
    :param input_file: sdf file containing compounds of which to predict LM stability
    :param species: H or R for which model to predict (human or rat)
    :param image: boolean for if images should be requested from the API
    :param version: either "latest" or "dev"
    :return: returns a list of dicts of the form
        {
        "name": str,
        "prediction": str(number),
        "measured": str(number) or "",
        "out_of_domain": str(bool),
        "latest_data_date": "YYYY-MM-DD",
        "model_version": str,
        "interp_image": str,
        "probs_image": str,
        "probs_list": list,
    }
    """
    url = "https://api.inductive.bio/0.1/properties/predict"
    headers = {
        "X-API-KEY": settings.INDUCTIVE_API_KEY,
        "X-CUSTOMER-ID": settings.INDUCTIVE_CUSTOMER_ID,
    }

    with open(input_file) as f:
        sdf_text = f.read()

        suppl = Chem.SDMolSupplier(input_file)
        ids = [m.GetProp("_Name") for m in suppl]

        body = {
            "model_id": f"denali_{species.lower()}lm",
            "model_version": version,
            "input_type": "SDF",
            "input": sdf_text,
            "molecule_ids": ids,
            "include_interpretations": image,
        }

        response = requests.post(url, json=body, headers=headers)
        response_dict = response.json()
        if not isinstance(response_dict, dict):
            response_dict = {}

        predictions = _parse_inductive_lm_response(response_dict)
        return predictions


def _parse_inductive_lm_response(response_dict):
    """
    Parse the InductiveBio predict response into Basechem preferred format
    :param response_dict: json response from Inductive predict
    :return: parsed response dict
    """
    response_preds = response_dict.get("predictions", {})
    if not response_preds:
        admin_email_message = f"InductiveBio API failed to generate LM predictions with a response of: \n {response_dict}"
        mail_admins(ADMIN_FAILURE, admin_email_message)

        return [
            {
                "name": "",
                "prediction": "InductiveBio API Error - LMs",
                "measured": "",
                "out_of_domain": "",
                "latest_data_date": "",
                "model_version": "",
                "interp_image": "",
                "probs_image": "",
                "probs_list": [],
            }
        ]

    predictions = []
    for cpd_data in response_preds:
        cpd_pred = {}
        cpd_pred["name"] = cpd_data["molecule_id"]
        cpd_pred["prediction"] = str(round(cpd_data["continuous_prediction"]))
        measured = cpd_data.get("measured_value", {}).get("continuous_value", "")
        if measured != "":
            measured = str(round(measured))
        cpd_pred["measured"] = measured
        cpd_pred["out_of_domain"] = str(cpd_data["out_of_domain_flag"])

        # Get base64 images if they exist
        interp_image = ""
        probs_image = ""
        if cpd_data.get("interpretation_images"):
            interp_image = cpd_data.get("interpretation_images")[0]["image"]
            probs_image = cpd_data.get("interpretation_images")[1]["image"]
        cpd_pred["interp_image"] = interp_image
        cpd_pred["probs_image"] = probs_image
        cpd_pred["probs_list"] = cpd_data.get("predicted_probabilities", [])

        cpd_pred["latest_data_date"] = response_dict.get("latest_data_date")
        cpd_pred["model_version"] = response_dict.get("model_version")

        predictions.append(cpd_pred)

    return predictions


##################
#  LOGD METHODS  #
##################


def run_inductive_alogd_predict(input_file, version=settings.INDUCTIVE_VERSION):
    """
    Send the given sdf file to the InductiveBio predict endpoint to retrieve LogD predictions
    :param input_file: sdf file containing compounds of which to predict LogD
    :param version: either "latest" or "dev"
    :return: returns a list of dicts of the form
    {
        "name": str,
        "prediction": str(number),
        "latest_data_date": "YYYY-MM-DD",
        "model_version": str,
    }
    """
    url = "https://api.inductive.bio/0.1/properties/predict"
    headers = {
        "X-API-KEY": settings.INDUCTIVE_API_KEY,
        "X-CUSTOMER-ID": settings.INDUCTIVE_CUSTOMER_ID,
    }

    with open(input_file) as f:
        sdf_text = f.read()

        suppl = Chem.SDMolSupplier(input_file)
        ids = [m.GetProp("_Name") for m in suppl]

        body = {
            "model_id": "denali_logd",
            "model_version": version,
            "input_type": "SDF",
            "input": sdf_text,
            "molecule_ids": ids,
        }

        response = requests.post(url, json=body, headers=headers)
        response_dict = response.json()
        if not isinstance(response_dict, dict):
            response_dict = {}

        predictions = _parse_inductive_logd_response(response_dict)
        return predictions


def _parse_inductive_logd_response(response_dict):
    """
    Parse the InductiveBio LogD predict response into Basechem preferred format
    :param response_dict: json response from Inductive predict
    :return: parsed response dict
    """
    response_preds = response_dict.get("predictions", {})
    if not response_preds:
        admin_email_message = f"InductiveBio API failed to generate aLodD predictions with a response of: \n {response_dict}"
        mail_admins(ADMIN_FAILURE, admin_email_message)

        return [
            {
                "name": "",
                "prediction": "InductiveBio API Error - LogD",
                "model_version": "",
                "latest_data_date": "",
            }
        ]

    predictions = []
    for cpd_data in response_preds:
        cpd_pred = {}
        cpd_pred["name"] = cpd_data["molecule_id"]
        cpd_pred["prediction"] = str(round(cpd_data["continuous_prediction"], 2))
        cpd_pred["model_version"] = response_dict.get("model_version")
        cpd_pred["latest_data_date"] = response_dict.get("latest_data_date")

        measured = cpd_data.get("measured_value", {}).get("continuous_value", "")
        if measured:
            measured = str(round(measured, 2))
        cpd_pred["measured"] = measured

    predictions.append(cpd_pred)

    return predictions


def update_inductive_logd_data(new_data_filepath):
    """
    PUT given file with logD data to Inductive endpoint
    :param new_data_filepath: path to file containing new data
    """
    url = "https://api.inductive.bio/0.1/data/sdf_upload"
    headers = {
        "X-API-KEY": settings.INDUCTIVE_API_KEY,
        "X-CUSTOMER-ID": settings.INDUCTIVE_CUSTOMER_ID,
    }

    with open(new_data_filepath) as f:
        sdf_text = f.read()

        if sdf_text:
            logd_body = {
                "data": sdf_text,
                "project_column": "project_code",
                "molecule_id_column": "dn_id",
                "properties_dict": {
                    "logd_avg": {
                        "date_column": "most_recent_date",
                        "model_ids": ["denali_logd"],
                    }
                },
            }

            response = requests.put(url, json=logd_body, headers=headers)
            r = response.json()
            if r.get("status") != "SUCCESS":
                admin_email_message = f"InductiveBio API failed to PUT new LogD data for the file: {new_data_filepath} with the response {r}"
                mail_admins(ADMIN_FAILURE, admin_email_message)
            elif r.get("status") == "SUCCESS":
                return "SUCCESS"
