import datetime
import logging
from collections import defaultdict

from django.conf import settings
from rdkit import Chem

#####################
###   CONSTANTS   ###
#####################

if "test" in settings.DTX_HOST:
    DTX_PROPCALC_EXP_ID = 256179
else:
    DTX_PROPCALC_EXP_ID = 272383
DTX_PROPCALC_SCRIPT_ID = 28106

DTX_LM_STABILITY_SCRIPT_ID = 14143
if "test" in settings.DTX_HOST:
    DTX_LM_STABILITY_EXP_ID = 268145
else:
    DTX_LM_STABILITY_EXP_ID = 269225

ASSAY_DATA_COLUMNS = [
    "DN_ID",
    "Stereo",
    "Target",
    "Experiment_ID",
    "Experiment_Date",
    "Assay_Name",
    "Analysis_Name",
    "Classification",
    "IC50_Value",
    "IC50_Curve",
]

logger = logging.getLogger("django")

# Denali Basechem makes use of an internal python library `dtxwrapper` that makes API calls
# to the Dotmatics API. The methods in the `try` block are internal implementations, methods
# in the `except` block are stubs that can be implemented based on a different organization's
# chemical registry. They may also be left as is, in which case Basechem will still operate
# but the assay emailer will not and DN ID will always be blank.
try:
    # If the required settings are blank, don't try to use Dotmatics API
    if not all([settings.DTX_HOST, settings.DTX_PASSWORD, settings.DTX_USER]):
        raise ImportError

    #################################
    #      DENALI IMPLEMENTATIONS      #
    #################################
    from dtxwrapper.browser import (
        get_parsed_den_studies_charts_merged_v6,
        get_parsed_den_studies_ic50_agg_merged_v2,
        get_parsed_denali_compound_assay_props,
        get_parsed_denali_logd_agg,
        get_parsed_reg_data_vw,
        get_parsed_structure,
        query_and_parse_dtx,
        query_and_parse_reg_data_vw_by_dn_id,
        query_and_parse_reg_data_vw_by_inchi,
        query_and_parse_reg_data_vw_by_reg_date,
        upload_file_to_dtx_experiment,
    )
    from dtxwrapper.browser.denali_cpds.den_studies_charts_merged_v6 import (
        DS_ID as DTX_CHARTS_V6_DS_ID,
    )
    from dtxwrapper.common import construct_query_data

    #################
    ###   UTILS   ###
    #################

    def check_dtx_for_inchi(inchi):
        """
        Query Dotmatics to determine if a compound already exists in registry
        :param inchi: the inchi of the molecule (generated with the '/suu' flag)
        :returns: a string, the DN ID of the matching compound
        """
        dn_ids = query_and_parse_reg_data_vw_by_inchi(inchi)
        if not dn_ids:
            logger.info(f"No Denali compounds found matching inchi: {inchi}")
            return ""
        elif len(dn_ids) == 1:
            return dn_ids[0]
        else:
            logger.info(
                f"Multiple Denali compounds ({', '.join(dn_ids)}) found matching inchi: {inchi}"
            )
            return choose_dn_id_by_stereo(dn_ids)

    def choose_dn_id_by_stereo(dn_ids):
        """
        Given a list of DN IDs that match a single compound, choose which one(s) should
        be assigned in basechem based on the stereo of the DN IDs. The rules are:
        - if any are a "mixture of diastereomers" or "mixture of enantionmers", assign all the "mixture" DNs
        - if only one "single unknown enantiomer" assign that DN
        - else, pick nothing
        :dn_ids: a list of DN IDs (strings of the form 'DN0000001')
        :returns: a string of DN IDs separated by a comma (usually only 1, can be multiple)
        """
        reg_data = get_parsed_reg_data_vw(dn_ids)

        mixtures, unknown = [], []
        for dn_id, data in reg_data.items():
            stereo = data["data"]["STEREO_COMMENTS"].lower()
            if "unknown" in stereo:
                unknown.append(dn_id)
            elif "mixture" in stereo:
                mixtures.append(dn_id)

        if mixtures:
            return ",".join(sorted(mixtures))
        elif len(unknown) == 1:
            return unknown[0]
        else:
            return ""

    def get_ic50_data(last_data_exp_id, target):
        """
        Wrapper to get the assay data of DNs for the given target that have new results
        :param last_data_exp_id: a 6 digit ELN id of the last sent data
        :param target: a Project's target name (as used in Dotmatics)
        :returns: a tuple (results, analyses) where
            - [0] results is a dictionary of tuples base64 encoded ic50 curves that looks like:
            {dn_id: [(stereochemistry, target, experiment id, experiment datetime object, assay name, ic50 value, base64 encoded ic50 curve)]}
            - [1] analyses is a dictionary of the form {assay_name: [analysis_names]} containing all the
            analyses whose data are in results
        """
        all_results = {}
        analyses = defaultdict(list)
        # Query DTX for DN IDs with recent data for the target
        query_params = [
            ("TARGET", DTX_CHARTS_V6_DS_ID, "equals", target),
            ("EXPERIMENT_ID", DTX_CHARTS_V6_DS_ID, "greaterthan", last_data_exp_id),
        ]
        query_data = construct_query_data(query_params)
        dn_ids = query_and_parse_dtx("10000", query_data)
        # Get registration info
        reg_data = get_parsed_reg_data_vw(dn_ids)
        # Get IC50 info
        assay_data = get_parsed_den_studies_charts_merged_v6(dn_ids)

        # Parse into tuples for easy dataframe generation
        for dn_id, dn_data in assay_data.items():
            dn_results = []
            for index, data in dn_data["data"].items():
                if (
                    (int(data.get("EXPERIMENT_ID")) > last_data_exp_id)
                    and (data.get("TARGET") == target)
                    and data.get("BASE64_GRAPH")
                ):
                    date = data["EXPERIMENT_DATE"]
                    date_str = datetime.datetime.strptime(date, "%m/%d/%Y %H:%M:%S")
                    stereo = (
                        reg_data.get(dn_id, {})
                        .get("data", {})
                        .get("STEREO_COMMENTS", "")
                    )
                    dn_results.append(
                        (
                            stereo,
                            target,
                            data["EXPERIMENT_ID"],
                            date_str,
                            data["ASSAY"],
                            data["ANALYSIS_NAME"],
                            data.get("CLASSIFICATION", ""),
                            data.get("RESULT", ""),
                            data["BASE64_GRAPH"],
                        )
                    )
                    if data["ANALYSIS_NAME"] not in analyses[data["ASSAY"]]:
                        analyses[data["ASSAY"]].append(data["ANALYSIS_NAME"])

            all_results[dn_id] = dn_results
        return all_results, analyses

    def get_agg_ic50_data(dn_ids):
        """
        Retrieves aggregate IC50 data for the given DN IDs
        :param dn_ids: a list of DN IDs
        :returns: a dictionary of the form {dn_id: aggregate_ic50_data}
        """
        return get_parsed_den_studies_ic50_agg_merged_v2(dn_ids)

    def get_logd_agg_data(date):
        """
        Wrapper to get the logD avg of DNs that have results since the given
        date from the LogD Agg table in Dotmatics
        :param date: date to query for more recent data in the form YYYYMMDD
        :returns: a list of mol objects with relevant properties
        """
        # Get new data since the given date
        query_params = [("EXP_CREATED_DATE_JMAX", "591", "greaterthan", date)]
        query_data = construct_query_data(query_params)

        mol_results = []
        dn_ids = query_and_parse_dtx("10000", query_data)
        if dn_ids:
            logd_data_dict = get_parsed_denali_logd_agg(dn_ids)
            structure_dict = get_parsed_structure(dn_ids)

            for dn, moltext in structure_dict.items():
                data = logd_data_dict[dn].get("data").get("1", "")
                if data:
                    # For some reason DTX has 4 compounds that don't follow date convention
                    if len(data.get("EXP_CREATED_DATE_JMAX", "")) > 8:
                        continue
                    try:
                        mol = Chem.MolFromMolBlock(moltext)
                        mol.SetProp("_Name", dn)
                        mol.SetProp("dn_id", dn)
                        mol.SetProp("logd_avg", data.get("LOG_D_JMEAN"))
                        mol.SetProp(
                            "exp_most_recent_date", data.get("EXP_CREATED_DATE_JMAX")
                        )
                        mol.SetProp("project_code", data.get("PROJECT_NAME"))
                        mol_results.append(mol)
                    except:
                        continue

        return mol_results

    def get_registered_structures(dn_ids):
        """
        Get moltext for the given DN IDs
        :param dn_ids: a list of DN IDs
        :returns: a dictionary of the form {dn_id: string of moltext}
        """
        return get_parsed_structure(dn_ids)

    def upload_file_to_dtx(file, exp_filename, exp_id, script_id=None):
        """
        Uploads a file to a Dotmatics experiment, optionally running a processing script upon upload.
        :param file: a file pointer to an open file w/ read permissions
        :param exp_filename: the name of the file as it should appear in the Dotmatics experiment
        :param exp_id: an int, the ID of a Dotmatics experiment to upload the file to
        :param script_id: (optional) an int, the ID of a Dotmatics script to run upon upload
        """
        return upload_file_to_dtx_experiment(file, exp_filename, exp_id, script_id)

    def get_all_dtx_properties():
        """
        Retrieve all assay properties for all compounds in Dotmatics
        :returns: a dictionary of the form {dn_id: {prop_name: prop_value}}
        """
        return get_parsed_denali_compound_assay_props([], all=True, first_only=True)

    def get_all_dn_after(dn_id=None, date=None):
        """
        Retrieve all DN IDs that were created after the given DN ID or date
        :param dn_id: a string, the DN ID to check against
        :param date: a datetime object, the date to check against
        :return: a list of DN IDs
        """
        if not dn_id and not date:
            raise ValueError("Must provide either a DN ID or a date")
        elif dn_id and date:
            raise ValueError("Must provide either a DN ID or a date, not both")
        elif dn_id:
            return query_and_parse_reg_data_vw_by_dn_id(dn_id)
        elif date:
            return query_and_parse_reg_data_vw_by_reg_date(date)

except ImportError:
    #################################
    #             STUBS             #
    #################################
    def check_dtx_for_inchi(inchi):
        """
        Query Dotmatics to determine if a compound already exists in registry
        :param inchi: the inchi of the molecule (generated with the '/suu' flag)
        :returns: a string, the DN ID of the matching compound
        """
        return ""

    def get_ic50_data(last_data_exp_id, target):
        """
        Wrapper to get the assay data of DNs for the given target that have new results
        :param last_data_exp_id: a 6 digit ELN id of the last sent data
        :param target: a Project's target name (as used in Dotmatics)
        :returns: a tuple (results, analyses) where
            - [0] results is a dictionary of tuples base64 encoded ic50 curves that looks like:
            {dn_id: [(stereochemistry, target, experiment id, experiment datetime object, assay name, ic50 value, base64 encoded ic50 curve)]}
            - [1] analyses is a dictionary of the form {assay_name: [analysis_names]} containing all the
            analyses whose data are in results
        """
        return {}, {}

    def get_agg_ic50_data(dn_ids):
        """
        Retrieves aggregate IC50 data for the given DN IDs
        :param dn_ids: a list of DN IDs
        :returns: a dictionary of the form {dn_id: aggregate_ic50_data}
        """
        return {dn_id: {} for dn_id in dn_ids}

    def get_logd_agg_data(date):
        """
        Wrapper to get the logD avg of DNs that have results since the given
        date from the LogD Agg table in Dotmatics
        :param date: date to query for more recent data in the form YYYYMMDD
        :returns: a list of mol objects with relevant properties
        """
        return []

    def get_registered_structures(dn_ids):
        """
        Get moltext for the given DN IDs
        :param dn_ids: a list of DN IDs
        :returns: a dictionary of the form {dn_id: string of moltext}
        """
        # The following string is what Dotmatics Register returns when the given DN ID does not exist
        no_structure_found = "  \r\n-ISIS-  03310609052D \r\n\r\n  0  0  0  0  0  0  0  0  0  0999 V2000 \r\nM  END\r\n\n"
        return {dn_id: no_structure_found for dn_id in dn_ids}

    def get_all_dtx_properties():
        """
        Retrieve all assay properties for all compounds in Dotmatics
        :returns: a dictionary of the form {dn_id: {prop_name: prop_value}}
        """
        return {}

    def get_all_dn_after(dn_id=None, date=None):
        """
        Retrieve all DN IDs that were created after the given DN ID or date
        :param dn_id: a string, the DN ID to check against
        :param date: a datetime object, the date to check against
        :return: a list of DN IDs
        """
        if not dn_id and not date:
            raise ValueError("Must provide either a DN ID or a date")
        elif dn_id and date:
            raise ValueError("Must provide either a DN ID or a date, not both")
        return []

    def upload_file_to_dtx(file, exp_filename, exp_id, script_id=None):
        """
        Uploads a file to a Dotmatics experiment, optionally running a processing script upon upload.
        :param file: a file pointer to an open file w/ read permissions
        :param exp_filename: the name of the file as it should appear in the Dotmatics experiment
        :param exp_id: an int, the ID of a Dotmatics experiment to upload the file to
        :param script_id: (optional) an int, the ID of a Dotmatics script to run upon upload
        """
        return
