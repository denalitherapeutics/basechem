from rdkit import Chem


def mock_run_inductive_lm_predict(input_file, species, image=False):
    """
    Mock for `run_inductive_lm_predict` to use in all test cases except those in `test_inductive_utils.py`
    """
    suppl = Chem.SDMolSupplier(input_file)
    return_data = []
    for mol in suppl:
        compound_data = {}
        if species == "R":
            compound_data = {
                "name": mol.GetProp("_Name"),
                "prediction": "54",
                "measured": "",
                "out_of_domain": "True",
                "latest_data_date": "2022-11-25",
                "model_version": "2.0.0",
                "interp_image": "",
                "probs_image": "",
                "probs_list": [0.049, 0.051, 0.9],
            }
        elif species == "H":
            compound_data = {
                "name": mol.GetProp("_Name"),
                "prediction": "20",
                "measured": "",
                "out_of_domain": "True",
                "latest_data_date": "2022-11-25",
                "model_version": "2.0.0",
                "interp_image": "",
                "probs_image": "",
                "probs_list": [0.2, 0.4, 0.4],
            }
        if image:
            compound_data[
                "interp_image"
            ] = "iVBORw0KGgoAAAANSUhEUgAAAAgAAAAIAQMAAAD+wSzIAAAABlBMVEX///+/v7+jQ3Y5AAAADklEQVQI12P4AIX8EAgALgAD/aNpbtEAAAAASUVORK5CYII=="
            compound_data[
                "probs_image"
            ] = "iVBORw0KGgoAAAANSUhEUgAAAAgAAAAIAQMAAAD+wSzIAAAABlBMVEX///+/v7+jQ3Y5AAAADklEQVQI12P4AIX8EAgALgAD/aNpbtEAAAAASUVORK5CYII=="
        return_data.append(compound_data)
    return return_data


def mock_run_inductive_alogd_predict(input_file):
    """
    Mock for `run_inductive_alogd_predict` to use in all test cases except those in `test_inductive_utils.py`
    """
    suppl = Chem.SDMolSupplier(input_file)
    return_data = []
    for mol in suppl:
        compound_data = {
            "name": mol.GetProp("_Name"),
            "prediction": "2.30",
            "latest_data_date": "2022-11-25",
            "model_version": "2.0.0",
            "measured": "2.22",
        }

        return_data.append(compound_data)

    return return_data


def mock_update_inductive_logd_data(input_file):
    """
    Mock for `update_inductive_logd_data` to use in all test cases
    Does nothing since we don't want to PUT any data
    """
    return
