import os

import pandas as pd
from rdkit import Chem

from basechem.common.mmpdb_utils import generate_mmpdb
from basechem.common.rdkit_utils import RDKitWrappers
from basechem.main.constants import (
    MMP_MAX_H_BOND_DONORS,
    MMP_MAX_LOGP,
    MMP_MAX_MW,
    MMP_MAX_TPSA,
)


def meets_mmp_prop_req(mol):
    """
    Checks if the compound is within the property criteria for a potential MMP
    :param mol: the mol of the compound
    :return: True if the compound is within property criteria
    """
    RDKitWrappers.generate_properties(mol)
    mw = float(mol.GetProp("mw"))
    logp = float(mol.GetProp("clogp"))
    tpsa = float(mol.GetProp("tpsa"))
    donors = float(mol.GetProp("donors"))
    return (
        mw < MMP_MAX_MW
        and logp < MMP_MAX_LOGP
        and tpsa < MMP_MAX_TPSA
        and donors < MMP_MAX_H_BOND_DONORS
    )


def generate_mmps(input_smiles, generate_path, radius=0):
    """
    Returns a dataframe with mmpdb generate output filtered to only include results that meet property criteria
    :param input_smiles: the smiles to generate mmps for
    :param generate_path: the path to the file where mmpdb generate output will be saved
    :param radius: the radius argument value for generate
    :return: a dataframe of generate output filtered
    """
    if not os.path.exists(generate_path):
        generated = generate_mmpdb(input_smiles, generate_path, radius)
        if not generated:
            return pd.DataFrame()
    generate_df = pd.read_csv(generate_path, sep="\s+")
    generate_df["final_mol"] = generate_df["final"].apply(
        lambda smiles: Chem.MolFromSmiles(smiles)
    )
    generate_df["meets_mmp_prop_req"] = generate_df["final_mol"].apply(
        lambda mol: meets_mmp_prop_req(mol)
    )

    return generate_df
