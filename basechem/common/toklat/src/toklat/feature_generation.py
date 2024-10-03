from collections import OrderedDict
from typing import Callable, Dict, List, Sequence

import numpy as np
import prolif as plf
from prolif.interactions.base import SingleAngle
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors
from scipy.spatial.distance import cdist
from sklearn.base import BaseEstimator, TransformerMixin

# _ATOM_TYPING_SMARTS and _DEFAULT_FP below can be modified
# to change the atom typing for distance features, and the
# multi-atom features extracted by prolif, respectively.
# When modifying _DEFAULT_FP, new interactions can be added
# but proceed with caution before modifying HBAcceptor or HBDonor
# features, as these are treated specially and used to update
# the distance features with hydrogen bond angles.

# atom typing smarts taken from prolif
_ATOM_TYPING_SMARTS = OrderedDict(
    {
        "arom": Chem.MolFromSmarts("[c,n,o,s]"),
        "hba": Chem.MolFromSmarts(
            (
                "[#7&!$([nX3])&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4&+1]),"
                "O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O),o+0,"
                "F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]"
            )
        ),
        "hbd": Chem.MolFromSmarts("[$([O,S;+0]),$([N;v3,v4&+1]),n+0]-[H]"),
        "hydrophobic": Chem.MolFromSmarts(
            ("[c,s,Br,I,S&H0&v2,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6X4H0]);+0]")
        ),
        "cation": Chem.MolFromSmarts("[+{1-},$([NX3&!$([NX3]-O)]-[C]=[NX3+])]"),
        "anion": Chem.MolFromSmarts("[-{1-},$(O=[C,S,P]-[O-])]"),
    }
)


class HBAcceptorRotatable(SingleAngle):
    """
    Custom prolif interaction for H-bonds involving OH, NH3, and SH groups.
    """

    def __init__(
        self,
        # standard prolif acceptor smarts
        acceptor=(
            "[#7&!$([nX3])&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4&+1]),"
            "O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O),o+0,"
            "F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]"
        ),
        # only terminal sp3 hybridized NH2/3, OH, SH donors
        donor="[!#1]-[$([N^3H2,N^3H3&+1]),$([OH1]),$([SH1])]",
        distance=3.5,
        DHA_angle=(0, 180),
    ):
        super().__init__(
            lig_pattern=acceptor,
            prot_pattern=donor,
            distance=distance,
            angle=DHA_angle,
            distance_atom="P2",
            metadata_mapping={"angle": "RDA_angle"},
        )


HBDonorRotatable = HBAcceptorRotatable.invert_role("HBDonorRotatable", "")

_DEFAULT_FP = plf.Fingerprint(
    interactions=[
        "HBAcceptor",
        "HBDonor",
        "HBAcceptorRotatable",
        "HBDonorRotatable",
        "PiStacking",
        "CationPi",
        "PiCation",
    ],
    parameters={
        "HBAcceptor": {
            "distance": 6,
            "DHA_angle": (0, 180),
            "donor": "[$([N;v3,v4&+1])&!$([N^3H2,N^3H3&+1]),n+0]-[H]",
        },
        "HBDonor": {
            "distance": 6,
            "DHA_angle": (0, 180),
            "donor": "[$([N;v3,v4&+1])&!$([N^3H2,N^3H3&+1]),n+0]-[H]",
        },
        "HBAcceptorRotatable": {
            "distance": 6,
        },
        "HBDonorRotatable": {
            "distance": 6,
        },
        "PiCation": {
            "distance": 6,
        },
        "CationPi": {
            "distance": 6,
        },
    },
    count=True,
    vicinity_cutoff=6,
)

# The following are examples
# of additional_ligand_features functions.
# any function can be used that takes an RDKit molecule
# as input and returns a float.


class NumHeavyAtomsGetter:
    def __call__(self, mol):
        return rdMolDescriptors.CalcNumHeavyAtoms(mol)


class CLogPGetter:
    def __call__(self, mol):
        return Crippen.MolLogP(mol)


class TPSAGetter:
    def __call__(self, mol):
        return rdMolDescriptors.CalcTPSA(mol)


class NumRotatableBondsGetter:
    def __call__(self, mol):
        return rdMolDescriptors.CalcNumRotatableBonds(mol)


class PropAsFloatGetter:
    def __init__(self, prop_name: str, default: float = 0):
        """
        Helper class to get a property from an RDKit molecule as a float.

        Can be passed as a additional ligand feature function to
        DistanceAngleTypingTransformer.
        :param prop_name: The name of the property to get.
        :param default: The default value to return if the property is not present
            or cannot be converted to a float.
        """
        self.prop_name = prop_name
        self.default = default

    def __call__(self, mol):
        try:
            return float(mol.GetProp(self.prop_name))
        except (KeyError, ValueError):
            print(
                f"WARNING: {self.prop_name} not found in molecule.",
                f"Setting to {self.default}.",
            )
            return self.default


def get_atom_types(mol):
    """
    Return a list of atom types for each atom in a molecule.

    Each atom's type is a string composed of the atom element and any additional
    matched _ATOM_TYPING_SMARTS patterns, separated by underscores.
    """
    types = [a.GetSymbol() for a in mol.GetAtoms()]
    for key, smarts in _ATOM_TYPING_SMARTS.items():
        has_type = [False] * len(types)
        for match in mol.GetSubstructMatches(smarts, maxMatches=100_000):
            for i in match:
                has_type[i] = True
        types = [t + "_" + key if has_type[i] else t for i, t in enumerate(types)]
    return types


def get_atom_types_protein(mol):
    """
    Hack to deal with bad typing of Ns in His residues.
    """
    types = get_atom_types(mol)
    resns = [
        atom.GetPDBResidueInfo().GetResidueName() if atom.GetPDBResidueInfo() else ""
        for atom in mol.GetAtoms()
    ]
    new_types = []
    for t, resn in zip(types, resns):
        if resn in ["HIS", "HID", "HIE"] and t in [
            "N_arom_cation",
            "N_cation",
            "N_arom_hba_cation",
            "N_arom_hbd_cation",
        ]:
            new_types.append("N_arom_hbd")
        elif resn in ["HIS", "HID", "HIE"] and t in [
            "N_arom_hba_anion",
            "N_anion",
            "N_arom",
            "N_hba_anion",
        ]:
            new_types.append("N_arom_hba")
        else:
            new_types.append(t)
    return new_types


class DistanceAngleTypingTransformer(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        receptor_path: str = None,
        additional_ligand_feature_functions: Dict[str, Callable] | None = None,
    ):
        """
        Extracts features from ligand-receptor complexes.

        :param receptor_path: The path to the receptor PDB file.
        :param additional_ligand_feature_functions: A dictionary of functions to
            extract additional ligand features from input ligands. The keys are
            the names of the features, and the values are functions that take an
            RDKit molecule as input and return a float.
        """
        self.receptor_path = receptor_path
        if receptor_path is not None:
            self.receptor = plf.Molecule(load_pdb_to_mol(receptor_path))
        else:
            self.receptor = None
        if additional_ligand_feature_functions is None:
            additional_ligand_feature_functions = {}
        self.additional_ligand_feature_functions = additional_ligand_feature_functions

    def fit(self, X: Sequence[Chem.Mol], y=None):
        return self

    def update_receptor(self, receptor_path: str):
        """
        Update the receptor used for feature extraction.

        :param receptor_path: The path to the new receptor PDB file.
        """
        self.receptor = plf.Molecule(load_pdb_to_mol(receptor_path))
        self.receptor_path = receptor_path
        return self

    def transform(self, X: Sequence[Chem.Mol]):
        """
        Given a sequence of ligands, extract features from each ligand-receptor complex.

        Returned dictionaries have the following keys:

        - 'distance_features': A list of dictionaries, where each dictionary contains
          information about a single ligand-receptor atom-atom interaction. It's keys
          are:
                - 'li': The index of the ligand atom.
                - 'pi': The index of the protein atom.
                - 'd': The distance between the atoms.
                - 'lt': The type of the ligand atom.
                - 'pt': The type of the protein atom.
                - 'ang': DHA angle, if the 2 atoms can form a hydrogen bond.
        - 'additional_prolif_features': A list of dictionaries, where each dictionary
          contains information about a single multi-atom interaction, such as pi
          stacking.
        - 'additional_ligand_features': A dictionary containing the additional ligand
            features for the ligand.

        :param X: A sequence of RDKit molecules representing ligands.
        :return: A list of dictionaries, where each dictionary contains the features
            for a single ligand-receptor complex.
        """
        if self.receptor is None:
            raise ValueError("Receptor must be set before transforming ligands.")
        additional_ligand_features = []
        for mol in X:
            additional_ligand_features.append(
                {k: f(mol) for k, f in self.additional_ligand_feature_functions.items()}
            )
        distance_features = calc_distance_interactions(
            X, self.receptor, get_atom_types_protein, get_atom_types
        )
        # get hbond angles and multi-atom interactions from prolif
        prolif_features = calc_prolif_interactions(X, self.receptor, _DEFAULT_FP)
        out_feats = []
        for d_feats, p_feats, g_feats in zip(
            distance_features, prolif_features, additional_ligand_features
        ):
            # add hbond angles to distance features
            hb_angles = p_feats["HB_angles"]
            for d_feat in d_feats:
                d_feat["hb_ang"] = hb_angles.get((d_feat["li"], d_feat["pi"]), None)

            out_feats.append(
                {
                    "distance_features": d_feats,
                    "additional_prolif_features": p_feats[
                        "additional_prolif_interactions"
                    ],
                    "additional_ligand_features": g_feats,
                }
            )
        return out_feats


def fix_anions_inplace(mol):
    """Rdkit treats anion from pdb files as a neutral with an implicit H."""
    for atom in mol.GetAtoms():
        if atom.GetNumImplicitHs() == 1 and atom.GetNumExplicitHs() == 0:
            atom.SetFormalCharge(-1)
    mol.UpdatePropertyCache(strict=False)


def cleanup_valence_violations(mol):
    rwmol = Chem.RWMol(mol)
    atoms_with_violations = [a for a in rwmol.GetAtoms() if a.HasValenceViolation()]
    for a in atoms_with_violations:
        bonds = a.GetBonds()
        for b in bonds:
            a2 = b.GetOtherAtom(a)
            if (
                a.GetPDBResidueInfo().GetResidueNumber()
                != a2.GetPDBResidueInfo().GetResidueNumber()
            ):
                rwmol.RemoveBond(a.GetIdx(), a2.GetIdx())
    return rwmol.GetMol()


def load_pdb_to_mol(pdb_path):
    """Load a PDB file into an RDKit molecule, fixing issues with anions."""
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
    mol = cleanup_valence_violations(mol)
    Chem.SanitizeMol(mol)
    fix_anions_inplace(mol)
    return mol


def calc_distance_interactions(
    ligand_mols: list,
    receptor_mol: Chem.Mol,
    protein_typing_fn: callable,
    ligand_typing_fn: callable,
    distance_cutoff: float = 6.0,
) -> List[dict]:
    """
    Calculates distance-based interactions between ligand and receptor atoms.

    Used by DistanceAngleTypingTransformer.
    """
    prot_types = protein_typing_fn(receptor_mol)
    prot_pos = receptor_mol.GetConformer().GetPositions()
    prot_elements = np.array([atom.GetAtomicNum() for atom in receptor_mol.GetAtoms()])
    ligand_features_list = []
    for lig in ligand_mols:
        lig_types = ligand_typing_fn(lig)
        lig_pos = lig.GetConformer().GetPositions()
        lig_elements = np.array([atom.GetAtomicNum() for atom in lig.GetAtoms()])
        lig_prot_distance = cdist(lig_pos, prot_pos)
        mask = lig_prot_distance < distance_cutoff
        edge_dist = lig_prot_distance[mask].tolist()
        i_indices, j_indices = np.where(mask)
        atom_interactions = []
        for i, j, d in zip(i_indices, j_indices, edge_dist):
            # don't include interactions with H
            if lig_elements[i] == 1 or prot_elements[j] == 1:
                continue
            atom_interactions.append(
                {
                    "li": int(i),
                    "pi": int(j),
                    "d": float(d),
                    "lt": lig_types[i],
                    "pt": prot_types[j],
                }
            )
        ligand_features_list.append(atom_interactions)
    return ligand_features_list


def _flatten_interactions(interactions):
    """
    Helper function to flatten prolif interactions into a list of dictionaries
    """
    flat_interactions = []
    for resid_interactions in interactions.values():
        for i_type, i_list in resid_interactions.items():
            for i in i_list:
                i = i.copy()
                i["ligand_idxs"] = i["parent_indices"]["ligand"]
                i["protein_idxs"] = i["parent_indices"]["protein"]
                del i["indices"]
                del i["parent_indices"]
                i["interaction_type"] = i_type
                flat_interactions.append(i)
    return flat_interactions


def calc_prolif_interactions(
    ligand_mols: list,
    receptor_mol: plf.Molecule,
    fp: plf.Fingerprint = _DEFAULT_FP,
) -> List[dict]:
    """
    Calculates H-bond angles and multi-atom interactions.

    Used by DistanceAngleTypingTransformer.
    """
    ligs = [plf.Molecule.from_rdkit(lig) for lig in ligand_mols]
    fp.run_from_iterable(ligs, receptor_mol, progress=False)
    all_plf_interactions = [_flatten_interactions(fp.ifp[i]) for i in range(len(ligs))]
    ligand_features_list = []
    for plf_interactions in all_plf_interactions:
        HB_angles = {}
        multi_atom_interactions = []
        for interact in plf_interactions:
            if interact["interaction_type"] in ["HBAcceptor", "HBDonor"]:
                HB_angles[
                    (interact["ligand_idxs"][0], interact["protein_idxs"][0])
                ] = interact["DHA_angle"]
            elif interact["interaction_type"] in [
                "HBAcceptorRotatable",
                "HBDonorRotatable",
            ]:
                # need to grab last element instead of first
                interact["DHA_angle"] = 180 - abs(109.5 - interact["RDA_angle"])
                HB_angles[
                    (interact["ligand_idxs"][-1], interact["protein_idxs"][-1])
                ] = interact["DHA_angle"]
            else:
                multi_atom_interactions.append(interact)
        ligand_features_list.append(
            {
                "HB_angles": HB_angles,
                "additional_prolif_interactions": multi_atom_interactions,
            }
        )
    return ligand_features_list
