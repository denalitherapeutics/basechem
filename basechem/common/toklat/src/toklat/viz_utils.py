import numpy as np
import py3Dmol
from rdkit import Chem

from .feature_generation import load_pdb_to_mol


def visualize_interpretation(
    mol: Chem.Mol,
    receptor_path: str,
    top_interactions: list,
    unsatisfied_ligand_atoms: list,
    unsatisfied_protein_atoms: list,
    interaction_scale: float = 0.5,
):
    """
    Create a py3Dmol visualization of a scored pose.

    :param mol: rdkit molecule
    :param receptor_path: path to the receptor pdb file
    :param top_interactions: list of dicts with keys "li" (ligand index)
        "pi" (protein index) and "summed_interaction", output
        by ScoringFunctionModel.predict_and_interpret
    :param unsatisfied_ligand_atoms: list of dicts with key "li",
        output by ScoringFunctionModel.predict_and_interpret
    :param unsatisfied_protein_atoms: list of dicts with key "pi",
        output by ScoringFunctionModel.predict_and_interpret
    :param interaction_scale: scaling factor for interaction line width.
    :return: py3Dmol.view object
    """
    mol = Chem.Mol(mol)
    pdb_mol = load_pdb_to_mol(receptor_path)
    p = py3Dmol.view(width=400, height=400)
    # add red and blue lines for interactions
    for contrib_dict in top_interactions:
        loc1 = mol.GetConformer().GetAtomPosition(contrib_dict["li"])
        loc2 = pdb_mol.GetConformer().GetAtomPosition(contrib_dict["pi"])
        val = contrib_dict["summed_interaction"]
        p.addCylinder(
            {
                "start": dict(x=loc1.x, y=loc1.y, z=loc1.z),
                "end": dict(x=loc2.x, y=loc2.y, z=loc2.z),
                "color": "blue" if val < 0 else "red",
                "radius": 0.15
                * np.power(np.abs(val), 0.5)
                / np.power(interaction_scale, 0.5),
                "dashed": True,
                "fromCap": 1,
                "toCap": 1,
            }
        )
    # add gray spheres for unsatisfied atoms
    for atom_dict in unsatisfied_ligand_atoms:
        loc = mol.GetConformer().GetAtomPosition(atom_dict["li"])
        p.addSphere(
            {
                "center": dict(x=loc.x, y=loc.y, z=loc.z),
                "color": "gray",
                "alpha": 0.5,
                "radius": 0.5,
                "wireframe": True,
            }
        )
    for atom_dict in unsatisfied_protein_atoms:
        loc = pdb_mol.GetConformer().GetAtomPosition(atom_dict["pi"])
        p.addSphere(
            {
                "center": dict(x=loc.x, y=loc.y, z=loc.z),
                "color": "gray",
                "alpha": 0.5,
                "radius": 0.5,
                "wireframe": True,
            }
        )
    p.addModel(Chem.MolToMolBlock(Chem.RemoveHs(mol)), "sdf")
    p.setStyle({"model": 0}, {"stick": {"colorscheme": "pinkCarbon"}})
    p.zoomTo({"model": 0})
    p.addModel(open(receptor_path).read(), "pdb")
    p.setStyle(
        {"model": 1},
        {
            "stick": {"colorscheme": "lightGrayCarbon"},
        },
    )
    return p
