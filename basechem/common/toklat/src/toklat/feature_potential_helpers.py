import numpy as np
from scipy.stats import norm

# This module contains helper functions for calculating interatomic potentials.
# Most existing functions are based on those from AutoDock Vina
# (https://doi.org/10.1002/jcc.21334). Each potential function provided to an
# instance of PotentialTransformer will be called to add one or more potential
# key-value pairs to the list of potentials for each atom pair.

# New potential functions can be added, and must be callable with the following
# signature:
# Arguments:
# - dist: 1D numpy array of atom surface distances
# - angle: 1D numpy array of atom angles
# - lt: 1D list of ligand atom types
# - pt: 1D list of protein atom types
# Returns:
# - 1D list of dictionaries, where each dictionary contains potential names as keys
#   and potential values as values. If a potential is not present (e.g. two atoms have
#   a repulsive potential of 0), the dictionary can be empty.
# Unneeded arguments can be ignored, but must be present in the function signature.


def _constrain(vals):
    # helper function to constrain values between 0 and 1
    return np.clip(vals, 0, 1)


class HBondPotential:
    def __init__(
        self,
        hbond_dist_lower: float,
        hbond_dist_upper: float,
        hbond_angle_lower: float,
        hbond_angle_upper: float,
        name_base: str,
    ):
        """
        Hydrogen bond potential function.

        Based on AutoDock Vina, with angle component based on GOLD PLP.
        Looks for "_hba" and "_hbd" suffixes in atom types to determine whether a
        hydrogen bond potential should be calculated.

        Calling this function will return two potentials, with keys {name_base}_dist
        and {name_base}_dist_angle.
        So, assuming name_base is "hbond", the output will look something like:
        {"hbond_dist": 0.5, "hbond_dist_angle": 0.3}
        The first potential is based only on distance, matching the AutoDock Vina
        definition, while the second potential is based on both distance and angle.

        :param hbond_dist_lower: Lower bound for maximum hydrogen bond potential
        :param hbond_dist_upper: Upper bound for 0 hydrogen bond potential
        :param hbond_angle_lower: Lower bound for 0 hydrogen bond angle potential
        :param hbond_angle_upper: Upper bound for maximum hydrogen bond DHA angle
            potential. (Note, 180 indicates "perfect" DHA alignment).
        :param name_base: Base name for potential output
        """
        self.hbond_slope = 1 / (hbond_dist_lower - hbond_dist_upper)
        self.hbond_intercept = -hbond_dist_upper / (hbond_dist_lower - hbond_dist_upper)
        self.hbond_angle_slope = 1 / (hbond_angle_upper - hbond_angle_lower)
        self.hbond_angle_intercept = -hbond_angle_lower / (
            hbond_angle_upper - hbond_angle_lower
        )
        self.name_base = name_base

    def __call__(self, dist, angle, lt, pt):
        is_hbond = 1 * np.array(
            [
                (("_hba" in lig and "_hbd" in prt) or ("_hbd" in lig and "_hba" in prt))
                for lig, prt in zip(lt, pt)
            ]
        )
        hbond_dist = (
            _constrain(self.hbond_slope * dist + self.hbond_intercept) * is_hbond
        )
        hbond_angle = hbond_dist * _constrain(
            self.hbond_angle_slope * angle + self.hbond_angle_intercept
        )
        out = [{} for _ in range(len(dist))]
        for i, (d, a) in enumerate(zip(hbond_dist, hbond_angle)):
            if d > 0:
                out[i][self.name_base + "_dist"] = float(d)
            if a > 0:
                out[i][self.name_base + "_dist_angle"] = float(a)
        return out


class HydrophobicPotential:
    def __init__(
        self,
        hydrophobic_dist_lower: float,
        hydrophobic_dist_upper: float,
        name_base: str,
    ):
        """
        Hydrophobic potential function.

        Based on AutoDock Vina, with a linear potential based on distance.
        Looks for ligand and protein atoms with the "_hydrophobic" suffix in the
        atom type to determine whether a hydrophobic potential should be calculated.

        :param hydrophobic_dist_lower: Lower bound for maximum hydrophobic potential
        :param hydrophobic_dist_upper: Upper bound for 0 hydrophobic potential
        :param name_base: Base name for potential output
        """
        self.hydrophobic_slope = 1 / (hydrophobic_dist_lower - hydrophobic_dist_upper)
        self.hydrophobic_intercept = -hydrophobic_dist_upper / (
            hydrophobic_dist_lower - hydrophobic_dist_upper
        )
        self.name_base = name_base

    def __call__(self, dist, angle, lt, pt):
        is_hydrophobic = 1 * np.array(
            [
                lig.split("_")[-1] == "hydrophobic"
                and prt.split("_")[-1] == "hydrophobic"
                for lig, prt in zip(lt, pt)
            ]
        )
        hydrophobic = (
            _constrain(self.hydrophobic_slope * dist + self.hydrophobic_intercept)
            * is_hydrophobic
        )
        return [{} if h <= 0 else {self.name_base: float(h)} for h in hydrophobic]


class VDWPotential:
    def __init__(self, vdw_feat_dist: float, vdw_feat_scale: float, name_base: str):
        """
        Van der Waals potential function.

        Based on AutoDock Vina, with a Gaussian potential based on distance.

        :param vdw_feat_dist: Distance for maximum value of the Gaussian potential
        :param vdw_feat_scale: Standard deviation for the Gaussian potential
        """
        self.vdw_feat_dist = vdw_feat_dist
        self.vdw_feat_scale = vdw_feat_scale
        self.name_base = name_base

    def __call__(self, dist, angle, lt, pt):
        vdw = norm.pdf(dist, loc=self.vdw_feat_dist, scale=self.vdw_feat_scale)
        return [{} if v <= 0 else {self.name_base: float(v)} for v in vdw]


class RepulsivePotential:
    def __init__(self, dist_start, name_base: str):
        """
        Repulsive potential function.

        Based on AutoDock Vina, with a quadratic potential based on degree
        of overlap between atoms.

        :param dist_start: Distance at which the repulsive potential starts to
            be applied
        :param name_base: Base name for potential output
        """
        self.dist_start = dist_start
        self.name_base = name_base

    def __call__(self, dist, angle, lt, pt):
        dist = dist - self.dist_start
        repulsion = np.minimum(dist, 0) ** 2
        return [{} if r <= 0 else {self.name_base: float(r)} for r in repulsion]


class PairwisePotential:
    def __init__(
        self,
        pairwise_dist_lower: float,
        pairwise_dist_upper: float,
        pairwise_max_val: float,
        name_base: str,
    ):
        """
        Pairwise potential function.

        This function calculates a distinctly named potential for each pair of atom
        types, with a linear potential based on distance.

        The potential output will have a key composed of the two atom types separated
        by a pipe character, e.g. "C|O", and preceded by any provided name_base.

        :param pairwise_dist_lower: Lower bound for maximum pairwise potential
        :param pairwise_dist_upper: Upper bound for 0 pairwise potential
        :param pairwise_max_val: Maximum value for the pairwise potential (achieved at
            pairwise_dist_lower)
        :param name_base: Base name for potential output
        """
        self.pairwise_slope = 1 / (pairwise_dist_lower - pairwise_dist_upper)
        self.pairwise_intercept = -pairwise_dist_upper / (
            pairwise_dist_lower - pairwise_dist_upper
        )
        self.pairwise_max_val = pairwise_max_val
        self.name_base = name_base

    def __call__(self, dist, angle, lt, pt):
        pairwise = self.pairwise_max_val * _constrain(
            self.pairwise_slope * dist + self.pairwise_intercept
        )
        return [
            {} if v <= 0 else {f"{self.name_base}{l}|{p}": float(v)}
            for v, l, p in zip(pairwise, lt, pt)
        ]
