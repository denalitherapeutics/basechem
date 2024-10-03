import json
import multiprocessing as mp
from copy import deepcopy
from typing import Sequence

import numpy as np
from scipy.sparse import lil_array
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

from .feature_potential_helpers import (
    HBondPotential,
    HydrophobicPotential,
    PairwisePotential,
    RepulsivePotential,
    VDWPotential,
)

_DEFAULT_POTENTIAL_FUNCTIONS = [
    VDWPotential(0, 0.5, "vdw"),
    HydrophobicPotential(0.5, 1.5, "hydrophobic"),
    HBondPotential(-0.7, 0, 100, 150, "hbond"),
    PairwisePotential(0, 1.5, 1.0, ""),
    RepulsivePotential(0, "repulsive"),
]


_VDW_RADII = {
    "C": 2.0,  # C radius accounts for implicit Hs
    "C_arom": 1.8,
    "N": 1.75,
    "O": 1.6,
    "S": 2.0,
    "P": 2.1,
    "F": 1.545,
    "Cl": 2.045,
    "Br": 2.165,
    "I": 2.36,
    "B": 2,  # Note: not in smina atom types, chosen somewhat arbitrarily
}


def _vdw_type(type_str):
    """
    Convert full atom type to type for VDW radius lookup.
    """
    type_parts = type_str.split("_")
    if type_parts[0] != "C":
        return type_parts[0]
    elif len(type_parts) > 1 and type_parts[1] == "arom":
        return "C_arom"
    else:
        return "C"


class PotentialVectorizer(BaseEstimator, TransformerMixin):
    def __init__(self, sort: bool = True):
        """
        Vectorize interatomic potentials and other pose features into a sparse array.

        :param sort: if True, sort feature names in alphabetical order
        """
        self.sort = sort

    def fit(self, X: Sequence[str] | Sequence[dict], y=None):
        vocab = {}
        feature_names = []
        for x in X:
            if type(x) == str:
                with open(x, "r") as f:
                    x = json.load(f)
            for feat_dict in x["distance_features"]:
                for k in feat_dict["potentials"].keys():
                    if k not in vocab:
                        vocab[k] = len(vocab)
                        feature_names.append(k)
            for k in x["additional_ligand_features"].keys():
                if k not in vocab:
                    vocab[k] = len(vocab)
                    feature_names.append(k)
            for feat_dict in x["additional_prolif_features"]:
                for k in feat_dict["potentials"].keys():
                    if k not in vocab:
                        vocab[k] = len(vocab)
                        feature_names.append(k)
        if self.sort:
            feature_names.sort()
            vocab = {f: i for i, f in enumerate(feature_names)}

        self.feature_names_ = feature_names
        self.vocabulary_ = vocab
        return self

    def transform(self, X):
        """
        Transform a list of pose feature dicts into a sparse matrix.
        """
        check_is_fitted(self, ["feature_names_", "vocabulary_"])
        out_matrix = lil_array((len(X), len(self.vocabulary_)))
        for i, x in enumerate(X):
            if type(x) == str:
                with open(x, "r") as f:
                    x = json.load(f)
            self._transform_instance(out_matrix, x, i)
        return out_matrix

    def _transform_instance(self, feat_matrix, x, idx):
        for feat_dict in x["distance_features"]:
            for k, v in feat_dict["potentials"].items():
                if k in self.vocabulary_:
                    feat_matrix[idx, self.vocabulary_[k]] += v
        for k, v in x["additional_ligand_features"].items():
            if k in self.vocabulary_:
                feat_matrix[idx, self.vocabulary_[k]] += v
        for feat_dict in x["additional_prolif_features"]:
            for k, v in feat_dict["potentials"].items():
                if k in self.vocabulary_:
                    feat_matrix[idx, self.vocabulary_[k]] += v


class PotentialTransformer(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        vdw_radii=None,
        potential_functions=None,
    ):
        """
        Calculates interatomic potentials based on distance features.

        :param vdw_radii: dictionary of atomic radii for potentials. Defaults to
            _VDW_RADII
        """
        if vdw_radii is None:
            vdw_radii = _VDW_RADII
        if potential_functions is None:
            potential_functions = _DEFAULT_POTENTIAL_FUNCTIONS
        self.vdw_radii = vdw_radii
        self.potential_functions = potential_functions

    def fit(self, X, y=None):
        return self

    def _transform_instance_on_disk(self, x):
        with open(x, "r") as f:
            x_dict = json.load(f)
        x_dict = self._transform_instance(x_dict)
        with open(x, "w") as f:
            json.dump(x_dict, f)

    def transform(self, X: Sequence[str] | Sequence[dict], n_cores=-1):
        """
        Calculate interatomic potentials.

        This transform function can either be passed input values in memory or
        can be used to modify files on disk in place if passed a list
        of file paths.

        Interatomic potentials will be added to each pairwise distance
        feature as a dictionary mapping potential names to potential values.

        :param X: list of file paths or dictionaries of features. If file paths, the
            files should be json files with the same format as the output of
            DistanceAngleTypingTransformer, and will be modified in place.
        :param n_cores: number of cores to use for parallel processing of files. If -1,
            will use all available cores.
        """
        if type(X[0]) == str:
            if n_cores == -1:
                n_cores = mp.cpu_count()
            with mp.Pool(n_cores) as pool:
                pool.map(self._transform_instance_on_disk, X)
            return X
        else:
            X = deepcopy(X)
            new_X = []
            for x in X:
                x = self._transform_instance(x)
                new_X.append(x)
            return new_X

    def _transform_instance(self, x):
        feats = x["distance_features"]
        len_feats = len(feats)
        potentials = [{} for _ in range(len_feats)]
        vdw_1 = np.array(
            [self.vdw_radii.get(_vdw_type(feat_dict["lt"]), 1.0) for feat_dict in feats]
        )
        vdw_2 = np.array(
            [self.vdw_radii.get(_vdw_type(feat_dict["pt"]), 1.0) for feat_dict in feats]
        )
        angle = np.array([feat_dict.get("hb_ang") or 0 for feat_dict in feats])
        d = np.array([feat_dict["d"] for feat_dict in feats])
        surface_dist = d - vdw_1 - vdw_2
        lt = np.array([feat_dict["lt"] for feat_dict in feats])
        pt = np.array([feat_dict["pt"] for feat_dict in feats])

        for potential_func in self.potential_functions:
            results = potential_func(surface_dist, angle, lt, pt)
            for i, result in enumerate(results):
                potentials[i].update(result)
        for i, feat_dict in enumerate(feats):
            feat_dict["potentials"] = potentials[i]
        feats = x["additional_prolif_features"]
        for feat_dict in feats:
            # additional prolif features like pi staking simply
            # encoded as present/absent
            feat_dict["potentials"] = {feat_dict["interaction_type"]: 1}
        return x
