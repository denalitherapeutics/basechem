import json
import logging
from typing import Sequence

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from scipy.sparse import csr_array, hstack
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .feature_potentials import PotentialVectorizer


class ScoringFunctionModel:
    def __init__(
        self,
        potential_generator_pipeline: Pipeline,
        C: float = 1,
        version: str = None,
        protein_unsatisfied_dist_cutoff: float = 3.3,
        unsatisfied_frequency_cutoff: int = 100,
    ):
        """
        Model to generate interpretable scores for protein-ligand complexes.

        :param potential_generator_pipeline: a pipeline that takes a list of rdkit
            molecules and returns a list of dictionaries containing potentials
            used to score the complex. This pipeline will likely be a
            DistanceAngleTypingTransformer followed by a PotentialTransformer.
        :param C: Inverse of regularization strength parameter used by the internal
            logistic regression model. Lower values specify stronger regularization.
        :param version: An optional name for the model version.
        :param protein_unsatisfied_dist_cutoff: When evaluating receptor atoms that
            may be unsatisfied, only considers atoms within this distance of a ligand
            atom. Keeping this low avoid flagging irrelevant atoms far from the ligand.
        :param unsatisfied_frequency_cutoff: When building a dictionary of expected
            contributions for protein and ligand atom types, which is then used for
            evaluating when atoms are unsatisfied, only includes atom types that
            occur at least this many times in the training dataset. This avoids
            flagging atoms where we don't have a good estimate of their expected
            contribution.
        """
        self.potential_generator_pipeline = potential_generator_pipeline
        self.C = C
        self.version = version
        self.protein_unsatisfied_dist_cutoff = protein_unsatisfied_dist_cutoff
        self.unsatisfied_frequency_cutoff = unsatisfied_frequency_cutoff

    def update_receptor(self, receptor_path: str):
        """
        Change the receptor used by the scoring function.
        :param receptor_path: path to the new receptor file
        """
        self.potential_generator_pipeline.steps[0][1].update_receptor(receptor_path)

    def fit_from_files(
        self, feat_files: Sequence[str], labels: Sequence[int], complexes: Sequence[str]
    ):
        """
        Method to fit the model using a list of pre-computed feature file paths.

        :param feat_files: list of paths to feature files. Each file should contain
            a dictionary output by PotentialVectorizer for a single pose.
        * can we better clarify what the expected format would look like? is there
            a better format for the input data?

        :param labels: list of pose labels. 0 indicates near-native (RMSD<2.0),
            1 indicates non-near-native.
        :param complexes: list of complex names. Length should be the same as labels
            and feat_files.
        """

        # create weights for training, that ensure all that for each
        # complex the sum of positive example weights is 0.5 and the
        # sum of negative example weights is 0.5
        weights = _create_weights(complexes, labels)
        # vectorize potentials to create feature array
        self.potential_vectorizer = PotentialVectorizer()
        train_array = self.potential_vectorizer.fit_transform(feat_files)
        # add one-hot encoded complex names to the feature array
        onehot = OneHotEncoder(sparse_output=True, handle_unknown="ignore")
        onehot_array = onehot.fit_transform(np.array(complexes).reshape(-1, 1))
        full_train_array = hstack([train_array, onehot_array])
        full_train_array = csr_array(full_train_array)
        # fit the model
        self.model = LogisticRegression(C=self.C, solver="liblinear", max_iter=1000)
        self.model.fit(full_train_array, labels, sample_weight=weights)
        # "throw out" the one-hot encoded complex names
        n_keep_feats = train_array.shape[1]
        self.model.coef_ = self.model.coef_[:, :n_keep_feats]
        self.model.n_features_in_ = n_keep_feats
        # calculate expected contributions (i.e., expected sum of potentials) for
        # ligand and protein atom types, which are used to create
        # visuals of "unsatisfied" atoms
        self._calculate_expected_contribution_dicts(feat_files, labels, complexes)
        return self

    def _calculate_expected_contribution_dicts(
        self,
        feat_files: Sequence[dict],
        labels: Sequence[int],
        complexes: Sequence[str],
    ):
        """
        Calculate expected contributions for atom types based on provided train data.

        This is later used in interpretability to flag unsatisfied atoms. We calculate
        the average contributions of ligand and protein atoms of each type across
        low-RMSD poses in the training set, and save this with the model. For protein
        atoms, we only consider those within protein_unsatisfied_dist_cutoff of a
        ligand atom (to avoid flagging irrelevant atoms far from the ligand).

        The arguments are the same as for fit_from_files, and this is generally only
        called from within that function.
        """
        # find all our examples of low-RMSD poses and load their potentials
        df = pd.DataFrame(
            {"feat_file": feat_files, "complex": complexes, "label": labels}
        )
        df = df.query("label == 0").drop_duplicates("complex")
        low_rmsd_files = df["feat_file"].tolist()
        potentials = [json.load(open(f, "r")) for f in low_rmsd_files]
        # coalesce potentials into the total positive or negative potential
        # for each atom pair interaction
        dist_feature_dfs = [
            self._get_summed_interactions_df_from_dist_features(p) for p in potentials
        ]
        # get total contribution for each ligand atom
        ligand_sums = [
            d.groupby(["li", "lt"])["summed_interaction"].sum().reset_index()
            for d in dist_feature_dfs
            if len(d) > 0
        ]
        # get total contribution for each protein atom
        protein_sums = [
            d.groupby(["pi", "pt"])[["d", "summed_interaction"]]
            .agg({"d": "min", "summed_interaction": "sum"})
            .reset_index()
            for d in dist_feature_dfs
            if len(d) > 0
        ]
        for d in protein_sums:
            d.columns = ["pi", "pt", "min_d", "summed_interaction"]
        ligand_atoms = pd.concat(ligand_sums, ignore_index=True)
        protein_atoms = pd.concat(protein_sums, ignore_index=True)
        # preserve only those protein atoms that are sufficiently
        # close to the ligand (as defined by
        # protein_unsatisfied_dist_cutoff)
        protein_atoms = protein_atoms.query(
            "min_d <= " + str(self.protein_unsatisfied_dist_cutoff)
        )
        # calculate the average contribution for each atom type,
        # limiting to common atom types (as defined by
        # unsatisfied_frequency_cutoff)
        ligand_type_data = (
            ligand_atoms.groupby("lt")["summed_interaction"]
            .agg(["count", "mean"])
            .query("count >= " + str(self.unsatisfied_frequency_cutoff))
        )
        ligand_dict = {
            lt: ligand_type_data.loc[lt]["mean"] for lt in ligand_type_data.index
        }
        protein_type_data = (
            protein_atoms.groupby("pt")["summed_interaction"]
            .agg(["count", "mean"])
            .query("count >= " + str(self.unsatisfied_frequency_cutoff))
        )
        protein_dict = {
            pt: protein_type_data.loc[pt]["mean"] for pt in protein_type_data.index
        }
        # add the resulting mapping from atom types to expected contributions
        # to the model so it can be used for flaging unsatisfied atoms
        self.protein_contribution_dict = protein_dict
        self.ligand_contribution_dict = ligand_dict

    def _get_summed_interactions_df_from_dist_features(self, mol_potentials: dict):
        """
        Convert potentials to a dataframe with each atom pair's summed effect on
        score.

        This function loops through all the potentials acting between each pair
        of atoms and determines the total summed effect of the atom pair on the
        model's score. E.g., atoms A and B might have some negative repulsive
        potential, and some positive h-bond potential, etc, with a final summed
        effect on toklat_score of -0.274.

        :param mol_potentials: dictionary of potentials for a single pose.
        :return: DataFrame with li, pi, lt, pt, summed_interaction columns.
        """
        distance_features = mol_potentials["distance_features"]
        for d in distance_features:
            potential_vals = d["potentials"]
            summed_interaction = sum(
                [
                    v * self.model.coef_[0][self.potential_vectorizer.vocabulary_[k]]
                    if k in self.potential_vectorizer.vocabulary_
                    else 0
                    for k, v in potential_vals.items()
                ]
            )
            d["summed_interaction"] = summed_interaction
        return pd.DataFrame(distance_features)

    def predict_from_files(self, feat_files: Sequence[str]):
        """
        Predict scores using a list of pre-computed feature file paths.

        :param feat_files: list of paths to feature files.
        """
        feat_array = self.potential_vectorizer.transform(feat_files)
        out = self.model.decision_function(feat_array)
        return pd.DataFrame({"toklat_score": out})

    def predict_sdf(
        self, sdf_in_path: str, sdf_out_path: str, sort_by_score: bool = False
    ) -> str:
        """
        Predict scores for molecules in an SDF file, outputting results to SDF.

        :param sdf_in_path: path to the input SDF file.
        :param sdf_out_path: path to the output SDF file.
        :param sort_by_score: if True, sort the output DataFrame by score, with
            best poses first. Defaults to False.
        :return: path to the output SDF file.
        """
        mols_df = PandasTools.LoadSDF(sdf_in_path, removeHs=False, embedProps=True)
        scores_df = self.predict(mols_df["ROMol"].tolist(), sort_by_score=False)
        scores_df = pd.concat([mols_df, scores_df.drop(columns="ROMol")], axis=1)
        if sort_by_score:
            scores_df = scores_df.sort_values("toklat_score").reset_index(drop=True)
        PandasTools.WriteSDF(scores_df, sdf_out_path, properties=scores_df.columns)
        return sdf_out_path

    def predict(
        self, mols: Sequence[Chem.Mol], sort_by_score: bool = False
    ) -> pd.DataFrame:
        """
        Predict scores for a list of rdkit molecules.

        :param mols: list of rdkit molecules.
        :param sort_by_score: if True, sort the output DataFrame by score, with
            best poses first. Defaults to False.
        :return: DataFrame containing:
            - a "toklat_score" column
            - an "ROMol" column with the original molecules.
            - a "toklat_error" column. If the prediction failed, this column contains
              the error message. Empty string otherwise.
        """
        outs, _, failures = self._predict_mols(mols)
        df = pd.DataFrame(
            {"toklat_score": outs, "ROMol": mols, "toklat_error": failures}
        )
        if sort_by_score:
            df = df.sort_values("toklat_score").reset_index(drop=True)
        return df

    def predict_and_interpret_sdf(
        self,
        sdf_in_path: str,
        sdf_out_path: str,
        n_top_interactions: int = 5,
        max_n_unsatisfied_ligand_atoms: int = 3,
        max_n_unsatisfied_protein_atoms: int = 3,
        unsatisfied_diff_from_expected: float = 0.3,
        sort_by_score: bool = False,
    ) -> str:
        """
        Predict scores for molecules in an SDF file, outputting results to SDF.

        :param sdf_in_path: path to the input SDF file.
        :param sdf_out_path: path to the output SDF file.
        :param n_top_interactions: number of top interactions to return for each
            pose.
        :param max_n_unsatisfied_ligand_atoms: maximum number of unsatisfied ligand
            atoms to return for each pose. Unsatisfied atoms are those whose
            interaction contribution differs from the expected contribution for that
            atom type by more than unsatisfied_diff_from_expected.
        :param max_n_unsatisfied_protein_atoms: maximum number of unsatisfied protein
            atoms to return for each pose. Unsatisfied atoms are those whose
            interaction contribution differs from the expected contribution for that
            atom type by more than unsatisfied_diff_from_expected.
        :param unsatisfied_diff_from_expected: threshold for determining whether an
            atom is unsatisfied.
        :param sort_by_score: if True, sort the output DataFrame by score, with
            best poses first. Defaults to False.
        :return: path to the output SDF file. See docstring for predict_and_interpret
            for details on the output columns.
        :return: DataFrame containing
            - "toklat_score" column
            - "ROMol" column with the original molecules
            - "toklat_top_interactions" column. Each entry is a list of dicts with
                indices and types of the ligand and protein atom for the interaction,
                as well as the interaction total contribution
            - "toklat_unsatisfied_ligand_atoms" column. Each entry is a list of dicts
                with the indices and types of each unsatisfied ligand atom and the
                difference between the expected and actual atom contribution.
            - "toklat_unsatisfied_protein_atoms" column. Each entry is a list of dicts
                with the indices and types of each unsatisfied protein atom and the
                difference between the expected and actual atom contribution.
            - "toklat_error" column. If the prediction failed, this column contains the
                error message. Empty string otherwise.

        """
        mols_df = PandasTools.LoadSDF(sdf_in_path, removeHs=False, embedProps=True)
        scores_df = self.predict_and_interpret(
            mols_df["ROMol"].tolist(),
            n_top_interactions,
            max_n_unsatisfied_ligand_atoms,
            max_n_unsatisfied_protein_atoms,
            unsatisfied_diff_from_expected,
            sort_by_score=False,
        )
        scores_df = pd.concat([mols_df, scores_df.drop(columns="ROMol")], axis=1)
        scores_df["toklat_top_interactions"] = scores_df[
            "toklat_top_interactions"
        ].apply(json.dumps)
        scores_df["toklat_unsatisfied_ligand_atoms"] = scores_df[
            "toklat_unsatisfied_ligand_atoms"
        ].apply(json.dumps)
        if sort_by_score:
            scores_df = scores_df.sort_values("toklat_score").reset_index(drop=True)
        scores_df["toklat_unsatisfied_protein_atoms"] = scores_df[
            "toklat_unsatisfied_protein_atoms"
        ].apply(json.dumps)
        PandasTools.WriteSDF(
            scores_df,
            sdf_out_path,
            properties=scores_df.columns,
        )
        return sdf_out_path

    def predict_and_interpret(
        self,
        mols: Sequence[Chem.Mol],
        n_top_interactions: int = 5,
        max_n_unsatisfied_ligand_atoms: int = 3,
        max_n_unsatisfied_protein_atoms: int = 3,
        unsatisfied_diff_from_expected: float = 0.3,
        sort_by_score: bool = False,
    ) -> pd.DataFrame:
        """
        Predict scores for a list of rdkit molecules.

        :param mols: list of rdkit molecules.
        :param n_top_interactions: number of top interactions to return for each
            pose.
        :param n_top_ligand_atoms: number of top ligand atoms to return for each
            pose.
        :param sort_by_score: if True, sort the output DataFrame by score, with
            best poses first. Defaults to False.
        :return: DataFrame containing
            - "toklat_score" column
            - "ROMol" column with the original molecules
            - "toklat_top_interactions" column. Each entry is a list of dicts with
                indices and types of the ligand and protein atom for the interaction,
                as well as the interaction total contribution
            - "toklat_unsatisfied_ligand_atoms" column. Each entry is a list of dicts
                with the indices and types of each unsatisfied ligand atom and the
                difference between the expected and actual atom contribution.
            - "toklat_unsatisfied_protein_atoms" column. Each entry is a list of dicts
                with the indices and types of each unsatisfied protein atom and the
                difference between the expected and actual atom contribution.
            - "toklat_error" column. If the prediction failed, this column contains the
                error message. Empty string otherwise.
        """
        outs = []
        potentials = []
        failures = []
        interpretations = []
        outs, potentials, failures = self._predict_mols(mols)
        for potential, failed in zip(potentials, failures):
            if not failed:
                interpretations.append(
                    self._get_interpretation(
                        potential,
                        n_top_interactions,
                        max_n_unsatisfied_ligand_atoms,
                        max_n_unsatisfied_protein_atoms,
                        unsatisfied_diff_from_expected,
                    )
                )
            else:
                interpretations.append(([], [], []))
        top_interactions = [i[0] for i in interpretations]
        top_ligand_atoms = [i[1] for i in interpretations]
        top_protein_atoms = [i[2] for i in interpretations]
        df = pd.DataFrame(
            {
                "toklat_score": outs,
                "ROMol": mols,
                "toklat_top_interactions": top_interactions,
                "toklat_unsatisfied_ligand_atoms": top_ligand_atoms,
                "toklat_unsatisfied_protein_atoms": top_protein_atoms,
                "toklat_error": failures,
            }
        )
        if sort_by_score:
            df = df.sort_values("toklat_score").reset_index(drop=True)
        return df

    def _predict_mols(self, mols: Sequence[Chem.Mol]):
        """
        Attempts to predict mol scores in a batch, falling back to one-by-one if error.

        Returns list of scores, list of potentials dicts, and list of errors.
        """
        hit_exception = False
        try:
            for mol in mols:
                _validate_mol_has_hs(mol)
            potentials = self.potential_generator_pipeline.transform(mols)
            feat_array = self.potential_vectorizer.transform(potentials)
            outs = self.model.decision_function(feat_array)
            failed = [""] * len(mols)
        except Exception:
            # if there's a failure, fall back to predicting one at a time,
            # which is slower.
            hit_exception = True
        if hit_exception:
            outs = []
            potentials = []
            failed = []
            for mol in mols:
                out, potential, fail = self._predict_mol(mol)
                outs.append(out)
                potentials.append(potential)
                failed.append(fail)
        return outs, potentials, failed

    def _predict_mol(self, mol: Chem.Mol):
        """
        Generate a score for a single molecule.

        Returns score, potentials dict, and error string.
        """
        try:
            _validate_mol_has_hs(mol)
            potential = self.potential_generator_pipeline.transform([mol])
            feat_array = self.potential_vectorizer.transform(potential)
            out = self.model.decision_function(feat_array)
            failed = ""
        except Exception as e:

            out = [999]
            potential = [[]]
            error_str = str(e)
            if len(error_str) > 200:
                failed = error_str[:200] + "..."
            else:
                failed = error_str
            logging.exception("Error predicting molecule", exc_info=e)
        return out[0], potential[0], failed

    def _get_interpretation(
        self,
        mol_potentials: dict,
        n_top_interactions: int,
        max_n_unsatisfied_ligand_atoms: int,
        max_n_unsatisfied_protein_atoms: int,
        unsatisfied_diff_from_expected: float,
    ):
        """
        Extracts top interactions and unsatisfied atoms from a single pose's potentials.

        :param mol_potentials: dictionary of potentials for a single pose.
        :param n_top_interactions: number of top interactions to return.
        :param max_n_unsatisfied_ligand_atoms: maximum number of unsatisfied ligand
            atoms to return.
        :param max_n_unsatisfied_protein_atoms: maximum number of unsatisfied protein
            atoms to return.
        :param unsatisfied_diff_from_expected: threshold for determining whether an
            atom is unsatisfied.

        :return: tuple of: dict of top interactions, dict of unsatisfied ligand atoms,
            dict of unsatisfied protein atoms.
        """
        if n_top_interactions is None:
            n_top_interactions = len(mol_potentials["distance_features"])
        # get the summed contribution of each atom pair to the score
        feat_df = self._get_summed_interactions_df_from_dist_features(mol_potentials)
        if len(feat_df) == 0:
            return ([], [], [])
        # order by magnitude of summed_interaction
        feat_df["abs_summed_interaction"] = np.abs(feat_df["summed_interaction"])
        # sort by magnitude of summed_interaction - we'll use top n for top interactions
        feat_df = feat_df.sort_values(
            "abs_summed_interaction", ascending=False
        ).reset_index(drop=True)
        # calculate the total contribution of each ligand and protein atom
        ligand_totals = (
            feat_df.groupby(["li", "lt"])["summed_interaction"].sum().reset_index()
        )
        protein_totals = (
            feat_df.groupby(["pi", "pt"])[["d", "summed_interaction"]]
            .agg({"d": "min", "summed_interaction": "sum"})
            .reset_index()
        )
        protein_totals.columns = ["pi", "pt", "min_d", "summed_interaction"]
        # don't consider protein atoms that are too far from the ligand
        # for unsatisfied atom detection
        protein_totals = protein_totals.query(
            "min_d <= " + str(self.protein_unsatisfied_dist_cutoff)
        ).reset_index(drop=True)
        # check how different each atom's contribution is from the expected
        # value for that type. If the type is not in the expected contribution
        # dict, we'll use a large default contrib to ensure it's never flagged.
        ligand_totals["contrib_diff_from_expected"] = ligand_totals.apply(
            lambda x: x["summed_interaction"]
            - self.ligand_contribution_dict.get(x["lt"], 100),
            axis=1,
        )
        protein_totals["contrib_diff_from_expected"] = protein_totals.apply(
            lambda x: x["summed_interaction"]
            - self.protein_contribution_dict.get(x["pt"], 100),
            axis=1,
        )
        # get the top unsatisfied atoms and ensure they're unsatisfied enough,
        # based on unsatisfied_diff_from_expected, to flag.
        ligand_totals = (
            ligand_totals.sort_values("contrib_diff_from_expected", ascending=False)
            .reset_index(drop=True)
            .iloc[:max_n_unsatisfied_ligand_atoms]
        )
        protein_totals = (
            protein_totals.sort_values("contrib_diff_from_expected", ascending=False)
            .reset_index(drop=True)
            .iloc[:max_n_unsatisfied_protein_atoms]
        )
        ligand_totals = ligand_totals.query(
            f"contrib_diff_from_expected > {unsatisfied_diff_from_expected}"
        )
        protein_totals = protein_totals.query(
            f"contrib_diff_from_expected > {unsatisfied_diff_from_expected}"
        )

        return (
            feat_df.head(n_top_interactions)[
                ["li", "pi", "pt", "lt", "summed_interaction"]
            ].to_dict(orient="records"),
            ligand_totals[["li", "lt", "contrib_diff_from_expected"]].to_dict(
                orient="records"
            ),
            protein_totals[["pi", "pt", "contrib_diff_from_expected"]].to_dict(
                orient="records"
            ),
        )


def _validate_mol_has_hs(mol):
    """
    Ensure that the molecule has explicit hydrogens.
    """
    mol_with_hs = Chem.AddHs(mol)
    if mol.GetNumAtoms() != mol_with_hs.GetNumAtoms():
        raise ValueError("Molecule must have explicit hydrogens.")


def _create_weights(complexes, labels, weight_sum=1):
    """
    Weights positive and negative examples within each complex equally.

    :param weight_sum: the sum of weights for the entire dataframe
    """
    label_counter = {}
    for label, complex in zip(labels, complexes):
        if complex not in label_counter:
            label_counter[complex] = [0, 0]
        label_counter[complex][label] += 1
    weights_one = {}
    weights_zero = {}
    for complex, label_counts in label_counter.items():
        if label_counts[0] == 0 or label_counts[1] == 0:
            weights_one[complex] = 0
            weights_zero[complex] = 0
            continue
        weights_one[complex] = (weight_sum / label_counts[1]) / 2
        weights_zero[complex] = (weight_sum / label_counts[0]) / 2
    weights = [0] * len(labels)
    for i, (label, complex) in enumerate(zip(labels, complexes)):
        if label == 1:
            weights[i] = weights_one[complex]
        else:
            weights[i] = weights_zero[complex]
    return weights


def sdf_load_util(sdf_path: str) -> pd.DataFrame:
    """
    Load an sdf file output by the model, parsing columns to appropriate types.

    :param sdf_path: path to the SDF file.
    :return: Pandas DataFrame containing the SDF data.
    """
    df = PandasTools.LoadSDF(sdf_path, removeHs=False)
    if "toklat_score" in df.columns:
        df["toklat_score"] = df["toklat_score"].astype(float)
    if "toklat_top_interactions" in df.columns:
        df["toklat_top_interactions"] = df["toklat_top_interactions"].apply(json.loads)
    if "toklat_unsatisfied_ligand_atoms" in df.columns:
        df["toklat_unsatisfied_ligand_atoms"] = df[
            "toklat_unsatisfied_ligand_atoms"
        ].apply(json.loads)
    if "toklat_unsatisfied_protein_atoms" in df.columns:
        df["toklat_unsatisfied_protein_atoms"] = df[
            "toklat_unsatisfied_protein_atoms"
        ].apply(json.loads)
    return df
