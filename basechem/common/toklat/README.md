# toklat
Inductive Bio and Denali collaboration

![image](https://github.com/denalitherapeutics/toklat/assets/9731569/5ecd5388-f86b-4ace-9a1b-e9051fb4bcd4)

## Installation

The package can be installed by running the following command in
the root directory of the repository. **Note: the code has been tested with python
3.10 and may not work in python 3.9 and below**.

```bash
pip install -e ."[dev]"
```

If you don't need to run the pre-commit hooks or use an ipython kernel, you can
omit the `[dev]` part of the command.

A full set of tested dependencies can be found in the `requirements.txt` file.

## Philosophy/goals of the toklat package

The toklat package aims to allow users to train and use interpretatable pose
scoring functions for protein-ligand interactions. The package is designed to:
1. Be highly configurable to support the extraction of a wide variety of
   features from the protein-ligand complex, add new feature extraction code
   with python, as well as allow incorporation of pre-calculated features from
   outside toklat.
   * (This contrasts with other docking programs, which generally offer a fixed
     set of partially configurable feature extractors and are not easily
     extensible.)
2. Enable the training of scoring functions from these features, using as
    training data a set of examples of near-native and non-near-native poses.
    * (This contrasts with other scoring functions, which are typically trained
    using affinity data rather than being trained to distinguish near-native and
    non-near-native poses. Training models for pose scoring based on affinity
    data has [various](https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c00049)
    [drawbacks]( https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155183).)
3. Provide a scoring function that is interpretable, meaning that the full set
    of input features is inspectable and attributable at the atomic level, and
    this can be converted into user-facing visualizations of the positive and
    negative contributions of a pose's score.
    * (This contrasts with other scoring functions, which typically provide a
        limited breakdown of contributions to a score, e.g. total contribution
        for each potential type or each residue.)

## Usage

Example usage of a trained SF can be found in the `example_usage` directory. This shows
how to load a trained model, use it to score poses, and generate a visualization of
contributions.

**Note**: By following the example notebook in the `example_usage` directory, you should be
able to make use of a trained model to score poses and generate visualizations. The detailed information
provided below in the "Package structure deep dive" section is helpful for diving deeper into the codebase,
but is not a prerequisite for general usage.

### Useful functions for general usage

Examples of these key functions are provided in the `example_usage/predict_and_visualize.ipynb` notebook.

* [`ScoringFunctionModel.update_receptor`]: Choose which receptor file to use for scoring poses. Receptor should have hydrogens and not have a ligand present. See docstring in [`toklat/model.py`](toklat/model.py).
* [`ScoringFunctionModel.predict_and_interpret_sdf`]: Generate scores and interpretability information for a set of poses in an SDF file. See docstring in [`toklat/model.py`](toklat/model.py). Options are provided for choosing how many top interactions or unsatisfied atoms to show, how strict to be in flagging unsatisfied atoms, and whether to sort the output SDF by score (which defaults to no sorting).
* [`sdf_load_util`]: Convenience function for loading SDF files created by `predict_and_interpret_sdf`. See docstring in [`toklat/model.py`](toklat/model.py).
* [`visualize_interpretation`]: Generate a visualization of the contributions to a pose's score. See docstring in [`toklat/viz_utils.py`](toklat/viz_utils.py).

## Package structure deep dive

The four primary modules within toklat are:
1. `toklat.feature_generation`: Contains code for extracting features from protein-ligand
    complexes.
2. `toklat.feature_potentials`: Contains code for converting features into a numerical
    potentials, and code for converting numerical potentials into a vector.
3. `toklat.model`: Contains code for `ScoringFunctionModel`, which is the central class for training
    and using a scoring function. Trained `ScoringFunctionModel` objects can score poses and can also
    output interpretability information for visualization.
4. `toklat.viz_utils`: Contains code for generating py3dmol visualizations.

### Feature generation (feature_generation module)

The main class here is `DistanceAngleTypingTransformer`, which essentially takes
a path to a protonated receptor (without a ligand present) and a protonated rdkit Mol ligand, and returns a dictionary containing all
of the extracted features for the protein-ligand complex.

The transformer supports three types of features:
* `additional_ligand_features`: arbitrary additional ligand features (e.g., torsion strain, LogD, etc)
  can be generated by passing functions to `additional_ligand_feature_functions` when
  constructing the `DistanceAngleTypingTransformer`. These are expected to take
  a `rdkit.Mol` and return a float. One notable type of additional ligand feature function is `PropAsFloatGetter`, which
  extracts a property that has already been added to the ligand (e.g., by inclusion in the SDF file). This is uses to extract
  pre-calculated torsion strain in the current model. `PropAsFloadGetter` will print a warning if an expected property is missing.
* `distance_features`: These store information about all pairwise protein atom/ligand atom interactions (up to a distance cutoffs, which defaults to 6A). They store the index, atom types, distance, and angle (for H bonds only) between the atoms.
* `additional_prolif_features`: These store information about additional prolif interactions, like pi stacking.

The atom typing for `distance_features` is defined by `_ATOM_TYPING_SMARTS` and `get_atom_types`, with some special code in `get_atom_types_protein` to deal with issues parsing pdb files in rdkit.

`additional_prolif_features`, as well as hydrogen bonds, are extracted by the `prolif` library. They are configured in `_DEFAULT_FP`.

Confusing quirk to watch out for: H bond angles are extracted by `prolif` but then are added in to the `distance_features` dictionary as an `ang` term for a given pair of atoms rather than being saved separately in `additional_prolif_features`.

Additional note: This code *in theory* supports the high configurability described in point (1) of the philosophy, but in practice it has been written with a focus on getting a single working version of the model in place. `_ATOM_TYPING_SMARTS` and `_DEFAULT_FP` are hardcoded in ways that you wouldn't want them to be in an easily configurable system, which could be fixed with some refactoring.

An example of what the output of `DistanceAngleTypingTransformer` looks like can be found below. This shows an example for one ligand. For multiple ligands, we would have a list with multiple of these dictionaries.

```python
{
    "distance_features": [
        # list of dictionaries for all pairwise protein-ligand interactions
        {
            "li": 1, # ligand atom index
            "pi": 123, # protein atom index
            "distance": 2.91 # distance between the atoms
            "ang": 161, # DHA angle for H-bonds
            "lt": "O_hba", # ligand atom type
            "pt": "N_hbd", # protein atom type
        }
        # ... additional pairwise distance features
    ],
    "additional_ligand_features": {
        # dictionary with any calculated additional ligand features for the pose
        "MaxSingleEnergy": 1.3,
        "TotalEnergy": 12
        # ... any additional additional ligand features
    },
    "additional_prolif_features": [
        # information about any features like pi stacking
        {
            "interaction_type": "PiStacking",
            "distance": 5.6,
            # other info output by prolif...
        }
        # ... any additional prolif interactions extracted
    ]
}
```

### Feature potentials (feature_potentials module)

The job of the main class in this module, `PotentialTransformer`, is to take the output
of `DistanceAngleTypingTransformer` which mostly contains just distances, atom type pairs, and angles, and
convert this into numeric potentials including repulsion, van der waal interactions, and h bonds.

It allows the user to control what potentials are extracted, and how they are parameterized, by passing a list of `potential_functions` to the `PotentialTransformer` constructor.

The default potentials are listed in `_DEFAULT_POTENTIAL_FUNCTIONS`, and include a VDW potential, a repulsive potential, an H bond potential, a hydrophobic potential,
and a pairwise potential. All but the last are designed similar to the potentials
in [AutoDock Vina](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21334), though the H bond potential has an added angle-dependent component. The pairwise potential is new and lets us learn a more complex scoring function: it lets us learn a different contribution for each pair of atom types. E.g., it would let us learn that Fluorines don't like to be next to Oxygen H-bond acceptors.

The potential functions are highly extensible, and are defined in `feature_potential_helpers.py`: one can create any object that has a `__call__` signature taking a distance, and angle, a protein atom type, and a ligand atom type as arguments and returning a dictionary. See `feature_potential_helpers.py` for details and examples.

When transforming data received from a `DistanceAngleTypingTransformer`, the `PotentialVectorizer` will add the calculated potentials into the data structure and return it. This results in a data structure like this:

```python
{
    "distance_features": [
        # list of dictionaries for all pairwise protein-ligand interactions
        {
            "li": 1, # ligand atom index
            "pi": 123, # protein atom index
            "distance": 2.91 # distance between the atoms
            "ang": 161, # DHA angle for H-bonds
            "lt": "O_hba", # ligand atom type
            "pt": "N_hbd", # protein atom type
            "potentials": {
                # these are the potentials one would expect to see if using _DEFAULT_POTENTIAL_FUNCTIONS
                # for this pair of atoms
                "vdw": 1, # there are potentials corresponding to potentials in vina
                "hbond_dist": 0.9, # there are angle-dependent and non-angle dependent H-bond potentials
                "hbond_dist_angle": 0.9,
                "repulsion": 0.43,
                "O_hba|N_hbd": 0.7 # there are potentials for specific atom type pairs
            }
        }
        # ... additional pairwise distance features
    ],
    "additional_ligand_features": {
        # dictionary with any calculated additional ligand features for the pose
        "MaxSingleEnergy": 1.3,
        "TotalEnergy": 12
        # ... any additional additional ligand features
    },
    "additional_prolif_features": [
        {
            "interaction_type": "PiStacking",
            "distance": 5.6,
            # other info output by prolif...
            "potentials": {
                "PiStacking": 1
            }
        }
        # ... any additional prolif interactions extracted
    ]
}
```

These data structures can then be passed to a `PotentialVectorizer`, which simply converts the potentials into a single vector for each pose, where each element is the total value of a potential for the pose (e.g., the summed repulsive potential across all interactions in the pose, or the summed H-bond potential).

### Model training and prediction (model module)

The primary class in this module, `ScoringFunctionModel`, has two jobs. It can be used to train a scoring function from a set of poses, via `fit_from_files`. And the trained model can be used to score poses via `predict` as well as interpretability information via `predict_and_interpret`.

#### Training

The training process in `fit_from_files` is mainly just a Logistic Regression on poses where the labels are 0 for near-native and 1 for non-near-native. However, along with the labels, we also pass a `complexes` vector with complex IDs for each pose, and these are used to modify the regression in two ways.

* First, the poses are given weights such that all the non-near-native poses have weights that add up to 0.5, and all the near-native poses have weights that add up to 0.5. In other words, the total weight of poses for the complex is 1, and non-near-native and near-native poses are weighted equally.

* Second, the regression is modified to include a complex-specific bias term. This is done by adding a column to the regression matrix for each complex which is "1" for poses from that complex and "0" for poses from any other complex. (These coefficients are thrown out when saving the model, since at prediction time we will have new complexes and will not provide a complex ID.)

These modifications are both designed to encourage the model to learn a scoring function that distinguishes near-native from non-near-native poses *within a complex* rather than "cheating" by learning what types of complexes are easier or more likely to have near-native poses. Both modifications appear to be essential to having well-performing models.

#### Prediction

The `predict` and `predict_sdf` methods are relatively straightforward, loading molecules and running them through the model pipeline to produce scores, with some error checking along the way.

The `predict_and_interpret` and `predict_and_interpret_sdf` methods are a bit more complex. Prediction on its own doesn't easily expose the exact interactions contributing to a score, so there is a special `_get_interpretation` function that takes the raw potentials for the pose (before they are vectorized) and uses them to calculate which pairwise interactions contribute most to the score.

### Visualizing predictions (viz_utils module)

This module contains a `visualize_interpretation` function that takes the interpretability output from `predict_and_interpret` and generates a visualization using py3Dmol.
