import logging
import os
import subprocess
from shutil import move

from django.conf import settings
from django.core.mail import mail_admins
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdDistGeom, rdFingerprintGenerator, rdFMCS

from basechem.common.constants import ADMIN_FAILURE

logger = logging.getLogger("django")


class RDKitWrappers:
    """
    This is a helper class that contains util functions for completing
    analyses for RDKit mol objects
    """

    @staticmethod
    def generate_properties(mol):
        """
        Update the given mol object to include properties
        :prop mol: an rdkit mol object
        :return: the rdkit mol object with properties attached
        """
        Chem.SanitizeMol(mol)

        # Counts Properties
        n_o = Descriptors.NOCount(mol)
        mol.SetProp("acceptors", "%d" % n_o)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        mol.SetProp("aromaticrings", "%d" % aromatic_rings)
        charge = Chem.GetFormalCharge(mol)
        mol.SetProp("charge", "%d" % charge)
        num_chiral_centers = len(
            Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
        )
        mol.SetProp("chiralcenters", "%d" % num_chiral_centers)
        nh_oh = Descriptors.NHOHCount(mol)
        mol.SetProp("donors", "%d" % nh_oh)
        frac_c_sp3 = Descriptors.FractionCSP3(mol)
        mol.SetProp("fractioncsp3", "%.2f" % frac_c_sp3)
        heavy_atoms = Descriptors.HeavyAtomCount(mol)
        mol.SetProp("heavyatoms", "%d" % heavy_atoms)
        rings = Descriptors.RingCount(mol)
        mol.SetProp("rings", "%d" % rings)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        mol.SetProp("rotatablebonds", "%d" % rot_bonds)

        # Physiochemical Properties
        # MolLogP is listed as cLogP, as in basechem1.0
        mLogP = Descriptors.MolLogP(mol)
        mol.SetProp("clogp", "%.1f" % mLogP)
        mw = Descriptors.MolWt(mol)
        mol.SetProp("mw", "%d" % mw)
        sol_idx = mLogP + aromatic_rings
        mol.SetProp("solubilityindex", "%.1f" % sol_idx)
        tpsa = Descriptors.TPSA(mol)
        mol.SetProp("tpsa", "%d" % tpsa)

        return mol

    @staticmethod
    def confgen(mol, output_path, prunermsthresh=0.75, numconfs=50):
        """
        Runs RDKit conformer generation using ETKDGv2 which takes into consideration experimental
        torsion preferences (ET=Experimental Torsions; DG=Distance Geometry). This function prunes
        conformers that are within prunermsthresh and runs an additional minimization step after
        refinement
        More information on method: https://doi.org/10.1021/acs.jcim.5b00654
        :param mol: rdkit mol object of the compound to generate conformers for
        :param output_path: path to file where conformers will be written
        :param prunermsthresh: RMSD similarity cutoff for too similar conformers
        :param numconfs: max number of conformers to generate, rdkit recommends 50 for <= 7 rot bonds
        :return: path to file with conformers
        """
        mol = Chem.AddHs(mol, addCoords=True)
        # Set parameters to use for confgen
        param = rdDistGeom.ETKDGv3()
        param.pruneRmsThresh = prunermsthresh
        param.onlyHeavyAtomsForRMS = True
        param.enforceChirality = True
        param.useExpTorsionAnglePrefs = True
        param.randomSeed = 0xF00D
        # Create conformers
        cids = rdDistGeom.EmbedMultipleConfs(mol, numconfs, param)

        # Set up forcefield to calculate energy
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant="MMFF94s")
        AllChem.AlignMolConformers(mol)

        # Calculate energy based on MMFF and organize conformers by energies
        energies = []
        for cid in cids:
            try:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
                e = ff.CalcEnergy()
            except:
                mail_admins(
                    ADMIN_FAILURE,
                    f"MMFF failed to calculate energy for conformer: {cid} in the file {output_path}",
                )
                e = 100000  # set to something absurd
            energies.append((cid, e))

        sorted_energies = sorted(energies, key=lambda x: x[1])

        w = Chem.SDWriter(output_path)
        for cid, e in sorted_energies:
            mol.SetProp("s_conf_id", str(cid))
            mol.SetProp("r_mmff_energy", str(e))
            w.write(mol, confId=cid)
        w.close()

        mol.SetProp("s_low_e_conf_id", str(sorted_energies[0][0]))
        mol.SetProp("r_mmff_low_energy", str(sorted_energies[0][1]))

        return output_path, mol

    def is_3d(mol):
        """
        Determine if the given mol object has 3D coordinates. `mol.GetConformer().Is3D()`
        works sometimes but is not always correct as it is just a flag that must be set.
        :prop mol: an rdkit mol object
        :return: a boolean, does this mol object have 3D coordinates?
        """
        x_coords = set()
        y_coords = set()
        z_coords = set()
        try:
            conf = mol.GetConformer()
        except ValueError:
            # Some mol objects do not have a conformer set (ex, mols generated from smiles),
            # which raises a ValueError. These mols are always 2D
            return False
        for atom_idx in range(mol.GetNumAtoms()):
            atom_pos = conf.GetAtomPosition(atom_idx)
            x_coords.add(atom_pos.x)
            y_coords.add(atom_pos.y)
            z_coords.add(atom_pos.z)
        return (len(x_coords) > 1) and (len(y_coords) > 1) and (len(z_coords) > 1)

    @staticmethod
    def clean_mol_object(mol):
        """
        Helper to clean a mol object to check if there are multiple compounds
        :param mol: ROMol Object
        :return: List of tuples of the form (mol, twoD) where
            - mol is an RDKit ROMol object
            - twoD is a boolean, is this mol object 2D?
        """
        romols = []
        # Splitting out disconnected fragments and treating each as their own mol object
        try:
            fragments = Chem.GetMolFrags(mol, asMols=True)
            for frag in fragments:
                twoD = not RDKitWrappers.is_3d(mol)
                # copy any additional tags to each new mol
                sd_tag_names = mol.GetPropNames()
                for sd_tag in sd_tag_names:
                    tag_val = mol.GetProp(sd_tag)
                    frag.SetProp(sd_tag, tag_val)
                romols.append((frag, twoD))
        except Exception as e:
            logger.error(f"ROMol object could not be fragmented: {e}")
        return romols

    @staticmethod
    def ligprep(mol, output_filepath):
        """
        RDKit based implementation of Schrodingers "LigPrep" leveraging Gypsum-DL by the Durant Lab
        `https://github.com/durrantlab/gypsum_dl.git`
        Generates ionization states, enumerates tautomers, provides different chiral states, adds Hs,
        and does basic confgen to generate a single 3D conformer for each enumerated structure
        :param mol: rdkit object of the compound to generate 3D structures for
        :param output_filepath: path to sdf file to write final structures
        :return: path to file with expanded structures
        """
        mol_name = mol.GetProp("_Name")
        dirname = os.path.dirname(output_filepath)
        # Gypsym always assigns the output name to includ `__input1`
        gypsum_outfile = f"{dirname}/{mol_name}__input1.sdf"

        # Write mol to smiles file for Gypsum-DL to use
        smiles_file = f"{dirname}/{mol_name}.smi"
        smiles = Chem.MolToSmiles(mol)
        with open(smiles_file, "w") as f:
            f.write(f"{smiles} \t{mol.GetProp('_Name')}")

        cmd = [
            "python",
            f"{settings.BASE_DIR}/common/gypsum_dl/run_gypsum_dl.py",
            "--source",
            smiles_file,
            "--use_durrant_lab_filters",
            "--output_folder",
            dirname,
            "--separate_output_files",
        ]

        process = subprocess.run(
            args=cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if process.returncode == 0:
            os.rename(gypsum_outfile, output_filepath)
            return output_filepath

    @staticmethod
    def superimpose(query_sdf, template_sdf, output_sdf):
        """
        Superimpose (rigid align) the given conformers onto a template molecule and set
        RMSD of each conformer to the template
        :param query_sdf: path to file with conformers of a single mol to superimpose
        :param template_sdf: path to file with template molecule
        :param output_sdf: path to file to write the superimposed conformers
        :return: path to file with superimposed conformers
        """
        query_mols = [m for m in Chem.SDMolSupplier(query_sdf, removeHs=False)]
        template = Chem.SDMolSupplier(template_sdf, removeHs=False)[0]

        mcs = rdFMCS.FindMCS(
            [template, query_mols[0]],  # all input mols should be the same structure
            threshold=0.8,
            completeRingsOnly=False,
            ringMatchesRingOnly=False,
        )
        mcs_patt = Chem.MolFromSmarts(mcs.smartsString)
        template_match = template.GetSubstructMatch(mcs_patt)

        w = Chem.SDWriter(output_sdf)
        save_template = True
        for mol in query_mols:
            # Try to align the entire Mol but if it fails, align based on MSC
            try:
                rmsd = AllChem.AlignMol(mol, template)
                mol.SetProp("r_bc_rmsd_to_lsalign", str(rmsd))
                if save_template:
                    # If successful, save the template to the file as well because it is also
                    # an exact conformer match to itself
                    template.SetProp("r_bc_rmsd_to_lsalign", "0")
                    template.SetProp("isTemplate", "True")
                    w.write(template)
                    save_template = False
            except RuntimeError:
                conf_match = mol.GetSubstructMatch(mcs_patt)
                rmsd = AllChem.AlignMol(
                    mol, template, atomMap=list(zip(conf_match, template_match))
                )
                mol.SetProp("r_bc_rmsd_to_lsalign", str(rmsd))
            w.write(mol)

        w.close()
        return output_sdf

    def pick_series(compound_mol, series_tuples):
        """
        Finds most similar series to the given compound
        :param compound_mol: RDKit mol object of the query compound
        :param series_tuples: list of tuples of (RDKit mol object, name) to use as references
        :return: name of series with the highest similarity to the query compound
        """
        rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
        comp_fp = rdkgen.GetFingerprint(compound_mol)
        series_fp = [rdkgen.GetFingerprint(s[0]) for s in series_tuples]

        try:
            sim = []
            for s_fp in series_fp:
                s = DataStructs.TanimotoSimilarity(comp_fp, s_fp)
                sim.append(s)
            i = sim.index(max(sim))
            return series_tuples[i][1]

        except:
            return ""

    def add_hydrogens_to_sdf(filepath):
        """
        Given an SDF file, add hydrogens to all molecules in the file. The original file is overwritten.
        :param filepath: a string, the path to the SDF file
        :returns: a string, the path to the SDF file with hydrogens added
        """
        mols = [m for m in Chem.SDMolSupplier(filepath, removeHs=False)]
        tmp_filepath = f"{os.path.splitext(filepath)[0]}_Hs.sdf"
        writer = Chem.SDWriter(tmp_filepath)
        for mol in mols:
            if mol:
                mol_w_hs = Chem.AddHs(mol, addCoords=True)
                if mol.GetNumAtoms() < mol_w_hs.GetNumAtoms():
                    writer.write(mol_w_hs)
                else:
                    writer.write(mol)
        writer.close()
        move(tmp_filepath, filepath)
        return filepath
