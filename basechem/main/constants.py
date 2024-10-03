MAX_N_CONFS_DISPLAY = 8
MIN_REL_ENERGY = 0.6
MAX_GEMS_PER_CO = 3
MAX_TORSION_ATTEMPTS = 3

RDKIT_2D_WARNING = "Warning: molecule is tagged as 3D, but all Z coords are zero"

# Project codes
AUTO = "AUTO"
TEST = "TEST"

# Types of Analyses
PROPCALC = "propcalc"
ALIGN = "align"
DOCK = "dock"
ESP = "esp"
TORSION = "torsion"
MMP = "mmp"
DTX_MMP = "dtx-mmp"
USER_UPLOADED = "user uploaded"

SAVED_FROM_OPTIONS = [
    (ALIGN, ALIGN),
    (DOCK, DOCK),
    (TORSION, TORSION),
    (USER_UPLOADED, USER_UPLOADED),
]

# Django Q task statuses
DROPPED = "dropped"
IN_PROGRESS = "in progress"
COMPLETE = "complete"
ERROR = "error"

# Propcalc Properties
ACCEPTORS = "Acceptors"
AROMATIC_RINGS = "AromaticRings"
CHARGE = "Charge"
DONORS = "Donors"
FRACTIONCSP3 = "FractionCSP3"
RINGS = "Rings"
ROTATABLEBONDS = "RotatableBonds"
HEAVYATOMS = "HeavyAtoms"
MW = "MW"
SOLUBILITYINDEX = "SolubilityIndex"
TPSA = "TPSA"
CLOGP = "cLogP"
ALOGD = "LogD"
RLM = "RLM"
HLM = "HLM"

ALL_PROPS = [
    ACCEPTORS,
    AROMATIC_RINGS,
    CHARGE,
    DONORS,
    FRACTIONCSP3,
    HEAVYATOMS,
    MW,
    RINGS,
    ROTATABLEBONDS,
    SOLUBILITYINDEX,
    TPSA,
    CLOGP,
    ALOGD,
    RLM,
    HLM,
]

# Max hours to wait for sending an assay email
MAX_HOURS_TO_WAIT = 4

# Compound DTX properties to save in compound measured_data field
DTX_PROPERTIES_TO_SAVE = ["HLM", "RLM", "LogD Avg", "Assay Result"]

# Max property values for Compound to be included as an MMP
MMP_MAX_MW = 550
MMP_MAX_LOGP = 5.5
MMP_MAX_TPSA = 150
MMP_MAX_H_BOND_DONORS = 4
