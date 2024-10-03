LOADING_MESSAGE = "If you leave this page you can come back later! You will get an email when this is done if it takes more than a minute."

# Slurm job states
PENDING = "PENDING"
RUNNING = "RUNNING"
CONFIGURING = "CONFIGURING"
COMPLETING = "COMPLETING"
FAILED = "FAILED"
COMPLETED = "COMPLETED"
UNKNOWN = "UNKNOWN"
# Slurm job state reasons
LAUNCH_FAILED_REASON = "JobHeldUser"
LAUNCH_FAILED_DESC = "launch failed requeued held"

SLURM_JOB_STATUS_CHECK_PERIOD = 60
SLURM_REST_API_USER = "basechem"

TORSION_ALERTS_TOTAL_ENERGY_THRESHOLD = 0.1
# This is the name of the property with torsion alerts data as returned by the Mayachemtools script
TORSION_ALERTS_PROP = "TorsionAlerts(RotBondIndices TorsionIndices TorsionAngle Energy EnergyLowerBoundCI EnergyUpperBoundCI HierarchyClass HierarchySubClass TorsionRule EnergyMethod AngleNotObserved MaxSingleEnergyAlert)"
TORSION_ALERTS_ENERGY_PROP = "TotalEnergy"
ADMIN_FAILURE = "Admin Failure"
ADMIN_NOTIFICATION = "Admin Notification"
