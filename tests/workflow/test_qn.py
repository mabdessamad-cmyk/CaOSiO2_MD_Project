from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.qn import QnAnalysis


filename = "simulations/data/glass_final_300K.data"


# ------------------------------------------------------------
# Read atomic system
# ------------------------------------------------------------

system = read_lammps(filename)


# ------------------------------------------------------------
# Build neighbor list
# ------------------------------------------------------------

neighbors = NeighborList(
    system,
    BOND_CUTOFFS
)

neighbors.build()


# ------------------------------------------------------------
# Qn analysis
# ------------------------------------------------------------

qn = QnAnalysis(
    system,
    neighbors
)


qn.report()
