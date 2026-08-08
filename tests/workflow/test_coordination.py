from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.coordination import CoordinationAnalysis
from workflow.parameters import BOND_CUTOFFS


filename = "simulations/data/glass_final_300K.data"


# ---------------------------------------------------------
# Read structure
# ---------------------------------------------------------

system = read_lammps(filename)


# ---------------------------------------------------------
# Build neighbors
# ---------------------------------------------------------

neighbors = NeighborList(
    system,
    BOND_CUTOFFS
)

neighbors.build()


# ---------------------------------------------------------
# Coordination analysis
# ---------------------------------------------------------

coord = CoordinationAnalysis(
    system,
    neighbors
)


coord.report()
