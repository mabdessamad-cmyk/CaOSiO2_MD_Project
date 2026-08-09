from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.bond_angles import BondAngleAnalysis
from workflow.parameters import BOND_CUTOFFS


FILENAME = "simulations/data/glass_final_300K.data"

system = read_lammps(FILENAME)

neighbors = NeighborList(
    system,
    BOND_CUTOFFS
)

neighbors.build()

angles = BondAngleAnalysis(
    system,
    neighbors
)

angles.report()
