from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.modifiers import ModifierAnalysis
from workflow.parameters import BOND_CUTOFFS


system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

neighbors = NeighborList(
    system,
    BOND_CUTOFFS
)

neighbors.build()

modifiers = ModifierAnalysis(
    system,
    neighbors
)

modifiers.report()
