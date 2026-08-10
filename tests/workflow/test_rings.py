from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.topology import TopologyAnalysis
from workflow.rings import RingAnalysis
from workflow.parameters import BOND_CUTOFFS


system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

neighbors = NeighborList(
    system,
    BOND_CUTOFFS
)

neighbors.build()

topology = TopologyAnalysis(
    system,
    neighbors
)

rings = RingAnalysis(
    system,
    topology
)

rings.report()
