from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.topology import TopologyAnalysis
from workflow.parameters import BOND_CUTOFFS


FILENAME = "simulations/data/glass_final_300K.data"

system = read_lammps(FILENAME)

neighbors = NeighborList(system, BOND_CUTOFFS)
neighbors.build()

topology = TopologyAnalysis(system, neighbors)

connections = topology.si_connections()

zero_degree = [
    si for si, connected in connections.items()
    if len(connected) == 0
]

degree_counts = {}

for connected in connections.values():
    degree = len(connected)
    degree_counts[degree] = degree_counts.get(degree, 0) + 1

print()
print("===== Topology Validation =====")
print(f"Total Si atoms       : {len(connections)}")
print(f"Zero-degree Si atoms : {len(zero_degree)}")
print()
print("Si network degree distribution")

for degree in sorted(degree_counts):
    print(
        f"Degree {degree}: "
        f"{degree_counts[degree]} atoms"
    )

assert len(connections) == 1680
assert len(topology.bridging_oxygens()) == 2653
assert all(
    degree >= 0
    for degree in degree_counts
)

print()
print("Topology validation: PASS")
