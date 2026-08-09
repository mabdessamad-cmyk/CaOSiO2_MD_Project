from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.qn import QnAnalysis


filename = "simulations/data/glass_final_300K.data"


# ------------------------------------------------------------
# Read system
# ------------------------------------------------------------

system = read_lammps(filename)


# ------------------------------------------------------------
# Build neighbors
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


# ------------------------------------------------------------
# Find anomalous oxygen atoms
# ------------------------------------------------------------

print()
print("===== Anomalous Oxygen Atoms =====")

oxygen_indices = [
    i for i in range(system.natoms)
    if system.types[i] == 2
]

other_oxygen = []

for i in oxygen_indices:

    si_neighbors = qn.si_neighbors_of_oxygen(i)

    if len(si_neighbors) not in (1, 2):

        other_oxygen.append(i)

        print()
        print(f"Atom index : {i}")
        print(f"Atom ID    : {system.ids[i]}")
        print(f"Position   : {system.positions[i]}")
        print(f"Si neighbors : {si_neighbors}")
        print(f"Number of Si neighbors : {len(si_neighbors)}")


# ------------------------------------------------------------
# Find anomalous Si atoms
# ------------------------------------------------------------

print()
print("===== Anomalous Si Atoms =====")

si_indices = [
    i for i in range(system.natoms)
    if system.types[i] == 1
]

invalid_si = []

for i in si_indices:

    oxygen_neighbors = [
        j
        for j in neighbors.neighbors_of(i)
        if system.types[j] == 2
    ]

    qn_value = qn.qn_of_si(i)

    if qn_value not in (0, 1, 2, 3, 4):

        invalid_si.append(i)

        print()
        print(f"Atom index : {i}")
        print(f"Atom ID    : {system.ids[i]}")
        print(f"Position   : {system.positions[i]}")
        print(f"O neighbors : {oxygen_neighbors}")
        print(f"Number of O neighbors : {len(oxygen_neighbors)}")
        print(f"Qn value : {qn_value}")


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print()
print("===== Diagnostic Summary =====")
print(f"Other oxygen atoms : {len(other_oxygen)}")
print(f"Invalid Si atoms   : {len(invalid_si)}")
