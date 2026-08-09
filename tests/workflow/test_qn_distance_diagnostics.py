import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.geometry import distance
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


qn = QnAnalysis(
    system,
    neighbors
)


# ------------------------------------------------------------
# Five-coordinate Si atoms
# ------------------------------------------------------------

print()
print("===== Five-coordinate Si atoms =====")

for si_id in [354, 489]:

    si_index = np.where(
        system.ids == si_id
    )[0][0]

    print()
    print(f"Si ID : {si_id}")
    print(f"Position : {system.positions[si_index]}")

    oxygen_neighbors = [
        j
        for j in neighbors.neighbors_of(si_index)
        if system.types[j] == 2
    ]

    distances = []

    for j in oxygen_neighbors:

        d = distance(
            system.positions[si_index],
            system.positions[j],
            system.box_lengths
        )

        distances.append(
            (system.ids[j], d)
        )

    distances.sort(key=lambda x: x[1])

    for oxygen_id, d in distances:

        print(
            f"O ID {oxygen_id:5d} : "
            f"{d:.6f} Angstrom"
        )


# ------------------------------------------------------------
# Isolated oxygen atoms
# ------------------------------------------------------------

print()
print("===== Isolated oxygen atoms =====")

isolated_oxygen_ids = [
    5095,
    2507,
    4764,
    4541,
    4375,
    1701,
    5202,
    3826,
]


for oxygen_id in isolated_oxygen_ids:

    oxygen_index = np.where(
        system.ids == oxygen_id
    )[0][0]

    print()
    print(f"O ID : {oxygen_id}")
    print(
        f"Position : "
        f"{system.positions[oxygen_index]}"
    )

    si_distances = []

    for j in range(system.natoms):

        if system.types[j] != 1:
            continue

        d = distance(
            system.positions[oxygen_index],
            system.positions[j],
            system.box_lengths
        )

        si_distances.append(
            (system.ids[j], d)
        )

    si_distances.sort(key=lambda x: x[1])

    print("Three nearest Si atoms:")

    for si_id, d in si_distances[:3]:

        print(
            f"Si ID {si_id:5d} : "
            f"{d:.6f} Angstrom"
        )
