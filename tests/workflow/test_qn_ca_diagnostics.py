import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS


filename = "simulations/data/glass_final_300K.data"


# ------------------------------------------------------------
# Read system
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
# Isolated oxygen atoms
# ------------------------------------------------------------

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


print()
print("===== Ca Environment of Isolated O Atoms =====")


for oxygen_id in isolated_oxygen_ids:

    oxygen_index = np.where(
        system.ids == oxygen_id
    )[0][0]

    ca_neighbors = [
        j
        for j in neighbors.neighbors_of(oxygen_index)
        if system.types[j] == 3
    ]

    print()
    print(f"O ID : {oxygen_id}")
    print(
        f"Ca neighbors : "
        f"{[int(system.ids[j]) for j in ca_neighbors]}"
    )
    print(
        f"Ca coordination : "
        f"{len(ca_neighbors)}"
    )
