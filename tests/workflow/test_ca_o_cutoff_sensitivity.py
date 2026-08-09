import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS


FILENAME = "simulations/data/glass_final_300K.data"

CA_O_CUTOFFS = [
    2.90, 3.00, 3.10, 3.20,
    3.30, 3.40, 3.50, 3.60
]

system = read_lammps(FILENAME)

print()
print("===== Ca-O Cutoff Sensitivity =====")
print()

print(
    f"{'Cutoff':>8} "
    f"{'Ca mean':>10} "
    f"{'Ca min':>8} "
    f"{'Ca max':>8}"
)

print("-" * 42)

for ca_o_cutoff in CA_O_CUTOFFS:

    cutoffs = dict(BOND_CUTOFFS)

    # Keep validated Si-O cutoff
    cutoffs[(1, 2)] = 2.30
    cutoffs[(2, 1)] = 2.30

    # Change only Ca-O
    cutoffs[(3, 2)] = ca_o_cutoff
    cutoffs[(2, 3)] = ca_o_cutoff

    neighbors = NeighborList(system, cutoffs)
    neighbors.build()

    ca_indices = np.where(system.types == 3)[0]

    ca_coordination = np.array([
        sum(
            system.types[j] == 2
            for j in neighbors.neighbors_of(i)
        )
        for i in ca_indices
    ])

    print(
        f"{ca_o_cutoff:8.2f} "
        f"{np.mean(ca_coordination):10.3f} "
        f"{np.min(ca_coordination):8d} "
        f"{np.max(ca_coordination):8d}"
    )
