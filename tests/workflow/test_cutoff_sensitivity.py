import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.qn import QnAnalysis
from workflow.parameters import BOND_CUTOFFS


FILENAME = "simulations/data/glass_final_300K.data"

# Si-O cutoffs to test
SI_O_CUTOFFS = [2.00, 2.10, 2.20, 2.30, 2.40, 2.50]


system = read_lammps(FILENAME)

print()
print("===== Si-O Cutoff Sensitivity =====")
print()

print(
    f"{'Cutoff':>8} "
    f"{'Si mean':>10} "
    f"{'Si min':>8} "
    f"{'Si max':>8} "
    f"{'Q0':>6} "
    f"{'Q1':>6} "
    f"{'Q2':>6} "
    f"{'Q3':>6} "
    f"{'Q4':>6} "
    f"{'Invalid':>8}"
)

print("-" * 80)


for si_o_cutoff in SI_O_CUTOFFS:

    cutoffs = dict(BOND_CUTOFFS)

    cutoffs[(1, 2)] = si_o_cutoff
    cutoffs[(2, 1)] = si_o_cutoff

    neighbors = NeighborList(
        system,
        cutoffs
    )

    neighbors.build()

    # --------------------------------------------------------
    # Si-O coordination
    # --------------------------------------------------------

    si_indices = np.where(system.types == 1)[0]

    si_coordination = np.array([
        sum(
            system.types[j] == 2
            for j in neighbors.neighbors_of(i)
        )
        for i in si_indices
    ])

    # --------------------------------------------------------
    # Qn analysis
    # --------------------------------------------------------

    qn = QnAnalysis(
        system,
        neighbors
    )

    distribution = qn.qn_distribution()

    q0 = distribution.get(0, 0)
    q1 = distribution.get(1, 0)
    q2 = distribution.get(2, 0)
    q3 = distribution.get(3, 0)
    q4 = distribution.get(4, 0)
    invalid = distribution.get("invalid", 0)

    print(
        f"{si_o_cutoff:8.2f} "
        f"{np.mean(si_coordination):10.3f} "
        f"{np.min(si_coordination):8d} "
        f"{np.max(si_coordination):8d} "
        f"{q0:6d} "
        f"{q1:6d} "
        f"{q2:6d} "
        f"{q3:6d} "
        f"{q4:6d} "
        f"{invalid:8d}"
    )
