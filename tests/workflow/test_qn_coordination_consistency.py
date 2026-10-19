from collections import Counter

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.qn import QnAnalysis


DATA_FILE = "simulations/data/glass_final_300K.data"


def test_silicon_coordination_and_qn_consistency():
    system = read_lammps(DATA_FILE)
    neighbors = NeighborList(system, BOND_CUTOFFS)
    neighbors.build()
    qn = QnAnalysis(system, neighbors)

    si_indices = [
        i for i, atom_type in enumerate(system.types)
        if atom_type == 1
    ]

    coordination_counts = Counter()
    qn_counts = Counter()

    for i in si_indices:
        oxygen_neighbors = [
            j for j in neighbors.neighbors_of(i)
            if system.types[j] == 2
        ]
        coordination_counts[len(oxygen_neighbors)] += 1
        qn_counts[qn.qn_of_si(i)] += 1

    # Reference composition and Si coordination populations.
    assert len(si_indices) == 1680
    assert coordination_counts == {4: 1675, 5: 5}

    # All-Si Qn values: Q0-Q4 are valid in this distribution;
    # Q5 is outside the defined Q0-Q4 range.
    assert sum(qn_counts[n] for n in range(5)) == 1678
    assert qn_counts[5] == 2

    # Non-tetrahedral Si must be counted by coordination,
    # not inferred solely from whether Qn is in the range 0-4.
    assert coordination_counts[5] == 5
