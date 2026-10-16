import numpy as np

from workflow.geometry import distance
from workflow.neighbors import NeighborList
from workflow.system import AtomicSystem


def brute_force_neighbors(system, cutoffs):
    result = [[] for _ in range(system.natoms)]

    for i in range(system.natoms):
        for j in range(i + 1, system.natoms):
            cutoff = cutoffs.get((system.types[i], system.types[j]))
            if cutoff is None:
                continue

            d = distance(
                system.positions[i],
                system.positions[j],
                system.box_lengths,
            )
            if d <= cutoff:
                result[i].append(j)
                result[j].append(i)

    for neighbors in result:
        neighbors.sort()

    return result


def test_cell_list_matches_brute_force():
    rng = np.random.default_rng(20261015)

    box = np.array([
        [2.0, 22.0],
        [-3.0, 24.0],
        [5.0, 31.0],
    ])
    lengths = box[:, 1] - box[:, 0]
    natoms = 180

    positions = box[:, 0] + rng.random((natoms, 3)) * lengths
    types = rng.choice([1, 2, 3], size=natoms)

    # Include atoms close across periodic boundaries.
    positions[0] = [2.05, -2.0, 10.0]
    positions[1] = [21.95, -2.0, 10.0]
    types[0] = 1
    types[1] = 2

    system = AtomicSystem(
        ids=np.arange(1, natoms + 1),
        types=types,
        charges=np.zeros(natoms),
        positions=positions,
        box=box,
    )

    cutoffs = {
        (1, 2): 2.30,
        (2, 1): 2.30,
        (3, 2): 3.10,
        (2, 3): 3.10,
    }

    expected = brute_force_neighbors(system, cutoffs)
    actual = NeighborList(system, cutoffs)
    actual.build()

    assert actual.list == expected

    # Check symmetry and coordination consistency.
    for i, neighbors in enumerate(actual.list):
        assert len(neighbors) == actual.coordination(i)
        for j in neighbors:
            assert i in actual.neighbors_of(j)
