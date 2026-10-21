from collections import Counter
from pathlib import Path

import pytest

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.qn import QnAnalysis
from workflow.topology import TopologyAnalysis


DATA_DIR = Path("simulations/data")

COMPOSITIONS = [
    {
        "x": 0.10,
        "file": "x010_glass_final_300K.data",
        "Si": 2016,
        "O": 4256,
        "Ca": 224,
    },
    {
        "x": 0.20,
        "file": "x020_glass_final_300K.data",
        "Si": 1920,
        "O": 4320,
        "Ca": 480,
    },
    {
        "x": 0.30,
        "file": "glass_final_300K.data",
        "Si": 1680,
        "O": 4080,
        "Ca": 720,
    },
    {
        "x": 0.40,
        "file": "x040_glass_final_300K.data",
        "Si": 1560,
        "O": 4160,
        "Ca": 1040,
    },
]


@pytest.mark.parametrize(
    "case",
    COMPOSITIONS,
    ids=lambda case: f"x={case['x']:.2f}",
)
def test_multicomposition_structural_consistency(case):
    # --------------------------------------------------
    # Load the final 300 K structure
    # --------------------------------------------------
    filename = DATA_DIR / case["file"]

    assert filename.is_file(), f"Missing structure: {filename}"

    system = read_lammps(str(filename))

    # --------------------------------------------------
    # Verify atom counts and composition
    # LAMMPS types: 1=Si, 2=O, 3=Ca
    # --------------------------------------------------
    type_counts = Counter(int(t) for t in system.types)

    assert type_counts[1] == case["Si"]
    assert type_counts[2] == case["O"]
    assert type_counts[3] == case["Ca"]

    expected_atoms = case["Si"] + case["O"] + case["Ca"]
    assert system.natoms == expected_atoms

    ca_fraction = case["Ca"] / (case["Ca"] + case["Si"])
    assert abs(ca_fraction - case["x"]) < 1.0e-12

    # --------------------------------------------------
    # Build the validated species-specific neighbor list
    # --------------------------------------------------
    neighbors = NeighborList(system, BOND_CUTOFFS)
    neighbors.build()

    # Check that neighbors are unique, exclude self, and
    # are symmetric in the periodic neighbor list.
    for i in range(system.natoms):
        adjacent = neighbors.neighbors_of(i)

        assert i not in adjacent
        assert len(adjacent) == len(set(adjacent))

        for j in adjacent:
            assert i in neighbors.neighbors_of(j)

    # --------------------------------------------------
    # Oxygen classification and Qn accounting
    # --------------------------------------------------
    qn = QnAnalysis(system, neighbors)

    oxygen_stats = qn.oxygen_statistics()

    assert oxygen_stats["total"] == case["O"]
    assert (
        oxygen_stats["BO"]
        + oxygen_stats["NBO"]
        + oxygen_stats["other"]
        == case["O"]
    )

    qn_counts = qn.qn_distribution()

    valid_qn_total = sum(qn_counts[n] for n in range(5))

    assert valid_qn_total + qn_counts["invalid"] == case["Si"]

    # Each BO classified by QnAnalysis must also be
    # recognized by the topology analysis.
    topology = TopologyAnalysis(system, neighbors)
    network_stats = topology.network_statistics()

    assert (
        network_stats["bridging_oxygen_count"]
        == oxygen_stats["BO"]
    )
    assert network_stats["number_si"] == case["Si"]

    # --------------------------------------------------
    # Network graph and connected-component accounting
    # --------------------------------------------------
    graph = topology.si_connections()
    components = topology.connected_components()
    component_stats = topology.component_statistics()

    si_indices = {
        i for i, atom_type in enumerate(system.types)
        if atom_type == 1
    }

    assert set(graph) == si_indices
    assert 0 < len(components) <= case["Si"]
    assert sum(len(component) for component in components) == case["Si"]

    # Every component is nonempty, and the reported
    # component statistics agree with the actual graph.
    assert all(len(component) > 0 for component in components)
    assert (
        component_stats["number_components"] == len(components)
    )
    assert (
        component_stats["largest_component"]
        == max(len(component) for component in components)
    )
    assert (
        component_stats["isolated_si"]
        == sum(len(component) == 1 for component in components)
    )
