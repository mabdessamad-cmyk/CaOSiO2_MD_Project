import numpy as np
import pytest

from workflow.periodic_rings import (
    bridging_oxygen_edges,
    classify_winding,
    directed_edge_displacement,
)


class FakeSystem:
    def __init__(self, types, positions, box_lengths):
        self.types = np.asarray(types, dtype=int)
        self.positions = np.asarray(positions, dtype=float)
        self.box_lengths = np.asarray(box_lengths, dtype=float)


class FakeNeighbors:
    def __init__(self, adjacency):
        self.adjacency = adjacency

    def neighbors_of(self, index):
        return self.adjacency[index]


def test_two_bridging_oxygens_between_same_si_pair_remain_distinct():
    system = FakeSystem(
        types=[1, 2, 2, 1],
        positions=[
            [1.0, 1.0, 1.0],  # Si, index 0
            [1.5, 1.0, 1.0],  # O,  index 1
            [1.5, 1.2, 1.0],  # O,  index 2
            [2.0, 1.0, 1.0],  # Si, index 3
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )
    neighbors = FakeNeighbors({
        0: [1, 2],
        1: [0, 3],
        2: [0, 3],
        3: [1, 2],
    })

    edges = bridging_oxygen_edges(system, neighbors)

    assert len(edges) == 2
    assert [edge["oxygen_index"] for edge in edges] == [1, 2]
    assert all((edge["si_i"], edge["si_j"]) == (0, 3) for edge in edges)
    np.testing.assert_allclose(
        edges[0]["displacement"], [1.0, 0.0, 0.0], atol=1e-12
    )
    np.testing.assert_allclose(
        edges[1]["displacement"], [1.0, 0.0, 0.0], atol=1e-12
    )


def test_edge_displacement_works_in_both_directions():
    edge = {
        "si_i": 2,
        "si_j": 7,
        "displacement": np.array([0.8, -0.2, 0.0]),
    }

    np.testing.assert_allclose(
        directed_edge_displacement(edge, 2, 7), [0.8, -0.2, 0.0]
    )
    np.testing.assert_allclose(
        directed_edge_displacement(edge, 7, 2), [-0.8, 0.2, 0.0]
    )


def test_edge_displacement_rejects_unrelated_si_pair():
    edge = {
        "si_i": 2,
        "si_j": 7,
        "displacement": np.array([0.8, 0.0, 0.0]),
    }

    with pytest.raises(ValueError, match="endpoints"):
        directed_edge_displacement(edge, 2, 9)


def test_boundary_crossing_bond_uses_minimum_image_displacements():
    system = FakeSystem(
        types=[1, 2, 1],
        positions=[
            [9.7, 1.0, 1.0],  # Si near upper x boundary
            [0.1, 1.0, 1.0],  # O across periodic boundary
            [0.5, 1.0, 1.0],  # Si
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )
    neighbors = FakeNeighbors({
        0: [1],
        1: [0, 2],
        2: [1],
    })

    edges = bridging_oxygen_edges(system, neighbors)

    assert len(edges) == 1
    np.testing.assert_allclose(
        edges[0]["displacement"], [0.8, 0.0, 0.0], atol=1e-12
    )


def test_zero_net_translation_is_contractible():
    result = classify_winding(
        [
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -2.0, 0.0],
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )

    assert result["classification"] == "contractible"
    np.testing.assert_array_equal(result["winding_numbers"], [0, 0, 0])
    np.testing.assert_allclose(result["net_displacement"], [0, 0, 0])


def test_one_box_length_net_translation_is_winding():
    result = classify_winding(
        [
            [3.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )

    assert result["classification"] == "winding"
    np.testing.assert_array_equal(result["winding_numbers"], [1, 0, 0])
    np.testing.assert_allclose(result["residual"], [0, 0, 0], atol=1e-12)


def test_non_lattice_net_translation_is_invalid():
    result = classify_winding(
        [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        box_lengths=[10.0, 10.0, 10.0],
    )

    assert result["classification"] == "invalid"


@pytest.mark.parametrize(
    "lengths",
    [
        [10.0, 10.0],
        [10.0, 0.0, 10.0],
        [10.0, float("nan"), 10.0],
    ],
)
def test_invalid_box_lengths_are_rejected(lengths):
    with pytest.raises(ValueError):
        classify_winding([], lengths)


def test_invalid_tolerance_is_rejected():
    with pytest.raises(ValueError, match="tolerance"):
        classify_winding([], [10.0, 10.0, 10.0], tolerance=0.0)


def _cycle_winding_from_coordinates(system, neighbors, oxygen_indices, si_cycle):
    """Build the edges from coordinates and classify one supplied Si cycle."""
    edges = bridging_oxygen_edges(system, neighbors)
    by_oxygen = {edge["oxygen_index"]: edge for edge in edges}

    assert len(by_oxygen) == len(edges)
    assert all(index in by_oxygen for index in oxygen_indices)

    displacements = []
    for k, oxygen_index in enumerate(oxygen_indices):
        start_si = si_cycle[k]
        end_si = si_cycle[(k + 1) % len(si_cycle)]
        displacements.append(
            directed_edge_displacement(
                by_oxygen[oxygen_index], start_si, end_si
            )
        )

    return classify_winding(displacements, system.box_lengths)


def test_complete_cycle_crossing_boundary_can_be_contractible():
    # Three Si atoms form a small triangle that straddles x = 0.
    # Its net displacement must be zero, not one box length.
    system = FakeSystem(
        types=[1, 1, 1, 2, 2, 2],
        positions=[
            [9.5, 5.0, 5.0],  # Si 0
            [0.5, 5.0, 5.0],  # Si 1
            [0.5, 6.0, 5.0],  # Si 2
            [0.0, 5.0, 5.0],  # O 3: Si 0--Si 1
            [0.5, 5.5, 5.0],  # O 4: Si 1--Si 2
            [0.0, 5.5, 5.0],  # O 5: Si 2--Si 0
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )
    neighbors = FakeNeighbors({
        0: [3, 5],
        1: [3, 4],
        2: [4, 5],
        3: [0, 1],
        4: [1, 2],
        5: [0, 2],
    })

    result = _cycle_winding_from_coordinates(
        system, neighbors,
        oxygen_indices=[3, 4, 5],
        si_cycle=[0, 1, 2],
    )

    assert result["classification"] == "contractible"
    np.testing.assert_array_equal(result["winding_numbers"], [0, 0, 0])
    np.testing.assert_allclose(result["net_displacement"], [0, 0, 0], atol=1e-12)


def test_complete_cycle_winds_once_around_periodic_x():
    # The Si positions advance by 3, 3, and 4 Angstrom around x.
    # The final edge crosses the boundary; the total is one box length.
    system = FakeSystem(
        types=[1, 1, 1, 2, 2, 2],
        positions=[
            [1.0, 5.0, 5.0],  # Si 0
            [4.0, 5.0, 5.0],  # Si 1
            [7.0, 5.0, 5.0],  # Si 2
            [2.5, 5.0, 5.0],  # O 3: Si 0--Si 1
            [5.5, 5.0, 5.0],  # O 4: Si 1--Si 2
            [9.0, 5.0, 5.0],  # O 5: Si 2--Si 0 across boundary
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )
    neighbors = FakeNeighbors({
        0: [3, 5],
        1: [3, 4],
        2: [4, 5],
        3: [0, 1],
        4: [1, 2],
        5: [0, 2],
    })

    result = _cycle_winding_from_coordinates(
        system, neighbors,
        oxygen_indices=[3, 4, 5],
        si_cycle=[0, 1, 2],
    )

    assert result["classification"] == "winding"
    np.testing.assert_array_equal(result["winding_numbers"], [1, 0, 0])
    np.testing.assert_allclose(
        result["net_displacement"], [10.0, 0.0, 0.0], atol=1e-12
    )


def test_cycle_winding_validates_edge_identity_and_closing_edge():
    from workflow.periodic_rings import classify_cycle_winding

    edges = [
        {
            "oxygen_index": 10,
            "si_i": 0,
            "si_j": 1,
            "displacement": np.array([1.0, 0.0, 0.0]),
        },
        {
            "oxygen_index": 11,
            "si_i": 1,
            "si_j": 2,
            "displacement": np.array([0.0, 1.0, 0.0]),
        },
        {
            "oxygen_index": 12,
            "si_i": 0,
            "si_j": 2,
            "displacement": np.array([1.0, 1.0, 0.0]),
        },
    ]

    result = classify_cycle_winding(
        si_cycle=[0, 1, 2],
        oxygen_indices=[10, 11, 12],
        edges=edges,
        box_lengths=[10.0, 10.0, 10.0],
    )

    assert result["classification"] == "contractible"
    np.testing.assert_array_equal(result["winding_numbers"], [0, 0, 0])


def test_cycle_winding_rejects_wrong_edge_for_cycle_pair():
    from workflow.periodic_rings import classify_cycle_winding

    edges = [
        {
            "oxygen_index": 10,
            "si_i": 0,
            "si_j": 2,
            "displacement": np.array([1.0, 0.0, 0.0]),
        },
        {
            "oxygen_index": 11,
            "si_i": 1,
            "si_j": 2,
            "displacement": np.array([0.0, 1.0, 0.0]),
        },
        {
            "oxygen_index": 12,
            "si_i": 0,
            "si_j": 1,
            "displacement": np.array([1.0, 1.0, 0.0]),
        },
    ]

    with pytest.raises(ValueError, match="does not connect"):
        classify_cycle_winding(
            si_cycle=[0, 1, 2],
            oxygen_indices=[10, 11, 12],
            edges=edges,
            box_lengths=[10.0, 10.0, 10.0],
        )


@pytest.mark.parametrize(
    "si_cycle, oxygen_indices, message",
    [
        ([0, 1], [10, 11], "at least three"),
        ([0, 1, 0], [10, 11, 12], "must not repeat vertices"),
        ([0, 1, 2], [10, 11], "exactly one"),
        ([0, 1, 2], [10, 10, 12], "must not reuse"),
    ],
)
def test_cycle_winding_rejects_malformed_cycle(
    si_cycle, oxygen_indices, message
):
    from workflow.periodic_rings import classify_cycle_winding

    with pytest.raises(ValueError, match=message):
        classify_cycle_winding(
            si_cycle=si_cycle,
            oxygen_indices=oxygen_indices,
            edges=[],
            box_lengths=[10.0, 10.0, 10.0],
        )


def test_cycle_winding_rejects_missing_oxygen_edge():
    from workflow.periodic_rings import classify_cycle_winding

    with pytest.raises(ValueError, match="No bridging-oxygen edge"):
        classify_cycle_winding(
            si_cycle=[0, 1, 2],
            oxygen_indices=[10, 11, 12],
            edges=[],
            box_lengths=[10.0, 10.0, 10.0],
        )


def test_cycle_winding_rejects_duplicate_edge_records():
    from workflow.periodic_rings import classify_cycle_winding

    edge = {
        "oxygen_index": 10,
        "si_i": 0,
        "si_j": 1,
        "displacement": np.array([1.0, 0.0, 0.0]),
    }

    with pytest.raises(ValueError, match="Duplicate edge records"):
        classify_cycle_winding(
            si_cycle=[0, 1, 2],
            oxygen_indices=[10, 11, 12],
            edges=[edge, edge],
            box_lengths=[10.0, 10.0, 10.0],
        )

@pytest.mark.parametrize(
    "adjacency, message",
    [
        ({0: [1], 1: [0, 0]}, "[Dd]uplicate"),
        ({0: [1], 1: []}, "symmetric"),
        ({0: [99], 1: [0]}, "index"),
        ({0: [0], 1: []}, "[Ss]elf"),
    ],
)
def test_bridging_oxygen_edges_rejects_invalid_neighbor_lists(
    adjacency, message
):
    system = FakeSystem(
        types=[1, 2],
        positions=[
            [1.0, 1.0, 1.0],
            [1.5, 1.0, 1.0],
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )
    neighbors = FakeNeighbors(adjacency)

    with pytest.raises(ValueError, match=message):
        bridging_oxygen_edges(system, neighbors)
