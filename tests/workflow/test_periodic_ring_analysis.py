import numpy as np

from workflow.periodic_rings import PeriodicRingAnalysis


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


def make_triangle_system():
    # Three Si atoms and three bridging oxygens form one triangle.
    system = FakeSystem(
        types=[1, 1, 1, 2, 2, 2],
        positions=[
            [1.0, 5.0, 5.0],
            [4.0, 5.0, 5.0],
            [2.5, 7.0, 5.0],
            [2.5, 5.0, 5.0],
            [3.25, 6.0, 5.0],
            [1.75, 6.0, 5.0],
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
    return system, neighbors


def test_periodic_analyzer_finds_contractible_triangle():
    system, neighbors = make_triangle_system()
    analysis = PeriodicRingAnalysis(system, neighbors)

    records = analysis.rings()

    assert len(records) == 1
    ring = records[0]
    assert ring["size"] == 3
    assert ring["classification"] == "contractible"
    assert len(ring["si_cycle"]) == 3
    assert len(ring["oxygen_indices"]) == 3
    assert len(set(ring["oxygen_indices"])) == 3
    assert analysis.distribution("contractible") == {3: 1}
    assert analysis.distribution("winding") == {}


def test_parallel_bridging_oxygen_edges_remain_distinct():
    # A triangle with two different bridging oxygens on edge 0--1.
    # The two oxygen choices create two distinct edge-labelled rings.
    system = FakeSystem(
        types=[1, 1, 1, 2, 2, 2, 2],
        positions=[
            [1.0, 1.0, 1.0],  # Si 0
            [3.0, 1.0, 1.0],  # Si 1
            [2.0, 3.0, 1.0],  # Si 2
            [2.0, 1.0, 1.0],  # O 3: edge 0--1
            [2.0, 1.2, 1.0],  # O 4: parallel edge 0--1
            [2.5, 2.0, 1.0],  # O 5: edge 1--2
            [1.5, 2.0, 1.0],  # O 6: edge 0--2
        ],
        box_lengths=[10.0, 10.0, 10.0],
    )
    neighbors = FakeNeighbors({
        0: [3, 4, 6],
        1: [3, 4, 5],
        2: [5, 6],
        3: [0, 1],
        4: [0, 1],
        5: [1, 2],
        6: [0, 2],
    })

    analysis = PeriodicRingAnalysis(system, neighbors)
    rings = analysis.rings()

    assert len(rings) == 2
    oxygen_paths = {ring["oxygen_indices"] for ring in rings}
    assert len(oxygen_paths) == 2
    assert all(ring["size"] == 3 for ring in rings)

    # The two oxygen-labelled paths remain distinct, but the alternative
    # bridging oxygens on the same Si-Si edge must preserve the ring's
    # periodic classification and net translation.
    assert len({ring["classification"] for ring in rings}) == 1

    winding_numbers = {
        ring["winding_numbers"] for ring in rings
    }
    assert len(winding_numbers) == 1

    np.testing.assert_allclose(
        rings[0]["net_displacement"],
        rings[1]["net_displacement"],
        atol=1e-12,
    )

def test_report_statistics_partition_ring_classifications():
    system, neighbors = make_triangle_system()
    analysis = PeriodicRingAnalysis(system, neighbors)

    total = analysis.statistics()["number"]
    contractible = analysis.statistics("contractible")["number"]
    winding = analysis.statistics("winding")["number"]
    invalid = analysis.statistics("invalid")["number"]

    assert total == contractible + winding + invalid

def make_winding_triangle_system():
    # Three Si atoms form a triangle winding once around periodic x.
    system = FakeSystem(
        types=[1, 1, 1, 2, 2, 2],
        positions=[
            [1.0, 5.0, 5.0],
            [4.0, 5.0, 5.0],
            [7.0, 5.0, 5.0],
            [2.5, 5.0, 5.0],
            [5.5, 5.0, 5.0],
            [9.0, 5.0, 5.0],
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
    return system, neighbors


def test_periodic_analyzer_finds_winding_triangle():
    system, neighbors = make_winding_triangle_system()
    analysis = PeriodicRingAnalysis(system, neighbors)

    records = analysis.rings()

    assert len(records) == 1
    ring = records[0]
    assert ring["size"] == 3
    assert ring["classification"] == "winding"
    np.testing.assert_array_equal(ring["winding_numbers"], [1, 0, 0])
    np.testing.assert_allclose(
        ring["net_displacement"], [10.0, 0.0, 0.0], atol=1e-12
    )
    np.testing.assert_allclose(
        ring["residual"], [0.0, 0.0, 0.0], atol=1e-12
    )
    assert analysis.distribution("contractible") == {}
    assert analysis.distribution("winding") == {3: 1}
    assert analysis.distribution("invalid") == {}


def test_periodic_analyzer_distinguishes_contractible_and_winding_rings():
    winding_system, winding_neighbors = make_winding_triangle_system()
    contractible_system, contractible_neighbors = make_triangle_system()

    winding_analysis = PeriodicRingAnalysis(
        winding_system, winding_neighbors
    )
    contractible_analysis = PeriodicRingAnalysis(
        contractible_system, contractible_neighbors
    )

    winding_ring = winding_analysis.rings()[0]
    contractible_ring = contractible_analysis.rings()[0]

    assert winding_ring["classification"] != contractible_ring["classification"]
    np.testing.assert_array_equal(
        winding_ring["winding_numbers"], [1, 0, 0]
    )
    np.testing.assert_array_equal(
        contractible_ring["winding_numbers"], [0, 0, 0]
    )
