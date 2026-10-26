from workflow.rings import RingAnalysis


class MockTopology:
    def __init__(self, graph):
        self._graph = graph

    def si_connections(self):
        return self._graph


def analyze_graph(graph):
    topology = MockTopology(graph)
    return RingAnalysis(None, topology)


def test_triangle_has_one_three_member_ring():
    graph = {
        0: {1, 2},
        1: {0, 2},
        2: {0, 1},
    }

    rings = analyze_graph(graph).primitive_cycles()

    assert len(rings) == 1
    assert sorted(map(len, rings)) == [3]


def test_square_has_one_four_member_ring():
    graph = {
        0: {1, 3},
        1: {0, 2},
        2: {1, 3},
        3: {0, 2},
    }

    rings = analyze_graph(graph).primitive_cycles()

    assert len(rings) == 1
    assert sorted(map(len, rings)) == [4]


def test_theta_graph_finds_all_three_four_member_rings():
    # Three internally disjoint paths connect vertices 0 and 1.
    graph = {
        0: {2, 3, 4},
        1: {2, 3, 4},
        2: {0, 1},
        3: {0, 1},
        4: {0, 1},
    }

    rings = analyze_graph(graph).primitive_cycles()

    assert len(rings) == 3
    assert sorted(map(len, rings)) == [4, 4, 4]


def test_square_with_diagonal_gives_two_triangles():
    graph = {
        0: {1, 2, 3},
        1: {0, 2},
        2: {0, 1, 3},
        3: {0, 2},
    }

    rings = analyze_graph(graph).primitive_cycles()

    assert len(rings) == 2
    assert sorted(map(len, rings)) == [3, 3]


def test_external_shortcut_rejects_six_member_cycle():
    # The external vertex 6 gives a two-edge path from 0 to 3.
    # The six-member cycle has a three-edge arc between them.
    graph = {
        0: {1, 5, 6},
        1: {0, 2},
        2: {1, 3},
        3: {2, 4, 6},
        4: {3, 5},
        5: {4, 0},
        6: {0, 3},
    }

    rings = analyze_graph(graph).primitive_cycles()

    assert not any(set(cycle) == set(range(6)) for cycle in rings)


def test_disconnected_triangle_and_square():
    graph = {
        0: {1, 2},
        1: {0, 2},
        2: {0, 1},
        3: {4, 6},
        4: {3, 5},
        5: {4, 6},
        6: {3, 5},
    }

    rings = analyze_graph(graph).primitive_cycles()

    assert sorted(map(len, rings)) == [3, 4]


def test_ring_size_limit_excludes_larger_cycle():
    # A chordless 13-member cycle must be excluded when the
    # maximum allowed ring size is 12.
    n = 13
    graph = {
        i: {(i - 1) % n, (i + 1) % n}
        for i in range(n)
    }

    rings = RingAnalysis(
        None, MockTopology(graph), max_ring_size=12
    ).primitive_cycles()

    assert rings == set()
