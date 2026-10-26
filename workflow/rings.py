"""
Primitive-ring analysis for the silicate Si-Si network.

Graph:
    Si atoms are vertices.
    Si-Si edges represent shared bridging oxygens.

Ring size:
    Number of Si vertices in the cycle.

Convention:
    Enumerate chordless cycles, then retain cycles for which
    the shortest graph distance between every pair of ring
    vertices equals their shorter distance along the ring.

This is a graph-topological analysis; periodic image information
is not explicitly represented in the Si-Si graph.
"""

from collections import deque
import numpy as np


SI = 1


class RingAnalysis:

    def __init__(self, system, topology, max_ring_size=12):
        if not isinstance(max_ring_size, int) or isinstance(max_ring_size, bool):
            raise TypeError("max_ring_size must be an integer.")
        if max_ring_size < 3:
            raise ValueError("max_ring_size must be at least 3.")

        self.system = system
        self.topology = topology
        self.max_ring_size = max_ring_size
        self._cycles_cache = None

    # --------------------------------------------------------
    # Graph
    # --------------------------------------------------------

    def graph(self):
        """Return and validate the undirected Si-Si network."""
        source = self.topology.si_connections()
        graph = {
            int(node): set(neighbors)
            for node, neighbors in source.items()
        }

        for node, neighbors in graph.items():
            for neighbor in neighbors:
                if neighbor not in graph:
                    raise ValueError(
                        f"Graph edge {node}-{neighbor} references "
                        "a missing vertex."
                    )
                if node not in graph[neighbor]:
                    raise ValueError(
                        f"Graph is not symmetric at edge "
                        f"{node}-{neighbor}."
                    )

        return graph

    # --------------------------------------------------------
    # Shortest paths
    # --------------------------------------------------------

    @staticmethod
    def _distances_from(graph, start):
        """Unweighted shortest-path distances from one vertex."""
        distances = {start: 0}
        queue = deque([start])

        while queue:
            current = queue.popleft()

            for neighbor in graph[current]:
                if neighbor not in distances:
                    distances[neighbor] = distances[current] + 1
                    queue.append(neighbor)

        return distances

    # --------------------------------------------------------
    # Canonical cycle representation
    # --------------------------------------------------------

    @staticmethod
    def _canonical_cycle(cycle):
        """Normalize rotation and traversal direction."""
        cycle = list(cycle)

        variants = []
        for sequence in (cycle, list(reversed(cycle))):
            for k in range(len(sequence)):
                variants.append(
                    tuple(sequence[k:] + sequence[:k])
                )

        return min(variants)

    # --------------------------------------------------------
    # Chordless-cycle enumeration
    # --------------------------------------------------------

    def _chordless_cycles(self, graph, max_size):
        """
        Enumerate chordless simple cycles up to max_size vertices.

        The smallest vertex in each cycle is used as its start
        vertex, avoiding rotational duplicates. Canonicalization
        removes duplicates caused by reverse traversal.
        """
        cycles = set()

        for start in sorted(graph):
            for first in sorted(graph[start]):
                if first <= start:
                    continue

                stack = [[start, first]]

                while stack:
                    path = stack.pop()
                    current = path[-1]

                    for neighbor in sorted(graph[current]):

                        # Close a cycle at its starting vertex.
                        if neighbor == start:
                            if 3 <= len(path) <= max_size:
                                cycles.add(
                                    self._canonical_cycle(path)
                                )
                            continue

                        # Do not extend paths beyond the size limit.
                        if len(path) >= max_size:
                            continue

                        # Keep start as the smallest cycle vertex.
                        if neighbor <= start or neighbor in path:
                            continue

                        # A new vertex must not connect to an
                        # earlier path vertex, except start when
                        # it closes the cycle.
                        earlier = path[:-1]
                        linked_earlier = [
                            vertex
                            for vertex in earlier
                            if vertex in graph[neighbor]
                        ]

                        if not linked_earlier:
                            stack.append(path + [neighbor])

                        elif linked_earlier == [start]:
                            cycles.add(
                                self._canonical_cycle(
                                    path + [neighbor]
                                )
                            )

        return cycles

    # --------------------------------------------------------
    # Primitive-ring criterion
    # --------------------------------------------------------

    @staticmethod
    def _is_primitive(cycle, distance_cache):
        """
        Check the shortest-path criterion.

        For every pair of cycle vertices, the shortest distance
        in the full graph must equal the shorter arc distance
        along the cycle. A strictly shorter external path rejects
        the candidate.
        """
        n = len(cycle)

        for i in range(n):
            start = cycle[i]

            if start not in distance_cache:
                # The cache is populated by primitive_cycles().
                raise RuntimeError(
                    "Missing shortest-path distances for a cycle vertex."
                )

            distances = distance_cache[start]

            for j in range(i + 1, n):
                target = cycle[j]
                ring_distance = min(j - i, n - (j - i))
                graph_distance = distances.get(target)

                if graph_distance is None:
                    return False

                if graph_distance < ring_distance:
                    return False

        return True

    # --------------------------------------------------------
    # Ring detection
    # --------------------------------------------------------

    def primitive_cycles(self):
        """
        Return unique size-limited primitive rings.

        Only rings containing 3 through max_ring_size Si vertices
        are considered. This is not exhaustive for larger rings.
        """
        if self._cycles_cache is not None:
            return set(self._cycles_cache)

        graph = self.graph()
        candidates = self._chordless_cycles(
            graph, self.max_ring_size
        )

        # Cache distances only for vertices appearing in candidates.
        distance_cache = {}
        for cycle in candidates:
            for vertex in cycle:
                if vertex not in distance_cache:
                    distance_cache[vertex] = (
                        self._distances_from(graph, vertex)
                    )

        accepted = {
            cycle
            for cycle in candidates
            if self._is_primitive(cycle, distance_cache)
        }
        self._cycles_cache = accepted
        return set(accepted)

    # --------------------------------------------------------
    # Ring sizes
    # --------------------------------------------------------

    def ring_sizes(self):
        cycles = self.primitive_cycles()

        return np.asarray(
            [len(cycle) for cycle in cycles],
            dtype=int
        )

    # --------------------------------------------------------
    # Distribution
    # --------------------------------------------------------

    def distribution(self):
        sizes = self.ring_sizes()

        if len(sizes) == 0:
            return {}

        return {
            int(size): int(np.sum(sizes == size))
            for size in sorted(set(sizes))
        }

    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    def statistics(self):
        sizes = self.ring_sizes()

        if len(sizes) == 0:
            return {
                "number": 0,
                "mean": np.nan,
                "std": np.nan,
                "minimum": None,
                "maximum": None,
            }

        return {
            "number": int(len(sizes)),
            "mean": float(np.mean(sizes)),
            "std": float(np.std(sizes)),
            "minimum": int(np.min(sizes)),
            "maximum": int(np.max(sizes)),
        }

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    def report(self):
        stats = self.statistics()
        distribution = self.distribution()

        print()
        print(
            "===== Silicate Size-Limited Primitive-Ring Topology "
            f"(3-{self.max_ring_size} Si) ====="
        )
        print(f"Number of rings : {stats['number']}")

        if stats["number"] == 0:
            return

        print(f"Mean ring size  : {stats['mean']:.3f} Si")
        print(f"Std deviation   : {stats['std']:.3f} Si")
        print(f"Minimum ring    : {stats['minimum']} Si")
        print(f"Maximum ring    : {stats['maximum']} Si")
        print()
        print("Ring-size distribution")

        for size, count in distribution.items():
            fraction = 100.0 * count / stats["number"]
            print(
                f"{size}-member : {count:5d} "
                f"({fraction:6.2f} %)"
            )
