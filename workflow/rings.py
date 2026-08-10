"""
Ring-topology analysis for the silicate Si-Si network.

The graph contains Si atoms as nodes and Si-O-Si bridging
connections as edges.

Ring size = number of Si atoms in the cycle.
"""

from collections import deque
import numpy as np


SI = 1


class RingAnalysis:

    def __init__(self, system, topology):
        self.system = system
        self.topology = topology

    # --------------------------------------------------------
    # Graph
    # --------------------------------------------------------

    def graph(self):
        """
        Return the validated Si-Si network graph.
        """

        return self.topology.si_connections()

    # --------------------------------------------------------
    # Shortest path
    # --------------------------------------------------------

    @staticmethod
    def shortest_path(graph, start, target, excluded_edge):
        """
        Shortest path between two Si atoms while excluding
        one Si-Si edge.
        """

        queue = deque([(start, [start])])
        visited = {start}

        while queue:

            current, path = queue.popleft()

            for neighbor in graph[current]:

                edge = tuple(
                    sorted((current, neighbor))
                )

                if edge == excluded_edge:
                    continue

                if neighbor == target:
                    return path + [neighbor]

                if neighbor not in visited:

                    visited.add(neighbor)

                    queue.append(
                        (neighbor, path + [neighbor])
                    )

        return None

    # --------------------------------------------------------
    # Cycle detection
    # --------------------------------------------------------

    def primitive_cycles(self):
        """
        Determine primitive shortest cycles.

        Each graph edge is temporarily removed and the
        shortest alternative path between its endpoints
        is searched.

        Duplicate cycles are removed using canonical
        rotation/reversal-independent representations.
        """

        graph = self.graph()

        cycles = set()

        for i in graph:

            for j in graph[i]:

                if i >= j:
                    continue

                edge = (i, j)

                path = self.shortest_path(
                    graph,
                    i,
                    j,
                    edge
                )

                if path is None:
                    continue

                # path includes both endpoints.
                cycle = path

                if len(cycle) < 3:
                    continue

                canonical = self._canonical_cycle(
                    cycle
                )

                cycles.add(canonical)

        return cycles

    # --------------------------------------------------------
    # Canonical cycle representation
    # --------------------------------------------------------

    @staticmethod
    def _canonical_cycle(cycle):
        """
        Canonical representation independent of starting
        atom and traversal direction.
        """

        cycle = list(cycle)

        variants = []

        for seq in (cycle, list(reversed(cycle))):

            for k in range(len(seq)):

                rotated = (
                    seq[k:] + seq[:k]
                )

                variants.append(
                    tuple(rotated)
                )

        return min(variants)

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

        result = {}

        for size in sorted(set(sizes)):

            result[int(size)] = int(
                np.sum(sizes == size)
            )

        return result

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
        print("===== Silicate Ring Topology =====")

        print(
            f"Number of rings : {stats['number']}"
        )

        if stats["number"] == 0:
            return

        print(
            f"Mean ring size  : "
            f"{stats['mean']:.3f} Si"
        )

        print(
            f"Std deviation   : "
            f"{stats['std']:.3f} Si"
        )

        print(
            f"Minimum ring    : "
            f"{stats['minimum']} Si"
        )

        print(
            f"Maximum ring    : "
            f"{stats['maximum']} Si"
        )

        print()
        print("Ring-size distribution")

        for size, count in distribution.items():

            fraction = (
                100.0 * count / stats["number"]
            )

            print(
                f"{size}-member : "
                f"{count:5d} "
                f"({fraction:6.2f} %)"
            )
