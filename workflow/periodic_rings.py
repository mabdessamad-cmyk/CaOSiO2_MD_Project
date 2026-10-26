"""
Periodic geometry primitives for silicate-ring analysis.

This module is intentionally independent of workflow/rings.py and
workflow/report.py. It preserves the identity of each bridging oxygen
and provides a periodic winding classification for supplied paths.

Si atoms and oxygen atoms are represented by zero-based array indices.
The graph edge orientation is from the lower Si array index to the higher
Si array index. Reverse traversal uses the negative displacement.
"""

import numpy as np

from workflow.geometry import vector


def bridging_oxygen_edges(system, neighbors):
    """
    Return one edge record per oxygen connected to exactly two Si atoms.

    Each record contains:
        oxygen_index : oxygen's zero-based array index
        si_i         : first Si array index
        si_j         : second Si array index
        displacement : periodic displacement from si_i to si_j via O

    Distinct bridging oxygens are retained as distinct edges, even when
    they connect the same pair of Si atoms.
    """
    types = np.asarray(system.types)
    positions = np.asarray(system.positions, dtype=float)
    box_lengths = np.asarray(system.box_lengths, dtype=float)

    if positions.shape != (len(types), 3):
        raise ValueError("system.positions must have shape (n_atoms, 3).")
    if box_lengths.shape != (3,):
        raise ValueError("system.box_lengths must contain three lengths.")

    if not np.all(np.isfinite(box_lengths)) or np.any(box_lengths <= 0):
        raise ValueError("All periodic box lengths must be finite and positive.")

    # Validate the complete neighbor list before constructing graph edges.
    n_atoms = len(types)
    adjacency = []

    for atom_index in range(n_atoms):
        try:
            atom_neighbors = list(neighbors.neighbors_of(atom_index))
        except (IndexError, KeyError, TypeError) as exc:
            raise ValueError(
                f"Could not read neighbors for atom index {atom_index}."
            ) from exc

        seen = set()

        for neighbor_index in atom_neighbors:
            if (
                isinstance(neighbor_index, (bool, np.bool_))
                or not isinstance(neighbor_index, (int, np.integer))
                or not 0 <= neighbor_index < n_atoms
            ):
                raise ValueError(
                    f"Invalid neighbor index {neighbor_index!r} "
                    f"for atom index {atom_index}."
                )

            if neighbor_index == atom_index:
                raise ValueError(
                    f"Self-neighbor detected for atom index {atom_index}."
                )

            if neighbor_index in seen:
                raise ValueError(
                    f"Duplicate neighbor {neighbor_index} "
                    f"for atom index {atom_index}."
                )

            seen.add(neighbor_index)

        adjacency.append(seen)

    for atom_index, atom_neighbors in enumerate(adjacency):
        for neighbor_index in atom_neighbors:
            if atom_index not in adjacency[neighbor_index]:
                raise ValueError(
                    f"Neighbor list is not symmetric: atom {atom_index} "
                    f"lists atom {neighbor_index}, but the reverse "
                    "relationship is missing."
                )

    edges = []

    for oxygen_index, atom_type in enumerate(types):
        if atom_type != 2:
            continue

        si_indices = sorted(
            index
            for index in adjacency[oxygen_index]
            if types[index] == 1
        )

        if len(si_indices) != 2:
            continue

        si_i, si_j = si_indices

        # Follow the local path Si_i -> O -> Si_j using minimum-image
        # vectors for both individual Si-O bonds.
        displacement = (
            vector(
                positions[si_i],
                positions[oxygen_index],
                box_lengths,
            )
            + vector(
                positions[oxygen_index],
                positions[si_j],
                box_lengths,
            )
        )

        edges.append({
            "oxygen_index": int(oxygen_index),
            "si_i": int(si_i),
            "si_j": int(si_j),
            "displacement": np.asarray(displacement, dtype=float),
        })

    return edges


def directed_edge_displacement(edge, start_si, end_si):
    """
    Return an edge displacement in the requested traversal direction.
    """
    si_i = edge["si_i"]
    si_j = edge["si_j"]

    if start_si == si_i and end_si == si_j:
        return np.array(edge["displacement"], dtype=float, copy=True)

    if start_si == si_j and end_si == si_i:
        return -np.array(edge["displacement"], dtype=float, copy=True)

    raise ValueError(
        "Requested Si pair does not match the endpoints of this edge."
    )


def classify_winding(edge_displacements, box_lengths, tolerance=1.0e-5):
    """
    Classify the net periodic translation of a closed path.

    Parameters
    ----------
    edge_displacements : iterable of 3-vectors
        Directed displacement for every edge in traversal order.
    box_lengths : array-like, shape (3,)
        Orthorhombic periodic box lengths.
    tolerance : float
        Maximum absolute component residual in Angstrom.

    Returns
    -------
    dict
        net_displacement : summed path displacement, Angstrom
        winding_numbers  : nearest integer lattice translation
        residual         : difference from that lattice translation
        classification   : 'contractible', 'winding', or 'invalid'

    'invalid' means the net displacement is not sufficiently close to an
    integer lattice translation. It can indicate an incorrectly assembled
    path, inconsistent edge directions, or an unsuitable tolerance.
    This function does not itself verify that the supplied path is a cycle.
    """
    lengths = np.asarray(box_lengths, dtype=float)

    if lengths.shape != (3,):
        raise ValueError("box_lengths must contain exactly three lengths.")
    if not np.all(np.isfinite(lengths)) or np.any(lengths <= 0):
        raise ValueError("Box lengths must be finite and positive.")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive.")

    vectors = [np.asarray(d, dtype=float) for d in edge_displacements]

    if any(d.shape != (3,) for d in vectors):
        raise ValueError("Every edge displacement must be a 3-vector.")
    if any(not np.all(np.isfinite(d)) for d in vectors):
        raise ValueError("Edge displacements must be finite.")

    if vectors:
        net = np.sum(vectors, axis=0)
    else:
        net = np.zeros(3, dtype=float)

    winding = np.rint(net / lengths).astype(int)
    residual = net - winding * lengths

    if np.any(np.abs(residual) > tolerance):
        classification = "invalid"
    elif np.any(winding != 0):
        classification = "winding"
    else:
        classification = "contractible"

    return {
        "net_displacement": net,
        "winding_numbers": winding,
        "residual": residual,
        "classification": classification,
    }


def classify_cycle_winding(
    si_cycle,
    oxygen_indices,
    edges,
    box_lengths,
    tolerance=1.0e-5,
):
    """Validate a Si cycle and classify its periodic winding."""
    si_cycle = list(si_cycle)
    oxygen_indices = list(oxygen_indices)

    if len(si_cycle) < 3:
        raise ValueError("A cycle must contain at least three Si vertices.")

    if any(
        not isinstance(index, (int, np.integer))
        or isinstance(index, (bool, np.bool_))
        for index in si_cycle
    ):
        raise TypeError("Si cycle indices must be integers.")

    if len(set(si_cycle)) != len(si_cycle):
        raise ValueError("A simple Si cycle must not repeat vertices.")

    if len(oxygen_indices) != len(si_cycle):
        raise ValueError(
            "A cycle must specify exactly one bridging oxygen per edge."
        )

    if any(
        not isinstance(index, (int, np.integer))
        or isinstance(index, (bool, np.bool_))
        for index in oxygen_indices
    ):
        raise TypeError("Oxygen indices must be integers.")

    if len(set(oxygen_indices)) != len(oxygen_indices):
        raise ValueError("A cycle must not reuse a bridging oxygen.")

    edge_by_oxygen = {}
    for edge in edges:
        oxygen_index = edge["oxygen_index"]
        if oxygen_index in edge_by_oxygen:
            raise ValueError(
                f"Duplicate edge records for oxygen {oxygen_index}."
            )
        edge_by_oxygen[oxygen_index] = edge

    displacements = []
    nvertices = len(si_cycle)

    for k, oxygen_index in enumerate(oxygen_indices):
        if oxygen_index not in edge_by_oxygen:
            raise ValueError(
                f"No bridging-oxygen edge record for oxygen {oxygen_index}."
            )

        start_si = si_cycle[k]
        end_si = si_cycle[(k + 1) % nvertices]
        edge = edge_by_oxygen[oxygen_index]

        if {edge["si_i"], edge["si_j"]} != {start_si, end_si}:
            raise ValueError(
                f"Oxygen {oxygen_index} does not connect the required "
                f"cycle edge {start_si}-{end_si}."
            )

        displacements.append(
            directed_edge_displacement(edge, start_si, end_si)
        )

    return classify_winding(
        displacements,
        box_lengths,
        tolerance=tolerance,
    )


class PeriodicRingAnalysis:
    """
    Oxygen-labelled periodic multigraph ring analysis.

    Ring candidates are chordless, primitive cycles in the underlying
    simple Si graph. Each distinct bridging oxygen on a cycle edge
    generates a distinct ring path.

    Ring size is the number of Si vertices. Each ring record preserves
    the ordered Si cycle, ordered oxygen-edge identities, winding vector,
    and winding classification.

    This analyzer does not modify the existing RingAnalysis results.
    """

    def __init__(
        self,
        system,
        neighbors,
        max_ring_size=12,
        tolerance=1.0e-5,
    ):
        from workflow.rings import RingAnalysis

        if (
            not isinstance(max_ring_size, (int, np.integer))
            or isinstance(max_ring_size, (bool, np.bool_))
        ):
            raise TypeError("max_ring_size must be an integer.")
        if max_ring_size < 3:
            raise ValueError("max_ring_size must be at least 3.")
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("tolerance must be finite and positive.")

        self.system = system
        self.neighbors = neighbors
        self.max_ring_size = int(max_ring_size)
        self.tolerance = float(tolerance)
        self.edges = bridging_oxygen_edges(system, neighbors)
        self._records_cache = None
        self._ring_analysis_class = RingAnalysis

        self.box_lengths = np.asarray(system.box_lengths, dtype=float)

        # Include all Si atoms, including isolated vertices.
        si_indices = np.where(np.asarray(system.types) == 1)[0]
        self._graph = {int(i): set() for i in si_indices}
        self._edges_by_pair = {}

        for edge in self.edges:
            i = edge["si_i"]
            j = edge["si_j"]

            self._graph[i].add(j)
            self._graph[j].add(i)

            pair = (min(i, j), max(i, j))
            self._edges_by_pair.setdefault(pair, []).append(edge)

        for pair in self._edges_by_pair:
            self._edges_by_pair[pair].sort(
                key=lambda edge: edge["oxygen_index"]
            )

    def graph(self):
        """Return the underlying simple Si graph."""
        return {
            node: set(adjacent)
            for node, adjacent in self._graph.items()
        }

    def _enumerate_records(self):
        ring_tools = self._ring_analysis_class

        # Reuse the existing graph algorithms through a small adapter
        # instance. The adapter supplies the same graph interface without
        # constructing a second atomic neighbor list.
        graph_adapter = object.__new__(ring_tools)

        candidates = graph_adapter._chordless_cycles(
            self._graph,
            self.max_ring_size,
        )

        distance_cache = {}
        for cycle in candidates:
            for vertex in cycle:
                if vertex not in distance_cache:
                    distance_cache[vertex] = (
                        ring_tools._distances_from(
                            self._graph, vertex
                        )
                    )

        primitive = [
            cycle
            for cycle in candidates
            if ring_tools._is_primitive(cycle, distance_cache)
        ]

        records = []

        for cycle in sorted(primitive):
            edge_choices = []

            for k, start_si in enumerate(cycle):
                end_si = cycle[(k + 1) % len(cycle)]
                pair = (min(start_si, end_si), max(start_si, end_si))
                choices = self._edges_by_pair.get(pair, [])

                if not choices:
                    raise RuntimeError(
                        f"Cycle edge {start_si}-{end_si} has no "
                        "bridging-oxygen edge record."
                    )

                edge_choices.append(choices)

            from itertools import product

            for selected_edges in product(*edge_choices):
                oxygen_indices = tuple(
                    edge["oxygen_index"]
                    for edge in selected_edges
                )

                winding = classify_cycle_winding(
                    si_cycle=cycle,
                    oxygen_indices=oxygen_indices,
                    edges=self.edges,
                    box_lengths=self.box_lengths,
                    tolerance=self.tolerance,
                )

                records.append({
                    "si_cycle": tuple(cycle),
                    "oxygen_indices": oxygen_indices,
                    "size": len(cycle),
                    "classification": winding["classification"],
                    "winding_numbers": tuple(
                        int(value)
                        for value in winding["winding_numbers"]
                    ),
                    "net_displacement": np.array(
                        winding["net_displacement"], copy=True
                    ),
                    "residual": np.array(
                        winding["residual"], copy=True
                    ),
                })

        return records

    def rings(self):
        """Return ring records; array-valued fields are copied."""
        if self._records_cache is None:
            self._records_cache = self._enumerate_records()

        result = []
        for record in self._records_cache:
            copied = dict(record)
            copied["net_displacement"] = np.array(
                record["net_displacement"], copy=True
            )
            copied["residual"] = np.array(
                record["residual"], copy=True
            )
            result.append(copied)

        return result

    def distribution(self, classification=None):
        """
        Return ring counts by size.

        classification may be 'contractible', 'winding', or 'invalid'.
        None returns all classifications combined.
        """
        allowed = {None, "contractible", "winding", "invalid"}
        if classification not in allowed:
            raise ValueError(
                "classification must be contractible, winding, invalid, "
                "or None."
            )

        records = [
            record for record in self.rings()
            if classification is None
            or record["classification"] == classification
        ]

        return {
            size: sum(record["size"] == size for record in records)
            for size in sorted({record["size"] for record in records})
        }

    def statistics(self, classification=None):
        """Summarize ring counts and sizes for one winding class or all."""
        distribution = self.distribution(classification)
        sizes = [
            record["size"]
            for record in self.rings()
            if classification is None
            or record["classification"] == classification
        ]

        if not sizes:
            return {
                "number": 0,
                "mean": float("nan"),
                "std": float("nan"),
                "minimum": None,
                "maximum": None,
                "distribution": distribution,
            }

        values = np.asarray(sizes, dtype=int)
        return {
            "number": int(values.size),
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "minimum": int(np.min(values)),
            "maximum": int(np.max(values)),
            "distribution": distribution,
        }

    def report(self):
        """Print ring-size distributions separated by winding class."""
        print()
        print(
            "===== Periodic Silicate Primitive Rings "
            f"(3-{self.max_ring_size} Si) ====="
        )

        for classification in (
            "contractible",
            "winding",
            "invalid",
        ):
            stats = self.statistics(classification)
            print()
            print(f"{classification.capitalize()} rings: {stats['number']}")

            for size, count in stats["distribution"].items():
                print(f"  {size}-member: {count}")
