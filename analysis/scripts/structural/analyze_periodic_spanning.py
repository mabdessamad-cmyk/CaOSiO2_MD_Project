#!/usr/bin/env python3
"""
Periodic spanning analysis of the BO-connected Si network.

Uses the existing BO definition and periodic geometry implementation.
Does not enumerate rings and does not modify existing project files.

A component winds along a periodic direction if a closed graph path
has a non-zero integer lattice translation along that direction.
"""

from collections import deque
from pathlib import Path

import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.periodic_rings import bridging_oxygen_edges


STRUCTURES = {
    0.10: "simulations/data/x010_glass_final_300K.data",
    0.20: "simulations/data/x020_glass_final_300K.data",
    0.30: "simulations/data/glass_final_300K.data",
    0.40: "simulations/data/x040_glass_final_300K.data",
}

EXPECTED_COMPONENTS = {
    0.10: (1, 2016, 0),
    0.20: (1, 1920, 0),
    0.30: (5, 1676, 4),
    0.40: (15, 1542, 12),
}

AXES = ("x", "y", "z")
INTEGER_TOLERANCE = 1.0e-6


def edge_translation(edge, positions, box_lengths):
    """Return integer image translation for the directed edge si_i -> si_j."""
    i = edge["si_i"]
    j = edge["si_j"]

    raw = positions[j] - positions[i]
    displacement = np.asarray(edge["displacement"], dtype=float)

    translation_float = (displacement - raw) / box_lengths
    translation = np.rint(translation_float).astype(int)

    residual = translation_float - translation

    if np.any(np.abs(residual) > INTEGER_TOLERANCE):
        raise RuntimeError(
            "Non-integer periodic translation for edge "
            f"{i}->{j}: {translation_float}; residual={residual}"
        )

    return translation


def analyze_component(nodes, adjacency, positions, box_lengths):
    """
    Assign image coordinates by BFS and identify cycle winding vectors.

    Each adjacency record is (neighbor_si, edge_translation), where
    translation is directed from the current Si to the neighbor.
    """
    root = min(nodes)
    image = {root: np.zeros(3, dtype=int)}
    queue = deque([root])

    winding_vectors = set()
    max_closure_residual = 0

    while queue:
        current = queue.popleft()

        for neighbor, translation in adjacency[current]:
            proposed = image[current] + translation

            if neighbor not in image:
                image[neighbor] = proposed
                queue.append(neighbor)
            else:
                mismatch = proposed - image[neighbor]

                if np.any(mismatch != 0):
                    winding_vectors.add(tuple(int(v) for v in mismatch))

    if len(image) != len(nodes):
        raise RuntimeError(
            "Internal BFS error: component traversal did not visit "
            "all vertices."
        )

    if not winding_vectors:
        return {
            "size": len(nodes),
            "winding_vectors": [],
            "winding_axes": (False, False, False),
            "rank": 0,
        }

    vectors = np.asarray(sorted(winding_vectors), dtype=int)
    axes = tuple(bool(np.any(vectors[:, k] != 0)) for k in range(3))
    rank = int(np.linalg.matrix_rank(vectors.astype(float)))

    return {
        "size": len(nodes),
        "winding_vectors": sorted(winding_vectors),
        "winding_axes": axes,
        "rank": rank,
    }


def analyze_structure(x, filename):
    print("\n" + "=" * 76)
    print(f"CaO fraction x = {x:.2f}")
    print(f"Structure: {filename}")
    print("=" * 76)

    system = read_lammps(filename)
    neighbors = NeighborList(system, BOND_CUTOFFS)
    neighbors.build()

    edges = bridging_oxygen_edges(system, neighbors)

    si_indices = [int(i) for i in system.silicon_indices()]
    positions = np.asarray(system.positions, dtype=float)
    box_lengths = np.asarray(system.box_lengths, dtype=float)

    adjacency = {i: [] for i in si_indices}

    for edge in edges:
        i = edge["si_i"]
        j = edge["si_j"]

        translation = edge_translation(
            edge, positions, box_lengths
        )

        adjacency[i].append((j, translation))
        adjacency[j].append((i, -translation))

    # Identify ordinary connected components.
    unseen = set(si_indices)
    components = []

    while unseen:
        root = min(unseen)
        queue = deque([root])
        unseen.remove(root)
        component = []

        while queue:
            current = queue.popleft()
            component.append(current)

            for neighbor, _ in adjacency[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)

        components.append(component)

    components.sort(key=len, reverse=True)

    component_results = [
        analyze_component(
            component, adjacency, positions, box_lengths
        )
        for component in components
    ]

    n_si = len(si_indices)
    largest_size = len(components[0])
    isolated = sum(len(component) == 1 for component in components)

    expected = EXPECTED_COMPONENTS[x]
    observed = (len(components), largest_size, isolated)

    print("\nCONNECTED-COMPONENT VALIDATION")
    print(f"  Number of components : {len(components)}")
    print(f"  Largest component    : {largest_size}")
    print(f"  Isolated Si atoms    : {isolated}")
    print(f"  Expected values      : {expected}")

    if observed != expected:
        raise RuntimeError(
            f"Component validation failed for x={x:.2f}: "
            f"observed {observed}, expected {expected}"
        )

    print("  Validation            : PASS")

    # Summarize periodic winding.
    print("\nPERIODIC WINDING RESULTS")

    for axis_index, axis in enumerate(AXES):
        winding_components = [
            result for result in component_results
            if result["winding_axes"][axis_index]
        ]

        winding_atoms = sum(
            result["size"] for result in winding_components
        )

        print(
            f"  {axis}-direction: "
            f"{len(winding_components)} components; "
            f"{winding_atoms}/{n_si} Si atoms "
            f"({100.0 * winding_atoms / n_si:.4f}%)"
        )

    any_winding = [
        result for result in component_results
        if any(result["winding_axes"])
    ]

    full_3d = [
        result for result in component_results
        if result["rank"] == 3
    ]

    any_winding_atoms = sum(result["size"] for result in any_winding)
    full_3d_atoms = sum(result["size"] for result in full_3d)

    print(
        f"\n  Components winding in >=1 direction: "
        f"{len(any_winding)}"
    )
    print(
        f"  Si in components winding in >=1 direction: "
        f"{any_winding_atoms}/{n_si} "
        f"({100.0 * any_winding_atoms / n_si:.4f}%)"
    )
    print(
        f"  Components with full 3D winding (rank 3): "
        f"{len(full_3d)}"
    )
    print(
        f"  Si in components with full 3D winding: "
        f"{full_3d_atoms}/{n_si} "
        f"({100.0 * full_3d_atoms / n_si:.4f}%)"
    )

    print("\nCOMPONENT-BY-COMPONENT DETAILS")
    print("  component | Si count | winding x/y/z | winding rank")

    for index, result in enumerate(component_results, start=1):
        axes = "".join(
            axis.upper() if result["winding_axes"][k] else "-"
            for k, axis in enumerate(AXES)
        )

        print(
            f"  {index:9d} | "
            f"{result['size']:8d} | "
            f"{axes:13s} | "
            f"{result['rank']:12d}"
        )

    print("\n  Number of BO-labelled Si-Si edges:", len(edges))
    print("  Orthorhombic box lengths (Å):", box_lengths)

    return component_results


def main():
    for x, filename in STRUCTURES.items():
        if not Path(filename).is_file():
            raise FileNotFoundError(filename)

        analyze_structure(x, filename)

    print("\n" + "=" * 76)
    print("SUCCESS: periodic-spanning analysis completed.")
    print("All four connected-component validations passed.")
    print("No ring enumeration was performed.")
    print("No existing project files were modified.")
    print("=" * 76)


if __name__ == "__main__":
    main()
