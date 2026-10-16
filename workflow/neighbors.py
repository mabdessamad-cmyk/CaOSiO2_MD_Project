"""
Periodic cell-list construction for atomic neighbor lists.

Preserves the NeighborList interface, pair-specific cutoffs,
minimum-image distance calculation, and symmetric neighbor lists.
"""

from itertools import product

import numpy as np

from workflow.geometry import distance


class NeighborList:
    """Neighbor list for an orthorhombic AtomicSystem."""

    def __init__(self, system, cutoffs):
        self.system = system
        self.cutoffs = cutoffs
        self.list = []

    def build(self):
        """Construct the neighbor list using periodic spatial cells."""
        natoms = self.system.natoms
        self.list = [[] for _ in range(natoms)]

        print(f"Building neighbor list ({natoms} atoms)...")

        if natoms < 2 or not self.cutoffs:
            print("Neighbor list completed.")
            return

        box = np.asarray(self.system.box, dtype=float)
        lengths = np.asarray(self.system.box_lengths, dtype=float)
        lower = box[:, 0]

        if np.any(lengths <= 0):
            raise ValueError("All simulation box lengths must be positive.")

        max_cutoff = max(float(c) for c in self.cutoffs.values())
        if max_cutoff <= 0:
            print("Neighbor list completed.")
            return

        # Each cell width is at least the largest cutoff.
        # Thus, searching the adjacent cells in each direction is sufficient.
        ncell = np.maximum(1, np.floor(lengths / max_cutoff).astype(int))
        cell_width = lengths / ncell

        positions = np.asarray(self.system.positions, dtype=float)
        types = np.asarray(self.system.types)

        cell_indices = np.floor(
            (positions - lower) / cell_width
        ).astype(int) % ncell

        cells = {}
        for atom_index, cell_index in enumerate(cell_indices):
            key = tuple(cell_index)
            cells.setdefault(key, []).append(atom_index)

        offsets = list(product((-1, 0, 1), repeat=3))

        for i in range(natoms):
            ci = cell_indices[i]
            candidate_cells = set()

            for offset in offsets:
                key = tuple((ci + offset) % ncell)
                candidate_cells.add(key)

            # Deduplicate cells when an axis has only one or two cells.
            candidates = set()
            for key in candidate_cells:
                candidates.update(cells.get(key, ()))

            for j in sorted(candidates):
                if j <= i:
                    continue

                cutoff = self.cutoffs.get((types[i], types[j]))
                if cutoff is None:
                    continue

                d = distance(positions[i], positions[j], lengths)
                if d <= cutoff:
                    self.list[i].append(j)
                    self.list[j].append(i)

        # Preserve the neighbor ordering of the original brute-force loop.
        for neighbors in self.list:
            neighbors.sort()

        print("Neighbor list completed.")

    def neighbors_of(self, atom):
        return self.list[atom]

    def coordination(self, atom):
        return len(self.list[atom])
