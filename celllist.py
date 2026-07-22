"""
celllist.py

Simple neighbor-grid (cell list) for periodic cubic boxes.
"""

import numpy as np

from geometry import distance_squared

class NeighborGrid:
    """
    Spatial grid for fast neighbor searching.
    """

    def __init__(self, box_length, cell_size):

        self.box = float(box_length)
        self.cell_size = float(cell_size)

        self.ncell = max(1, int(np.floor(self.box / self.cell_size)))
        self.cell_size = self.box / self.ncell

        self.grid = {}

    def clear(self):
        self.grid.clear()

    def cell_index(self, position):
        """
        Return integer cell coordinates.
        """

        c = np.floor(position / self.cell_size).astype(int)

        c %= self.ncell

        return tuple(c)

    def add_atom(self, atom_id, position):
        """
        Store one atom in the grid.
        """

        cell = self.cell_index(position)

        if cell not in self.grid:
            self.grid[cell] = []

        self.grid[cell].append(atom_id)

    def nearby_atoms(self, position):
        """
        Return indices of atoms in the 27 neighboring cells.
        """

        cx, cy, cz = self.cell_index(position)

        neigh = []

        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):

                    key = (
                        (cx + dx) % self.ncell,
                        (cy + dy) % self.ncell,
                        (cz + dz) % self.ncell,
                    )

                    if key in self.grid:
                        neigh.extend(self.grid[key])
        return neigh

    def can_insert(self, species, position, cutoff=3.0):
        """
        Return True if the new atom is at least cutoff Å away
        from all nearby atoms.
        """

        cutoff2 = cutoff * cutoff

        for atom_id in self.nearby_atoms(position):

            other = self.positions[atom_id]

            if distance_squared(position, other, self.box) < cutoff2:
                return False

        return True


    def insert(self, species, position):
        """
        Insert a newly accepted atom.
        """

        if not hasattr(self, "positions"):
            self.positions = []

        atom_id = len(self.positions)

        self.positions.append(position.copy())

        self.add_atom(atom_id, position)


