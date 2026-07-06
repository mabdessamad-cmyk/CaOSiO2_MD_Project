"""
celllist.py

Simple neighbor-grid (cell list) for periodic cubic boxes.
"""

import numpy as np


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
