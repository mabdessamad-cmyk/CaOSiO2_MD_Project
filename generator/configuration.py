"""
configuration.py

Container for a glass configuration.
"""

import numpy as np


class GlassConfiguration:
    """
    Stores atom species and Cartesian coordinates.
    """

    def __init__(self, box_length):

        self.box = float(box_length)

        self.positions = []
        self.species = []

    def __len__(self):
        return len(self.species)

    def add_atom(self, species, position):
        """
        Add one atom.
        """

        self.species.append(species)
        self.positions.append(np.asarray(position, dtype=float))

    @property
    def natoms(self):
        return len(self.species)

    def get_position(self, i):
        return self.positions[i]

    def get_species(self, i):
        return self.species[i]

    def positions_array(self):
        return np.asarray(self.positions)

    def species_list(self):
        return list(self.species)
