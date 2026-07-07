"""
generator.py

Random Sequential Addition (RSA) generator
for x(CaO)–(1-x)SiO2 glasses.
"""

import numpy as np

from composition import calculate_composition
from configuration import GlassConfiguration
from celllist import NeighborGrid
from geometry import distance_squared
from constants import cutoff

class GlassGenerator:
    """
    Glass generator.
    """

    def __init__(
        self,
        x=0.30,
        formula_units=2400,
        density=1.50,
        seed=12345,
    ):
        """
        Initialize the generator.
        """

        self.seed = seed
        self.rng = np.random.default_rng(seed)

        self.comp = calculate_composition(
            x=x,
            formula_units=formula_units,
            density=density,
        )

        self.box = self.comp["box_length"]

        self.cfg = GlassConfiguration(self.box)

        self.grid = NeighborGrid(
            box_length=self.box,
            cell_size=3.0,
        )

        self.species_lists = {}

        self.make_species_lists()

    def summary(self):
        """
        Print a summary of the system.
        """

        print("=" * 60)
        print("CaO–SiO2 Glass Generator")
        print("=" * 60)

        print(f"x(CaO)         : {self.comp['x']:.2f}")
        print(f"Formula units  : {self.comp['formula_units']}")
        print(f"Si atoms       : {self.comp['Si']}")
        print(f"Ca atoms       : {self.comp['Ca']}")
        print(f"O atoms        : {self.comp['O']}")
        print(f"Total atoms    : {self.comp['natoms']}")
        print(f"Density        : {self.comp['density']:.3f} g/cm³")
        print(f"Box length     : {self.box:.3f} Å")
        print(f"Random seed    : {self.seed}")

    def make_species_lists(self):
        """
        Build and shuffle the atom insertion lists.
        """

        self.species_lists = {
            "Si": ["Si"] * self.comp["Si"],
            "Ca": ["Ca"] * self.comp["Ca"],
            "O":  ["O"]  * self.comp["O"],
        }

        for atoms in self.species_lists.values():
            self.rng.shuffle(atoms)

    def check_position(self, species, position):
        """
        Check whether an atom can be placed at the proposed position.

        Parameters
        ----------
        species : str
        Species to insert ("Si", "Ca", or "O").

        position : numpy.ndarray
        Trial Cartesian coordinates.

        Returns
        -------
        bool
            True if the position is acceptable.
        """

        neighbors = self.grid.nearby_atoms(position)

        for atom_id in neighbors:

            other_species = self.cfg.get_species(atom_id)

            other_position = self.cfg.get_position(atom_id)

            r2 = distance_squared(
                position,
                other_position,
                self.box,

            )

            rcut = cutoff(
                species,
                other_species,
            )

            if r2 < rcut * rcut:
                return False

        return True
