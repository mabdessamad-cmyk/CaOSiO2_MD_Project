"""
generator.py

Random Sequential Addition (RSA) generator
for x(CaO)–(1-x)SiO2 glasses.
"""

import numpy as np

from composition import calculate_composition
from configuration import GlassConfiguration
from celllist import NeighborGrid


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
