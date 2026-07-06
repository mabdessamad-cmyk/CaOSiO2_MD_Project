"""
generator.py

Random Sequential Addition (RSA) generator
for CaO-SiO2 glasses.
"""

import random
import numpy as np

from composition import calculate_composition
from configuration import GlassConfiguration
from geometry import random_position, distance_squared
from constants import cutoff
from celllist import NeighborGrid


class GlassGenerator:
    """
    Publication-grade glass generator.
    """
def make_species_list(self):
    """
    Build and randomly shuffle the atom insertion list.
    """

    species = []

    species.extend(["Si"] * self.comp["Si"])
    species.extend(["O"]  * self.comp["O"])
    species.extend(["Ca"] * self.comp["Ca"])

    self.rng.shuffle(species)

    return species

    def __init__(
        self,
        x=0.30,
        total_atoms=6480,
        density=1.50,
        seed=12345,
    ):

        self.seed = seed
        self.rng = random.Random(seed)

        self.comp = calculate_composition(
            x=x,
            total_atoms=total_atoms,
            density=density,
        )

        self.box = self.comp["box_length"]

        self.cfg = GlassConfiguration(self.box)

        self.grid = NeighborGrid(
            box_length=self.box,
            cell_size=3.0,
        )

        print("=" * 60)
        print("CaO-SiO2 Glass Generator")
        print("=" * 60)

        print(f"x(CaO)      : {x:.2f}")
        print(f"Atoms       : {total_atoms}")
        print(f"Density     : {density:.3f} g/cm3")
        print(f"Box length  : {self.box:.3f} Å")
        print()
