#!/usr/bin/env python3
"""
generate_glass.py

Publication-grade random structure generator for
x(CaO)-(1-x)SiO2 glasses.

Author : Mouaffak + ChatGPT
Version : 1.0
"""

from dataclasses import dataclass, field
from time import perf_counter
import random
import numpy as np

from constants import TYPE, CHARGE
from composition import calculate_composition
from geometry import random_position
from celllist import NeighborGrid


# ============================================================
# Configuration object
# ============================================================

@dataclass
class GlassConfiguration:
    """
    Complete atomic configuration.
    """

    box_length: float

    positions: list = field(default_factory=list)
    species: list = field(default_factory=list)
    atom_types: list = field(default_factory=list)
    charges: list = field(default_factory=list)

    def add_atom(self, species, position):
        """
        Add one atom.
        """

        self.species.append(species)
        self.positions.append(np.asarray(position, dtype=float))
        self.atom_types.append(TYPE[species])
        self.charges.append(CHARGE[species])

    @property
    def natoms(self):
        return len(self.positions)


# ============================================================
# Generator
# ============================================================

class GlassGenerator:

    def __init__(self,
                 x=0.30,
                 total_atoms=6480,
                 density=1.50,
                 seed=12345):

        self.x = x
        self.total_atoms = total_atoms
        self.density = density
        self.seed = seed

        self.rng = random.Random(seed)
        self.nprng = np.random.default_rng(seed)
        # Composition dictionary
        self.comp = calculate_composition(
            x=x,
            total_atoms=total_atoms,
            density=density
        )

        self.box = self.comp["box_length"]

               self.cfg = GlassConfiguration(self.box)

        # Largest minimum cutoff = 3.0 Å
        self.grid = NeighborGrid(
            box_length=self.box,
            cell_size=3.0
        )

        self.start_time = None

        self.accepted = 0
        self.rejected = 0
        self.max_trials = 0

        self.accepted = 0
        self.rejected = 0
        self.max_trials = 0

        print("=" * 60)
        print(" CaO-SiO2 Glass Generator")
        print("=" * 60)

        print()

        print(f"x(CaO)            : {x:.2f}")
        print(f"Total atoms       : {total_atoms}")
        print(f"Density           : {density:.3f} g/cm3")
        print(f"Box length        : {self.box:.3f} Å")

        print()

        print("Composition")

        print("--------------------------")

        print(f"Si : {self.comp['Si']}")
        print(f"Ca : {self.comp['Ca']}")
        print(f"O  : {self.comp['O']}")

        print()

    # --------------------------------------------------------

    def make_species_list(self):
        """
        Create shuffled species list.
        """

        species = (
            ["Si"] * self.comp["Si"] +
            ["Ca"] * self.comp["Ca"] +
            ["O"]  * self.comp["O"]
        )

        random.shuffle(species)

        return species
# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    gen = GlassGenerator(
        x=0.30,
        total_atoms=6480,
        density=1.50,
        seed=12345
    )

    species = gen.make_species_list()

    print(f"Species list created: {len(species)} atoms")

    print("First 20 atoms:")
    print(species[:20])
