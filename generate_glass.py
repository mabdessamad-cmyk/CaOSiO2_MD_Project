#!/usr/bin/env python3
"""
generate_glass.py

Publication-grade random structure generator for
x(CaO)-(1-x)SiO2 glasses.

Author : Mouaffak + ChatGPT
Version : 1.0
"""

from writer import write_lammps_data
from dataclasses import dataclass, field
from time import perf_counter
import random
import numpy as np

from constants import TYPE, CHARGE
from composition import calculate_composition
from geometry import random_position, distance
from celllist import NeighborGrid
from validate import validate_composition, validate_distances
from topology import (
    find_si_neighbors,
    connect_tetrahedra,
    validate_coordination,
    complete_tetrahedra,
    merge_oxygens,
    classify_oxygens,
    check_neighbor_symmetry,
)

SIO_BOND = 1.62  # Å

TETRA = np.array([
    [ 1,  1,  1],
    [ 1, -1, -1],
    [-1,  1, -1],
    [-1, -1,  1],
], dtype=float)

TETRA /= np.linalg.norm(TETRA[0])

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

    def get_species(self, i):
        return self.species[i]

    def get_position(self, i):
        return self.positions[i]
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

        formula_units = round(total_atoms / (3.0 - x))

        self.comp = calculate_composition(
        x=x,
        formula_units=formula_units,
        density=density
        )
        self.box = self.comp["box_length"]

        self.cfg = GlassConfiguration(self.box) 

# -----------------------------
# Network data
# -----------------------------

        self.si_neighbors = {}
        self.bridging_oxygens = []
        self.nonbridging_oxygens = []

        # Largest minimum cutoff = 3.0 Å
        self.grid = NeighborGrid(
            box_length=self.box,
            cell_size=3.0
        )

        self.start_time = None
        self.min_dist = {
            ("Si", "Si"): 3.00,
            ("Si", "O"):  1.45,
            ("Si", "Ca"): 2.80,
            ("O",  "O"):  2.20,
            ("O",  "Ca"): 2.00,
            ("Ca", "Ca"): 2.80,
        }

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

    def pair_key(self, species1, species2):
        """
        Return the canonical key for a pair of species.
        """

        order = {"Si": 0, "O": 1, "Ca": 2}

        if order[species1] <= order[species2]:
            return (species1, species2)

        return (species2, species1)

    def tetrahedron(self, si_pos):
        """
        Return the four oxygen positions around one Si atom.
        """

        oxygens = []

        for direction in TETRA:
            pos = si_pos + SIO_BOND * direction
            pos %= self.box
            oxygens.append(pos)

        return oxygens

    def place_atoms(self):
        """
        Place only the silicon atoms.
        """

        n_si = self.comp["Si"]

        nside = int(np.ceil(n_si ** (1.0 / 3.0)))
        spacing = self.box / nside

        positions = []

        for ix in range(nside):
            for iy in range(nside):
                for iz in range(nside):

                    if len(positions) >= n_si:
                        break

                    pos = np.array([
                        (ix + 0.5) * spacing,
                        (iy + 0.5) * spacing,
                        (iz + 0.5) * spacing
                    ])

                    pos += self.nprng.uniform(-0.25, 0.25, 3)
                    pos %= self.box

                    positions.append(pos)

                if len(positions) >= n_si:
                    break

            if len(positions) >= n_si:
                break

        self.nprng.shuffle(positions)

        for i, pos in enumerate(positions, start=1):

            if self.can_place("Si", pos):
                self.cfg.add_atom("Si", pos)
            else:
                raise RuntimeError("Unable to place Si atom.")

            if i % 500 == 0:
                print(f"{i}/{n_si} atoms placed")

        self.si_indices = list(range(self.cfg.natoms))

        # Initialize the Si neighbor dictionary
        for i in range(self.comp["Si"]):
            self.si_neighbors[i] = []

        print()
        print("All Si atoms successfully placed.")

    def build_si_network(self):
        """
        Build a realistic Si neighbor graph.
        """

        print("Building Si network...")

        cutoff = 5.0
        max_neighbors = 4

        nsi = len(self.si_indices)

        # Initialize the neighbor dictionary
        for i in range(nsi):
            self.si_neighbors[i] = []

        si_positions = [
            self.cfg.get_position(i)
            for i in self.si_indices
        ]

        neighbors = find_si_neighbors(
            si_positions,
            self.box,
            cutoff=5.0
        )

        self.si_neighbors = neighbors

        for i in range(10):
            print(f"Si {i}: {neighbors[i]}")

        bad = check_neighbor_symmetry(neighbors)
        print(f"Asymmetric links : {bad}")

        counts = [len(n) for n in neighbors]

        print()
        print("Neighbor statistics")
        print("--------------------------")
        print(f"Average neighbors : {sum(counts)/len(counts):.2f}")
        print(f"Minimum neighbors : {min(counts)}")
        print(f"Maximum neighbors : {max(counts)}")

        target_nbo = 2 * self.comp["Ca"]
        target_bo = self.comp["O"] - target_nbo

        oxygen_positions, oxygen_owners = connect_tetrahedra(
            si_positions,
            neighbors,
            self.box,
            target_bo=target_bo
        )

        coordination = validate_coordination(
        oxygen_owners,
        len(si_positions)
        )

        print()
        print("Coordination before completion")
        print("-----------------------------")
        print(f"Average : {sum(coordination)/len(coordination):.2f}")
        print(f"Minimum : {min(coordination)}")
        print(f"Maximum : {max(coordination)}")

        oxygen_positions, oxygen_owners = complete_tetrahedra(
        si_positions,
        oxygen_positions,
        oxygen_owners,
        coordination
        )

        print()
        print("Before merging")
        print("--------------------------")
        print(f"Total oxygens : {len(oxygen_positions)}")

        bo_before = sum(1 for x in oxygen_owners if len(x) == 2)
        nbo_before = sum(1 for x in oxygen_owners if len(x) == 1)

        print(f"BO  : {bo_before}")
        print(f"NBO : {nbo_before}")

        coordination = validate_coordination(
        oxygen_owners,
        len(si_positions)
        )

        print()
        print("Coordination after completion")
        print("----------------------------")
        print(f"Average : {sum(coordination)/len(coordination):.2f}")
        print(f"Minimum : {min(coordination)}")
        print(f"Maximum : {max(coordination)}")

        merged_oxygens = oxygen_positions
        owners = oxygen_owners

        self.oxygen_positions = merged_oxygens
        self.oxygen_owners = owners

        bo, nbo = classify_oxygens(owners)

        print()
        print("Topology statistics")
        print("--------------------------")
        print(f"Temporary oxygens : {len(oxygen_positions)}")
        print(f"Merged oxygens    : {len(merged_oxygens)}")
        print(f"BO                : {bo}")
        print(f"NBO               : {nbo}")

    def analyze_si_o_distances(self):

        from geometry import distance

        bo = []
        nbo = []

        for k, owners in enumerate(self.oxygen_owners):

            o = self.oxygen_positions[k]

            if len(owners) == 2:

                for si in owners:
                    d = distance(
                        self.cfg.get_position(self.si_indices[si]),
                        o,
                        self.box
                    )
                    bo.append(d)

            else:

                si = owners[0]

                d = distance(
                    self.cfg.get_position(self.si_indices[si]),
                    o,
                    self.box
                )
                nbo.append(d)

        print()
        print("Si-O distance analysis")
        print("--------------------------")
        print(f"BO : min={min(bo):.3f}  avg={sum(bo)/len(bo):.3f}")
        print(f"NBO: min={min(nbo):.3f}  avg={sum(nbo)/len(nbo):.3f}")

    # ==========================================================
    # Validation
    # ==========================================================

    def place_oxygen(self):
        """
        Generate the oxygen atoms.
        """

        print("Generating oxygen atoms...")

        # To be implemented
        target = self.comp["O"]

        for i in range(self.comp["O"]):

            trials = 0

            while True:

                trials += 1

                if trials > 10000:
                    raise RuntimeError("Unable to place O atom")

                pos = random_position(self.box, self.nprng)

                if self.can_place("O", pos):

                    self.cfg.add_atom("O", pos)
                    break

            if (i + 1) % 500 == 0:
                print(f"{i + 1}/{target} oxygen atoms")

    def place_calcium(self):
        """
        Generate the calcium atoms.
        """

        print("Generating calcium atoms...")

        for i in range(self.comp["Ca"]):

            trials = 0

            while True:

                trials += 1

                if trials > 10000:
                    raise RuntimeError("Unable to place Ca atom")

                pos = random_position(self.box, self.nprng)

                if self.can_place("Ca", pos):

                    self.cfg.add_atom("Ca", pos)
                    break

            if (i + 1) % 100 == 0:
                print(f"{i+1}/{self.comp['Ca']} calcium atoms")

    def generate(self):
        """
        Build the complete configuration.
        """

        self.place_atoms()
        self.build_si_network()
        self.place_oxygen()
        self.analyze_si_o_distances()
        self.place_calcium()

        return self.cfg

    def can_place(self, species, position):
        """
        Return True if the atom can be placed without violating
        any minimum-distance constraint.
        """

        for i in range(self.cfg.natoms):

            other_species = self.cfg.get_species(i)
            other_pos = self.cfg.get_position(i)

            pair = self.pair_key(species, other_species)
            cutoff = self.min_dist[pair]

            if distance(position, other_pos, self.box) < cutoff:
                return False

        return True
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

    cfg = gen.generate()

    write_lammps_data(cfg, "output/glass_6480.data")
    validate_composition(cfg)
    validate_distances(cfg)

    print()
    print("LAMMPS data file written to output/glass_6480.data")
    print()
    print("=" * 60)
    print(f"Generated {cfg.natoms} atoms")
    print("=" * 60)
