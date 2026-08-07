from dataclasses import dataclass
import numpy as np


@dataclass
class AtomicSystem:
    """
    Atomic configuration loaded from a LAMMPS data file.
    """

    ids: np.ndarray
    types: np.ndarray
    charges: np.ndarray
    positions: np.ndarray
    box: np.ndarray

    # ---------------------------------------------------------
    # Basic properties
    # ---------------------------------------------------------

    @property
    def natoms(self):
        return len(self.ids)

    @property
    def box_lengths(self):
        return self.box[:, 1] - self.box[:, 0]

    # ---------------------------------------------------------
    # Chemical information
    # ---------------------------------------------------------

    def atom_name(self, atom_type):
        """
        Return the chemical symbol corresponding to a LAMMPS atom type.
        """
        names = {
            1: "Si",
            2: "O",
            3: "Ca",
        }
        return names.get(atom_type, "Unknown")

    # ---------------------------------------------------------
    # Atom selections
    # ---------------------------------------------------------


    def silicon_indices(self):
        return np.where(self.types == 1)[0]


    def oxygen_indices(self):
        return np.where(self.types == 2)[0]


    def calcium_indices(self):
        return np.where(self.types == 3)[0]

    # ---------------------------------------------------------
    # Statistics
    # ---------------------------------------------------------


    def number_of(self, atom_type):
        """
        Return the number of atoms of a given LAMMPS atom type.
        """
        return int(np.sum(self.types == atom_type))

    # ---------------------------------------------------------
    # Display
    # ---------------------------------------------------------

    def summary(self):
        print("===== Atomic System =====")
        print(f"Atoms      : {self.natoms}")
        print(f"Si atoms   : {self.number_of(1)}")
        print(f"O atoms    : {self.number_of(2)}")
        print(f"Ca atoms   : {self.number_of(3)}")
        print()

        print("Simulation box")
        print(self.box)
        print()

        print("Box lengths (Å)")
        print(self.box_lengths)
