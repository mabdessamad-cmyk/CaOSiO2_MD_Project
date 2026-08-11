"""
Density analysis for CaO-SiO2 glasses.
"""

import numpy as np


ATOMIC_MASSES = {
    1: 28.0855,   # Si
    2: 15.9994,   # O
    3: 40.078,    # Ca
}

AMU_TO_GRAM = 1.66053906660e-24
ANGSTROM3_TO_CM3 = 1.0e-24


class DensityAnalysis:
    """
    Calculate mass, volume, and density of an AtomicSystem.
    """

    def __init__(self, system):
        self.system = system

    def total_mass_amu(self):
        """
        Total system mass in atomic mass units.
        """
        return float(
            sum(
                ATOMIC_MASSES[int(atom_type)]
                for atom_type in self.system.types
            )
        )

    def volume_angstrom3(self):
        """
        Simulation volume in Å^3.
        """
        lengths = self.system.box_lengths

        if np.any(lengths <= 0):
            raise ValueError("Invalid simulation box dimensions.")

        return float(np.prod(lengths))

    def density_g_cm3(self):
        """
        Mass density in g/cm^3.
        """
        mass_g = self.total_mass_amu() * AMU_TO_GRAM
        volume_cm3 = self.volume_angstrom3() * ANGSTROM3_TO_CM3

        return mass_g / volume_cm3

    def composition(self):
        """
        Return atom counts by chemical species.
        """
        return {
            "Si": int(np.sum(self.system.types == 1)),
            "O": int(np.sum(self.system.types == 2)),
            "Ca": int(np.sum(self.system.types == 3)),
        }

    def report(self):
        """
        Print density information.
        """
        lengths = self.system.box_lengths
        composition = self.composition()

        print("===== Density Analysis =====")
        print(f"Box lengths : {lengths}")
        print(f"Volume      : {self.volume_angstrom3():.3f} Å^3")
        print(f"Composition : {composition}")
        print(f"Mass        : {self.total_mass_amu():.4f} amu")
        print(f"Density     : {self.density_g_cm3():.4f} g/cm^3")
