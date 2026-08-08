"""
workflow/coordination.py

Coordination-number analysis for CaO-SiO2 glasses.

This module uses the validated NeighborList infrastructure and reports
species-specific coordination numbers.
"""

import numpy as np

from workflow.parameters import ATOM_TYPES


class CoordinationAnalysis:
    """
    Calculate species-specific coordination numbers.

    Parameters
    ----------
    system : AtomicSystem
        Atomic configuration.
    neighbors : NeighborList
        Previously constructed neighbor list.
    """

    def __init__(self, system, neighbors):

        self.system = system
        self.neighbors = neighbors

    # ------------------------------------------------------------
    # Coordination of one atom
    # ------------------------------------------------------------

    def coordination(self, atom_index, neighbor_type):
        """
        Count neighbors of a specified atom type.

        Parameters
        ----------
        atom_index : int
            Index of the central atom.

        neighbor_type : int
            LAMMPS atom type of the neighbors.

        Returns
        -------
        int
            Number of neighbors of the requested type.
        """

        count = 0

        for j in self.neighbors.neighbors_of(atom_index):

            if self.system.types[j] == neighbor_type:
                count += 1

        return count

    # ------------------------------------------------------------
    # Coordination array
    # ------------------------------------------------------------

    def coordination_array(self, central_type, neighbor_type):
        """
        Calculate coordination numbers for all atoms of a given type.

        Parameters
        ----------
        central_type : int
            Atom type of the central atoms.

        neighbor_type : int
            Atom type being counted.

        Returns
        -------
        numpy.ndarray
            Coordination number for each central atom.
        """

        indices = np.where(self.system.types == central_type)[0]

        values = np.zeros(len(indices), dtype=int)

        for k, i in enumerate(indices):

            values[k] = self.coordination(
                i,
                neighbor_type
            )

        return values

    # ------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------

    def statistics(self, central_type, neighbor_type):
        """
        Calculate basic coordination statistics.

        Returns
        -------
        dict
            Number of atoms, mean, standard deviation,
            minimum, and maximum coordination.
        """

        values = self.coordination_array(
            central_type,
            neighbor_type
        )

        if len(values) == 0:

            return {
                "number": 0,
                "mean": np.nan,
                "std": np.nan,
                "min": None,
                "max": None,
            }

        return {
            "number": len(values),
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "min": int(np.min(values)),
            "max": int(np.max(values)),
        }

    # ------------------------------------------------------------
    # Human-readable report
    # ------------------------------------------------------------

    def report(self):
        """
        Print coordination statistics for the CaO-SiO2 system.
        """

        si = ATOM_TYPES
        print()
        print("===== Coordination Analysis =====")

        analyses = [
            ("Si-O", 1, 2),
            ("O-Si", 2, 1),
            ("Ca-O", 3, 2),
            ("O-Ca", 2, 3),
        ]

        for name, central, neighbor in analyses:

            stats = self.statistics(
                central,
                neighbor
            )

            print()
            print(name)

            print(
                f"Number of atoms : {stats['number']}"
            )

            print(
                f"Mean coordination : "
                f"{stats['mean']:.3f}"
            )

            print(
                f"Std deviation : "
                f"{stats['std']:.3f}"
            )

            print(
                f"Minimum : {stats['min']}"
            )

            print(
                f"Maximum : {stats['max']}"
            )
