"""
workflow/qn.py

Qn, bridging-oxygen (BO), and non-bridging-oxygen (NBO)
analysis for CaO-SiO2 glasses.

Definitions
-----------
BO:
    Oxygen connected to two Si atoms.

NBO:
    Oxygen connected to one Si atom.

Qn:
    A SiO4 tetrahedron containing n bridging oxygens.
"""

import numpy as np

from workflow.parameters import ATOM_TYPES


SI = 1
O = 2
CA = 3


class QnAnalysis:
    """
    Analyze Qn species and BO/NBO populations.

    Parameters
    ----------
    system : AtomicSystem
        Atomic configuration.

    neighbors : NeighborList
        Previously constructed and built neighbor list.
    """

    def __init__(self, system, neighbors):

        self.system = system
        self.neighbors = neighbors

    # --------------------------------------------------------
    # Si neighbors of an oxygen
    # --------------------------------------------------------

    def si_neighbors_of_oxygen(self, oxygen_index):
        """
        Return Si neighbors of one oxygen atom.

        Parameters
        ----------
        oxygen_index : int
            Array index of the oxygen atom.

        Returns
        -------
        list
            Indices of neighboring Si atoms.
        """

        return [
            j
            for j in self.neighbors.neighbors_of(oxygen_index)
            if self.system.types[j] == SI
        ]

    # --------------------------------------------------------
    # Oxygen classification
    # --------------------------------------------------------

    def oxygen_class(self, oxygen_index):
        """
        Classify one oxygen as BO or NBO.

        Returns
        -------
        str
            "BO", "NBO", or "other".
        """

        n_si = len(self.si_neighbors_of_oxygen(oxygen_index))

        if n_si == 2:
            return "BO"

        if n_si == 1:
            return "NBO"

        return "other"

    # --------------------------------------------------------
    # Oxygen statistics
    # --------------------------------------------------------

    def oxygen_statistics(self):
        """
        Calculate BO/NBO statistics.

        Returns
        -------
        dict
            Counts and fractions of BO, NBO, and other oxygen atoms.
        """

        oxygen_indices = np.where(
            self.system.types == O
        )[0]

        bo = 0
        nbo = 0
        other = 0

        for i in oxygen_indices:

            classification = self.oxygen_class(i)

            if classification == "BO":
                bo += 1

            elif classification == "NBO":
                nbo += 1

            else:
                other += 1

        total = len(oxygen_indices)

        return {
            "total": total,
            "BO": bo,
            "NBO": nbo,
            "other": other,
            "BO_fraction": bo / total if total else np.nan,
            "NBO_fraction": nbo / total if total else np.nan,
            "other_fraction": other / total if total else np.nan,
        }

    # --------------------------------------------------------
    # Qn of one Si atom
    # --------------------------------------------------------

    def qn_of_si(self, si_index):
        """
        Determine Qn value of one Si atom.

        Parameters
        ----------
        si_index : int
            Array index of the Si atom.

        Returns
        -------
        int
            Number of bridging oxygens around the Si atom.
        """

        n_bo = 0

        for j in self.neighbors.neighbors_of(si_index):

            if self.system.types[j] != O:
                continue

            if self.oxygen_class(j) == "BO":
                n_bo += 1

        return n_bo

    # --------------------------------------------------------
    # Qn distribution
    # --------------------------------------------------------

    def qn_distribution(self):
        """
        Calculate the Qn distribution of all Si atoms.

        Returns
        -------
        dict
            Counts for Q0 through Q4.
        """

        si_indices = np.where(
            self.system.types == SI
        )[0]

        distribution = {
            0: 0,
            1: 0,
            2: 0,
            3: 0,
            4: 0,
        }

        invalid = 0

        for i in si_indices:

            qn = self.qn_of_si(i)

            if qn in distribution:
                distribution[qn] += 1
            else:
                invalid += 1

        distribution["invalid"] = invalid

        return distribution

    # --------------------------------------------------------
    # Qn fractions
    # --------------------------------------------------------

    def qn_fractions(self):
        """
        Calculate normalized Qn fractions.

        Returns
        -------
        dict
            Fractions of Q0 through Q4.
        """

        distribution = self.qn_distribution()

        total = sum(
            distribution[n]
            for n in range(5)
        )

        if total == 0:
            return {
                n: np.nan
                for n in range(5)
            }

        return {
            n: distribution[n] / total
            for n in range(5)
        }

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    def report(self):
        """
        Print a human-readable Qn / BO / NBO report.
        """

        print()
        print("===== Qn / BO / NBO Analysis =====")

        # Oxygen statistics

        oxygen = self.oxygen_statistics()

        print()
        print("Oxygen classification")

        print(
            f"Total O atoms : {oxygen['total']}"
        )

        print(
            f"BO           : {oxygen['BO']}"
            f" ({100 * oxygen['BO_fraction']:.2f} %)"
        )

        print(
            f"NBO          : {oxygen['NBO']}"
            f" ({100 * oxygen['NBO_fraction']:.2f} %)"
        )

        print(
            f"Other        : {oxygen['other']}"
            f" ({100 * oxygen['other_fraction']:.2f} %)"
        )

        # Qn distribution

        distribution = self.qn_distribution()
        fractions = self.qn_fractions()

        print()
        print("Qn distribution")

        for n in range(5):

            print(
                f"Q{n} : "
                f"{distribution[n]:4d} atoms "
                f"({100 * fractions[n]:6.2f} %)"
            )

        print(
            f"Invalid : {distribution['invalid']}"
        )
