"""
Modifier-cation topology analysis for CaO-SiO2 glasses.

Analyzes the oxygen environments coordinated to Ca atoms,
using the validated Ca-O neighbor list and the QnAnalysis
BO/NBO classification.
"""

import numpy as np

from workflow.parameters import ATOM_TYPES
from workflow.qn import QnAnalysis


SI = 1
O = 2
CA = 3


class ModifierAnalysis:

    def __init__(self, system, neighbors):
        self.system = system
        self.neighbors = neighbors
        self.qn = QnAnalysis(system, neighbors)

    # --------------------------------------------------------
    # Ca atoms
    # --------------------------------------------------------

    def ca_atoms(self):
        """Return indices of all Ca atoms."""

        return [
            i
            for i, atom_type in enumerate(self.system.types)
            if atom_type == CA
        ]

    # --------------------------------------------------------
    # O neighbors of Ca
    # --------------------------------------------------------

    def oxygen_neighbors_of_ca(self, ca_index):
        """Return oxygen neighbors of one Ca atom."""

        return [
            j
            for j in self.neighbors.neighbors_of(ca_index)
            if self.system.types[j] == O
        ]

    # --------------------------------------------------------
    # Ca-O environment
    # --------------------------------------------------------

    def environment_of_ca(self, ca_index):
        """
        Classify the oxygen neighbors of one Ca atom.

        Returns
        -------
        dict
            Counts of BO, NBO, and other oxygen neighbors.
        """

        counts = {
            "BO": 0,
            "NBO": 0,
            "other": 0,
        }

        for oxygen in self.oxygen_neighbors_of_ca(ca_index):

            classification = self.qn.oxygen_class(oxygen)

            counts[classification] += 1

        return counts

    # --------------------------------------------------------
    # Global Ca environment
    # --------------------------------------------------------

    def global_environment(self):
        """
        Count Ca-O contacts by oxygen type.

        Returns
        -------
        dict
            Total BO, NBO, and other Ca-O contacts.
        """

        totals = {
            "BO": 0,
            "NBO": 0,
            "other": 0,
        }

        for ca in self.ca_atoms():

            environment = self.environment_of_ca(ca)

            for key in totals:
                totals[key] += environment[key]

        return totals

    # --------------------------------------------------------
    # Per-Ca statistics
    # --------------------------------------------------------

    def coordination_statistics(self):

        ca_indices = self.ca_atoms()

        total = []
        bo = []
        nbo = []
        other = []

        for ca in ca_indices:

            environment = self.environment_of_ca(ca)

            bo.append(environment["BO"])
            nbo.append(environment["NBO"])
            other.append(environment["other"])

            total.append(
                environment["BO"]
                + environment["NBO"]
                + environment["other"]
            )

        return {
            "number_ca": len(ca_indices),
            "mean_ca_o": float(np.mean(total)),
            "std_ca_o": float(np.std(total)),
            "mean_bo": float(np.mean(bo)),
            "mean_nbo": float(np.mean(nbo)),
            "mean_other": float(np.mean(other)),
        }

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    def report(self):

        totals = self.global_environment()
        stats = self.coordination_statistics()

        total_contacts = sum(totals.values())

        print()
        print("===== Ca Modifier Environment =====")

        print(
            f"Ca atoms          : "
            f"{stats['number_ca']}"
        )

        print(
            f"Mean Ca-O         : "
            f"{stats['mean_ca_o']:.3f}"
        )

        print(
            f"Std deviation     : "
            f"{stats['std_ca_o']:.3f}"
        )

        print()
        print("Ca-O contact classification")

        for key in ("BO", "NBO", "other"):

            count = totals[key]

            fraction = (
                100.0 * count / total_contacts
                if total_contacts
                else 0.0
            )

            print(
                f"{key:6s} : "
                f"{count:5d} "
                f"({fraction:6.2f} %)"
            )

        print()
        print("Mean oxygen contacts per Ca")

        print(
            f"BO     : "
            f"{stats['mean_bo']:.3f}"
        )

        print(
            f"NBO    : "
            f"{stats['mean_nbo']:.3f}"
        )

        print(
            f"Other  : "
            f"{stats['mean_other']:.3f}"
        )
