"""
Neighbor-list construction for atomic systems.
"""

import numpy as np

from workflow.geometry import distance


class NeighborList:
    """
    Neighbor list for an AtomicSystem.
    """

    def __init__(self, system, cutoffs):
        """
        Parameters
        ----------
        system : AtomicSystem

        cutoffs : dict
            Dictionary containing cutoff distances
            for each atom-type pair.
        """


        self.system = system
        self.cutoffs = cutoffs
        self.list = []

    def build(self):
        """
        Construct the neighbor list using a brute-force search.
        """

        natoms = self.system.natoms

        self.list = [[] for _ in range(natoms)]

        print(f"Building neighbor list ({natoms} atoms)...")

        for i in range(natoms):

            for j in range(i + 1, natoms):

                d = distance(
                    self.system.positions[i],
                    self.system.positions[j],
                    self.system.box_lengths,
                )

                type_i = self.system.types[i]
                type_j = self.system.types[j]

                cutoff = self.cutoffs.get((type_i, type_j), None)

                if cutoff is None:
                    continue

                if d <= cutoff:

                    self.list[i].append(j)
                    self.list[j].append(i)


        print("Neighbor list completed.")

    def neighbors_of(self, atom):
        return self.list[atom]

    def coordination(self, atom):
        return len(self.list[atom])
