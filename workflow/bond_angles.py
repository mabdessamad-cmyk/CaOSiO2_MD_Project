"""
Bond-angle analysis for CaO-SiO2 glasses.
"""

import numpy as np


class BondAngleAnalysis:
    """
    Calculate bond-angle distributions for selected triplets.

    Angles are calculated using minimum-image vectors.
    """

    def __init__(self, system, neighbors):
        self.system = system
        self.neighbors = neighbors

    def _vector(self, i, j):
        """
        Minimum-image vector from atom i to atom j.
        """

        dr = (
            self.system.positions[j]
            - self.system.positions[i]
        )

        box_lengths = self.system.box_lengths

        dr -= box_lengths * np.rint(dr / box_lengths)

        return dr

    def angle(self, i, j, k):
        """
        Calculate angle i-j-k in degrees.
        """

        v1 = self._vector(j, i)
        v2 = self._vector(j, k)

        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)

        if n1 == 0.0 or n2 == 0.0:
            raise ValueError("Zero-length bond encountered.")

        cosine = np.dot(v1, v2) / (n1 * n2)

        cosine = np.clip(cosine, -1.0, 1.0)

        return float(np.degrees(np.arccos(cosine)))

    def angles_for_triplet(self, outer_type, center_type):
        """
        Calculate all angles outer-center-outer.

        Parameters
        ----------
        outer_type : int
            Atom type at both ends.

        center_type : int
            Atom type at the center.

        Returns
        -------
        numpy.ndarray
            Angles in degrees.
        """

        angles = []

        for j in np.where(
            self.system.types == center_type
        )[0]:

            neighbors = [
                i
                for i in self.neighbors.neighbors_of(j)
                if self.system.types[i] == outer_type
            ]

            for a in range(len(neighbors)):
                for b in range(a + 1, len(neighbors)):

                    i = neighbors[a]
                    k = neighbors[b]

                    angles.append(
                        self.angle(i, j, k)
                    )

        return np.asarray(angles, dtype=float)

    def statistics(self, outer_type, center_type):
        """
        Calculate basic statistics for a bond-angle distribution.
        """

        values = self.angles_for_triplet(
            outer_type,
            center_type
        )

        if len(values) == 0:
            return {
                "number": 0,
                "mean": np.nan,
                "std": np.nan,
                "min": np.nan,
                "max": np.nan,
            }

        return {
            "number": len(values),
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
        }

    def report(self):
        """
        Print statistics for selected bond-angle distributions.
        """

        analyses = [
            ("O-Si-O", 2, 1),
            ("Si-O-Si", 1, 2),
            ("O-Ca-O", 2, 3),
        ]

        print()
        print("===== Bond Angle Analysis =====")

        for name, outer, center in analyses:

            stats = self.statistics(
                outer,
                center
            )

            print()
            print(name)
            print(f"Number of angles : {stats['number']}")
            print(f"Mean angle       : {stats['mean']:.3f} deg")
            print(f"Std deviation    : {stats['std']:.3f} deg")
            print(f"Minimum          : {stats['min']:.3f} deg")
            print(f"Maximum          : {stats['max']:.3f} deg")
