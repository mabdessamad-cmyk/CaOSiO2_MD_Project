"""
Topology construction for CaO-SiO2 glasses.
"""

import numpy as np

# Ideal tetrahedral directions
TETRA = np.array([
    [ 1.0,  1.0,  1.0],
    [ 1.0, -1.0, -1.0],
    [-1.0,  1.0, -1.0],
    [-1.0, -1.0,  1.0],
], dtype=float)

# Normalize the vectors
for i in range(4):
    TETRA[i] /= np.linalg.norm(TETRA[i])


def tetrahedron(si_position, bond_length):
    """
    Return the four oxygen positions around one silicon atom.
    """

    oxygens = []

    for direction in TETRA:
        oxygens.append(si_position + bond_length * direction)

    return oxygens

def generate_tetrahedra(si_positions, bond_length=1.62):

    oxygen_positions = []
    oxygen_owner = []

    for si_id, si in enumerate(si_positions):

        oxygens = tetrahedron(si, bond_length)

        for o in oxygens:
            oxygen_positions.append(o)
            oxygen_owner.append(si_id)

    return oxygen_positions, oxygen_owner


import numpy as np
def merge_oxygens(oxygen_positions, oxygen_owner, tolerance=0.25):

    """
    Merge oxygen atoms that are closer than 'tolerance'.

    Parameters
    ----------
    oxygen_positions : list of ndarray
    tolerance : float
        Merge distance in Å.

    Returns
    -------
    unique_oxygens : list
    mapping : list
        mapping[i] = index of the merged oxygen corresponding
        to temporary oxygen i.
    """

    unique_oxygens = []
    mapping = []
    owners = []

    tol2 = tolerance * tolerance

    for k, pos in enumerate(oxygen_positions):

        found = False

        for i, u in enumerate(unique_oxygens):

            if np.sum((pos - u)**2) < tol2:

                # average position
                unique_oxygens[i] = 0.5 * (u + pos)

                mapping.append(i)

            if oxygen_owner[k] not in owners[i]:
                owners[i].append(oxygen_owner[k])

                found = True
                break

        if not found:

            unique_oxygens.append(pos.copy())

            mapping.append(len(unique_oxygens) - 1)

            owners.append([oxygen_owner[k]])

    return unique_oxygens, mapping, owners
