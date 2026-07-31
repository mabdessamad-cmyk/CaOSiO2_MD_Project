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

def classify_oxygens(owners):
    """
    Count bridging and non-bridging oxygens.
    """

    bo = 0
    nbo = 0

    for o in owners:

        if len(o) == 2:
            bo += 1

        elif len(o) == 1:
            nbo += 1

    return bo, nbo
from geometry import distance

def find_si_neighbors(si_positions, box_length, cutoff):
    """
    Find neighboring Si atoms within a cutoff distance.

    Returns
    -------
    neighbors : list of lists
        neighbors[i] contains the indices of Si atoms
        neighboring Si atom i.
    """

    n = len(si_positions)

    neighbors = [[] for _ in range(n)]

    for i in range(n - 1):

        for j in range(i + 1, n):

            d = distance(
                si_positions[i],
                si_positions[j],
                box_length
            )

            if d <= cutoff:

                neighbors[i].append((j, d))
                neighbors[j].append((i, d))
    for i in range(n):

        neighbors[i].sort(key=lambda x: x[1])

        neighbors[i] = [
            index
            for index, dist in neighbors[i][:4]
        ]

    mutual = [[] for _ in range(n)]

    for i in range(n):

        for j in neighbors[i]:

            if i in neighbors[j]:
                mutual[i].append(j)

    return mutual


def check_neighbor_symmetry(neighbors):
    """
    Count asymmetric neighbor relationships.
    """

    bad = 0

    for i, lst in enumerate(neighbors):

        for j in lst:

            if i not in neighbors[j]:
                bad += 1

    return bad

def connect_tetrahedra(
    si_positions,
    neighbors,
    box_length,
    target_bo=None
):
    """
    Create one bridging oxygen (BO) for each mutual Si-Si pair.

    Returns
    -------
    oxygen_positions : list of ndarray
    oxygen_owners : list of list
    """

    oxygen_positions = []
    oxygen_owners = []

    for i, lst in enumerate(neighbors):

        for j in lst:

            if j <= i:
                continue

            vec = si_positions[j] - si_positions[i]

            # Minimum-image convention
            vec -= box_length * np.round(vec / box_length)

            d = np.linalg.norm(vec)

            direction = vec / d

            pos = si_positions[i] + 1.62 * direction
            pos %= box_length

            oxygen_positions.append(pos)
            oxygen_owners.append([i, j])

    if target_bo is not None and len(oxygen_positions) > target_bo:

        oxygen_positions = oxygen_positions[:target_bo]
        oxygen_owners = oxygen_owners[:target_bo]

    return oxygen_positions, oxygen_owners

def validate_coordination(oxygen_owners, n_si):
    """
    Compute the oxygen coordination number of each Si atom.
    """

    coordination = [0] * n_si

    for owners in oxygen_owners:

        for si in owners:

            coordination[si] += 1

    return coordination

def complete_tetrahedra(
    si_positions,
    oxygen_positions,
    oxygen_owners,
    coordination
):
    tetra_dirs = [
        np.array([ 1,  1,  1], dtype=float),
        np.array([ 1, -1, -1], dtype=float),
        np.array([-1,  1, -1], dtype=float),
        np.array([-1, -1,  1], dtype=float),
    ]

    tetra_dirs = [v / np.linalg.norm(v) for v in tetra_dirs]

    for i in range(len(si_positions)):

        missing = 4 - coordination[i]

        for k in range(missing):

            direction = tetra_dirs[coordination[i]]
            offset = 1.62 * direction

            new_position = si_positions[i] + offset

            oxygen_positions.append(new_position)

            oxygen_owners.append([i])

    return oxygen_positions, oxygen_owners
