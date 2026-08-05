import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# ==========================================================
# Parameters
# ==========================================================

DATA_FILE = "glass_final_300K.data"

SI_TYPE = 1
O_TYPE = 2
CA_TYPE = 3

SIO_CUTOFF = 2.325

# ==========================================================
# Read LAMMPS data file
# ==========================================================

box = np.zeros((3,2))

atoms = []

with open(DATA_FILE) as f:

    lines = f.readlines()

inside_atoms = False

for line in lines:

    if "xlo xhi" in line:
        s = line.split()
        box[0] = [float(s[0]), float(s[1])]

    elif "ylo yhi" in line:
        s = line.split()
        box[1] = [float(s[0]), float(s[1])]

    elif "zlo zhi" in line:
        s = line.split()
        box[2] = [float(s[0]), float(s[1])]

    elif line.startswith("Atoms"):
        inside_atoms = True
        continue

    elif inside_atoms:

        if len(line.strip()) == 0:
            continue
        s = line.split()

        if len(s) < 9:
            continue

        atom_id = int(s[0])
        atom_type = int(s[1])
        charge = float(s[2])

        x = float(s[3])
        y = float(s[4])
        z = float(s[5])

        atoms.append([atom_id, atom_type, x, y, z])

atoms = np.array(atoms)

print(f"Read {len(atoms)} atoms")
# ==========================================================
# Build atom lists
# ==========================================================

Lx = box[0,1] - box[0,0]
Ly = box[1,1] - box[1,0]
Lz = box[2,1] - box[2,0]

positions = {}
types = {}

Si_atoms = []
O_atoms = []

for atom in atoms:

    aid = int(atom[0])
    atype = int(atom[1])

    positions[aid] = np.array(atom[2:5])
    types[aid] = atype

    if atype == SI_TYPE:
        Si_atoms.append(aid)

    elif atype == O_TYPE:
        O_atoms.append(aid)

print(f"Si atoms : {len(Si_atoms)}")
print(f"O atoms  : {len(O_atoms)}")

# ==========================================================
# Minimum image convention
# ==========================================================

def minimum_image(dr, box_length):

    if dr > 0.5 * box_length:
        dr -= box_length
    elif dr < -0.5 * box_length:
        dr += box_length

    return dr

def distance_vector(r1, r2):

    dx = minimum_image(r2[0] - r1[0], Lx)
    dy = minimum_image(r2[1] - r1[1], Ly)
    dz = minimum_image(r2[2] - r1[2], Lz)

    return np.array([dx, dy, dz])

# ==========================================================
# Build Si-O neighbor list
# ==========================================================

Si_neighbors = defaultdict(list)
O_neighbors = defaultdict(list)

for si in Si_atoms:

    r_si = positions[si]

    for o in O_atoms:

        r_o = positions[o]

        dv = distance_vector(r_si, r_o)
        d = np.linalg.norm(dv)

        if d <= SIO_CUTOFF:

            Si_neighbors[si].append(o)
            O_neighbors[o].append(si)

print("Average Si neighbors:",
      np.mean([len(Si_neighbors[s]) for s in Si_atoms]))

# ==========================================================
# BO / NBO classification
# ==========================================================

BO = 0
NBO = 0
FREE = 0
TRICLUSTER = 0

Qn = defaultdict(int)

for o in O_atoms:

    nSi = len(O_neighbors[o])

    if nSi == 2:
        BO += 1

    elif nSi == 1:
        NBO += 1

    elif nSi == 0:
        FREE += 1

    else:
        TRICLUSTER += 1


# ==========================================================
# Qn distribution
# ==========================================================

for si in Si_atoms:

    q = 0

    for o in Si_neighbors[si]:

        if len(O_neighbors[o]) == 2:
            q += 1

    Qn[q] += 1

# ==========================================================
# Consistency check
# ==========================================================

bad_si = 0

for si in Si_atoms:

    total = len(Si_neighbors[si])

    bo = 0
    nbo = 0

    for o in Si_neighbors[si]:

        if len(O_neighbors[o]) == 2:
            bo += 1
        elif len(O_neighbors[o]) == 1:
            nbo += 1

    if bo + nbo != total:
        bad_si += 1

print("\n========== Consistency Check ==========")
print(f"Inconsistent Si atoms : {bad_si}")
print(f"Total Si atoms        : {len(Si_atoms)}")

# ==========================================================
# Print results
# ==========================================================

print("\n========== BO / NBO ==========")
print(f"BO          : {BO}")
print(f"NBO         : {NBO}")
print(f"Free O      : {FREE}")
print(f"Tri-cluster : {TRICLUSTER}")

print("\n========== Qn ==========")

for n in range(5):
    print(f"Q{n} : {Qn[n]}")

