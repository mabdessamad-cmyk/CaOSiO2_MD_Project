import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.ndimage import gaussian_filter1d
# ==========================================================
# Parameters
# ==========================================================

DATA_FILE = "glass_final_300K.data"

SI_TYPE = 1
O_TYPE = 2
CA_TYPE = 3

SIO_CUTOFF = 2.325      # Å

NBINS = 180

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
# O-Si-O angles
# ==========================================================

osi_angles = []

for si in Si_atoms:

    neigh = Si_neighbors[si]

    if len(neigh) < 2:
        continue

    r_si = positions[si]

    for i in range(len(neigh)-1):
        for j in range(i+1, len(neigh)):

            r1 = positions[neigh[i]]
            r2 = positions[neigh[j]]

            v1 = distance_vector(r_si, r1)
            v2 = distance_vector(r_si, r2)

            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)

            c = np.dot(v1, v2)
            c = np.clip(c, -1.0, 1.0)

            angle = np.degrees(np.arccos(c))

            osi_angles.append(angle)

osi_angles = np.array(osi_angles)

print(f"O-Si-O angles computed : {len(osi_angles)}")
print(f"Average O-Si-O angle   : {osi_angles.mean():.2f}°")

# ==========================================================
# Si-O-Si angles
# ==========================================================

siosi_angles = []

for o in O_atoms:

    neigh = O_neighbors[o]

    if len(neigh) < 2:
        continue

    r_o = positions[o]

    for i in range(len(neigh)-1):
        for j in range(i+1, len(neigh)):

            r1 = positions[neigh[i]]
            r2 = positions[neigh[j]]

            v1 = distance_vector(r_o, r1)
            v2 = distance_vector(r_o, r2)

            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)

            c = np.dot(v1, v2)
            c = np.clip(c, -1.0, 1.0)

            angle = np.degrees(np.arccos(c))

            siosi_angles.append(angle)

siosi_angles = np.array(siosi_angles)

print(f"Si-O-Si angles computed : {len(siosi_angles)}")
print(f"Average Si-O-Si angle   : {siosi_angles.mean():.2f}°")

# ==========================================================
# Save angle data
# ==========================================================

np.savetxt("osi_o_angles.dat", osi_angles, fmt="%.6f")
np.savetxt("si_o_si_angles.dat", siosi_angles, fmt="%.6f")

# ==========================================================
# O-Si-O distribution
# ==========================================================

hist, bins = np.histogram(
    osi_angles,
    bins=360,
    range=(0,180),
    density=True
)

centers = 0.5 * (bins[:-1] + bins[1:])

hist = gaussian_filter1d(hist, sigma=2)

plt.figure(figsize=(7,5))

plt.plot(centers, hist, linewidth=2.5)

plt.xlabel("O–Si–O angle (°)", fontsize=16)
plt.ylabel("Probability density", fontsize=16)

plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

plt.xlim(85,135)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig("osi_o_distribution.png", dpi=600)

# ==========================================================
# Si-O-Si distribution
# ==========================================================

hist, bins = np.histogram(
    siosi_angles,
    bins=360,
    range=(0,180),
    density=True
)

centers = 0.5 * (bins[:-1] + bins[1:])

hist = gaussian_filter1d(hist, sigma=2)

plt.figure(figsize=(7,5))

plt.plot(centers, hist, linewidth=2.5)

plt.xlabel("Si–O–Si angle (°)", fontsize=16)
plt.ylabel("Probability density", fontsize=16)

plt.xticks(fontsize=13)
plt.yticks(fontsize=13)

plt.xlim(110,180)

plt.grid(alpha=0.3)

plt.tight_layout()

plt.savefig("si_o_si_distribution.png", dpi=600)
