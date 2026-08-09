"""
Find the first minimum after the main RDF peak
for Si-O and Ca-O correlations.
"""

from pathlib import Path

import numpy as np


# ------------------------------------------------------------
# Locate RDF data relative to the project root
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[3]

RDF_FILE = PROJECT_ROOT / "analysis" / "tables" / "rdf_numeric.dat"


# ------------------------------------------------------------
# Read RDF data
# ------------------------------------------------------------

rows = []

with RDF_FILE.open() as f:

    for line in f:

        s = line.split()

        if len(s) == 14:
            rows.append([float(x) for x in s])


data = np.array(rows)


if data.size == 0:
    raise RuntimeError(
        f"No valid RDF rows found in {RDF_FILE}"
    )


# ------------------------------------------------------------
# Extract columns
# ------------------------------------------------------------

r = data[:, 1]

g_sio = data[:, 4]     # Si-O
g_cao = data[:, 10]    # Ca-O


# ------------------------------------------------------------
# Si-O
# ------------------------------------------------------------

imax = np.argmax(g_sio)

imin = imax + np.argmin(
    g_sio[imax:]
)

print("===== RDF Minima =====")

print()
print("Si-O")
print(
    f"Peak    : r = {r[imax]:.3f} Å   "
    f"g(r) = {g_sio[imax]:.3f}"
)
print(
    f"Minimum : r = {r[imin]:.3f} Å   "
    f"g(r) = {g_sio[imin]:.3f}"
)


# ------------------------------------------------------------
# Ca-O
# ------------------------------------------------------------

imax = np.argmax(g_cao)

imin = imax + np.argmin(
    g_cao[imax:]
)

print()
print("Ca-O")
print(
    f"Peak    : r = {r[imax]:.3f} Å   "
    f"g(r) = {g_cao[imax]:.3f}"
)
print(
    f"Minimum : r = {r[imin]:.3f} Å   "
    f"g(r) = {g_cao[imin]:.3f}"
)
