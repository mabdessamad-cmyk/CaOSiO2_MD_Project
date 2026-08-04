import numpy as np

# Read only RDF rows (14 columns)
rows = []

with open("rdf_numeric.dat") as f:
    for line in f:
        s = line.split()
        if len(s) == 14:
            rows.append([float(x) for x in s])

data = np.array(rows)

r = data[:,1]

g_sio = data[:,4]   # Si-O
g_cao = data[:,10]  # Ca-O

# ----- Si-O -----
imax = np.argmax(g_sio)
imin = imax + np.argmin(g_sio[imax:])

print("Si-O")
print(f"Peak    : r = {r[imax]:.3f} Å   g(r) = {g_sio[imax]:.3f}")
print(f"Minimum : r = {r[imin]:.3f} Å   g(r) = {g_sio[imin]:.3f}")

# ----- Ca-O -----
imax = np.argmax(g_cao)
imin = imax + np.argmin(g_cao[imax:])

print("\nCa-O")
print(f"Peak    : r = {r[imax]:.3f} Å   g(r) = {g_cao[imax]:.3f}")
print(f"Minimum : r = {r[imin]:.3f} Å   g(r) = {g_cao[imin]:.3f}")
