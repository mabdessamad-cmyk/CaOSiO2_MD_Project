import numpy as np
import matplotlib.pyplot as plt

rows = []

with open("partial_rdf.dat") as f:
    for line in f:
        if line.startswith("#"):
            continue

        s = line.split()

        # Keep only RDF rows (14 columns)
        if len(s) == 14:
            rows.append([float(x) for x in s])

data = np.array(rows)

r = data[:,1]

# g(r) columns
g_sisi = data[:,2]
g_sio  = data[:,4]
g_sica = data[:,6]
g_oo   = data[:,8]
g_oca  = data[:,10]
g_caca = data[:,12]

plt.figure(figsize=(8,6))

plt.plot(r, g_sisi, label="Si-Si", lw=2)
plt.plot(r, g_sio,  label="Si-O",  lw=2)
plt.plot(r, g_sica, label="Si-Ca", lw=2)
plt.plot(r, g_oo,   label="O-O",   lw=2)
plt.plot(r, g_oca,  label="O-Ca",  lw=2)
plt.plot(r, g_caca, label="Ca-Ca", lw=2)

plt.xlim(0,8)
plt.ylim(bottom=0)

plt.xlabel("r (Å)", fontsize=13)
plt.ylabel("g(r)", fontsize=13)
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("partial_rdf.png", dpi=300)

plt.show()
