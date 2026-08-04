import numpy as np

reading = False
atoms = []

with open("coordination.dump") as f:
    for line in f:
        if line.startswith("ITEM: ATOMS"):
            reading = True
            atoms = []
            continue

        if line.startswith("ITEM:"):
            reading = False
            continue

        if reading:
            atoms.append([float(x) for x in line.split()])

# Keep only the LAST snapshot
data = np.array(atoms)

atom_type = data[:,1].astype(int)
si_o = data[:,2]
ca_o = data[:,3]

print(f"Average Si-O CN = {si_o[atom_type==1].mean():.3f}")
print(f"Average Ca-O CN = {ca_o[atom_type==3].mean():.3f}")

print("\nSi coordination distribution")
vals, cnt = np.unique(si_o[atom_type==1], return_counts=True)
for v, c in zip(vals, cnt):
    print(f"CN={int(v)} : {c}")

print("\nCa coordination distribution")
vals, cnt = np.unique(ca_o[atom_type==3], return_counts=True)
for v, c in zip(vals, cnt):
    print(f"CN={int(v)} : {c}")
