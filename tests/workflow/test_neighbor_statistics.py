from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS


system = read_lammps("simulations/data/glass_final_300K.data")

neighbors = NeighborList(system, BOND_CUTOFFS)
neighbors.build()

print()
print("===== Neighbor Statistics =====")

for atom_type in [1, 2, 3]:

    indices = [i for i in range(system.natoms)
               if system.types[i] == atom_type]

    coord = [neighbors.coordination(i) for i in indices]

    print()
    print(f"Atom type {atom_type}")
    print(f"Number of atoms : {len(indices)}")
    print(f"Average coordination : {sum(coord)/len(coord):.3f}")
    print(f"Minimum : {min(coord)}")
    print(f"Maximum : {max(coord)}")
