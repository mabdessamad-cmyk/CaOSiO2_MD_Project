from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS

system = read_lammps("simulations/data/glass_final_300K.data")


neighbors = NeighborList(system, BOND_CUTOFFS)

neighbors.build()

print()

print("Atom 0")
print("ID          :", system.ids[0])
print("Type        :", system.types[0])
print("Charge      :", system.charges[0])
print("Position    :", system.positions[0])
print("Neighbors :", neighbors.neighbors_of(0))
print("Coordination :", neighbors.coordination(0))
