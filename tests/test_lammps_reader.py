from workflow.lammps_reader import read_lammps

system = read_lammps("simulations/data/glass_final_300K.data")

print()

print(system.natoms)
print(system.box)

print(system.positions.shape)
print(system.types.shape)
