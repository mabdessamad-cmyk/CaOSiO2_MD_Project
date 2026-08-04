from workflow.lammps_reader import read_lammps

system = read_lammps("simulations/data/glass_final_300K.data")

print("\n===== Summary =====")
print(f"Atoms      : {system.natoms}")
print(f"Box shape  : {system.box.shape}")
print(f"Positions  : {system.positions.shape}")
print(f"Types      : {system.types.shape}")
print(f"Charges    : {system.charges.shape}")

print("\nFirst atom")
print(f"ID         : {system.ids[0]}")
print(f"Type       : {system.types[0]}")
print(f"Charge     : {system.charges[0]}")
print(f"Position   : {system.positions[0]}")
