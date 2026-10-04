from workflow.lammps_reader import read_lammps

system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

system.summary()

print()
print("First atom")
print("ID       :", system.ids[0])
print("Type     :", system.types[0])
print("Charge   :", system.charges[0])
print("Position :", system.positions[0])
