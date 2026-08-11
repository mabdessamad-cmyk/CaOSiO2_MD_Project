from workflow.lammps_reader import read_lammps
from workflow.density import DensityAnalysis


system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

analysis = DensityAnalysis(system)

composition = analysis.composition()
volume = analysis.volume_angstrom3()
mass = analysis.total_mass_amu()
density = analysis.density_g_cm3()

ca_fraction = composition["Ca"] / (
    composition["Ca"] + composition["Si"]
)

assert system.natoms == 6480

assert composition["Si"] == 1680
assert composition["O"] == 4080
assert composition["Ca"] == 720

assert abs(ca_fraction - 0.30) < 1.0e-12

assert volume > 0.0
assert mass > 0.0
assert density > 0.0

assert abs(volume - 85722.161) < 1.0
assert abs(density - 2.7375) < 1.0e-4

print("===== Density Validation =====")
print(f"CaO fraction : {ca_fraction:.3f}")
print(f"Volume       : {volume:.3f} Å^3")
print(f"Mass         : {mass:.4f} amu")
print(f"Density      : {density:.4f} g/cm^3")
print("Density validation: PASS")
