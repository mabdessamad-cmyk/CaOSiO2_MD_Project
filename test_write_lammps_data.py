from generator import GlassGenerator
from writer import write_lammps_data

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
    seed=12345,
)

cfg = g.generate()

write_lammps_data(
    cfg,
    "output/glass_54.data",
)

print("LAMMPS data file successfully written.")
