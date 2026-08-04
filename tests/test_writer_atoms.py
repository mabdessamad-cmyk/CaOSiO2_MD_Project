from generator import GlassGenerator
from writer import write_atoms

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
)

cfg = g.generate()

with open("output/atoms_test.data", "w") as f:
    write_atoms(f, cfg)

print("Atoms section written successfully.")
