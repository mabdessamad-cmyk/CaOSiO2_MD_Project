from generator import GlassGenerator

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
)

cfg = g.generate()

print()
print("Generated atoms:", cfg.natoms)
