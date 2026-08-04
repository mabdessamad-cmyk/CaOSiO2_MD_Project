from generator import GlassGenerator

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
)

g.place_species("Si")

print()

print("Atoms inserted:", g.cfg.natoms)
print("Si stored:", g.cfg.species.count("Si"))
