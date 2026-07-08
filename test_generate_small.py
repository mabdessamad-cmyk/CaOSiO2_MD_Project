from generator import GlassGenerator

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
)

g.place_species("Si")
g.place_species("Ca")
g.place_species("O")

print()
print("Final statistics")
print("----------------")
print("Total atoms :", g.cfg.natoms)

print("Si =", g.cfg.species.count("Si"))
print("Ca =", g.cfg.species.count("Ca"))
print("O  =", g.cfg.species.count("O"))
