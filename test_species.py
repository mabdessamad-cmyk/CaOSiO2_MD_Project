from generator import GlassGenerator

g = GlassGenerator(
    x=0.30,
    formula_units=2400,
    density=1.50,
    seed=12345,
)

species = g.make_species_list()

print()

print("Total species :", len(species))

print("Si =", species.count("Si"))
print("O  =", species.count("O"))
print("Ca =", species.count("Ca"))

print()

print("First 20 atoms:")

print(species[:20])
