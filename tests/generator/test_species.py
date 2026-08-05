from generator import GlassGenerator

g = GlassGenerator(
    x=0.30,
    formula_units=2400,
    density=1.50,
    seed=12345,
)

print("Species list sizes")
print("------------------")

print("Si =", len(g.species_lists["Si"]))
print("Ca =", len(g.species_lists["Ca"]))
print("O  =", len(g.species_lists["O"]))

print()

print("First 10 Si atoms:")
print(g.species_lists["Si"][:10])

print()

print("First 10 Ca atoms:")
print(g.species_lists["Ca"][:10])

print()

print("First 10 O atoms:")
print(g.species_lists["O"][:10])
