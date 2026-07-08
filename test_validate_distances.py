from generator import GlassGenerator
from geometry import distance
from constants import cutoff

# ------------------------------------------------------------
# Generate a small test system
# ------------------------------------------------------------

g = GlassGenerator(
    x=0.30,
    formula_units=20,
    density=1.50,
    seed=12345,
)

g.place_species("Si")
g.place_species("Ca")
g.place_species("O")

print()
print("Checking minimum-distance constraints...")
print()

natoms = g.cfg.natoms

violations = 0

for i in range(natoms):

    species_i = g.cfg.get_species(i)
    pos_i = g.cfg.get_position(i)

    for j in range(i + 1, natoms):

        species_j = g.cfg.get_species(j)
        pos_j = g.cfg.get_position(j)

        r = distance(
            pos_i,
            pos_j,
            g.box,
        )

        rcut = cutoff(
            species_i,
            species_j,
        )

        if r < rcut:

            violations += 1

            print(
                f"Violation: "
                f"{i} ({species_i}) - "
                f"{j} ({species_j})   "
                f"distance = {r:.4f} Å   "
                f"cutoff = {rcut:.4f} Å"
            )

print()

if violations == 0:
    print("SUCCESS")
    print("No minimum-distance violations detected.")
else:
    print(f"FAILED: {violations} violations detected.")
