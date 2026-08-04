import numpy as np

from generator import GlassGenerator

g = GlassGenerator()

print("Initially:", g.cfg.natoms)

p = np.array([10.0, 10.0, 10.0])

g.insert_atom("Si", p)

print("After insertion:", g.cfg.natoms)

print("Species:", g.cfg.get_species(0))

print("Position:", g.cfg.get_position(0))

print("Neighbor search:", g.grid.nearby_atoms(p))
