import numpy as np

from generator import GlassGenerator

g = GlassGenerator()

# First atom: the box is empty, so any position should be accepted.
p1 = np.array([10.0, 10.0, 10.0])

print("Empty box:")
print(g.check_position("Si", p1))

# Insert the atom manually.
g.cfg.add_atom("Si", p1)
g.grid.add_atom(0, p1)

# Second atom: too close to the first.
p2 = np.array([10.5, 10.0, 10.0])

print("\nToo close:")
print(g.check_position("Si", p2))

# Third atom: far enough away.
p3 = np.array([20.0, 20.0, 20.0])

print("\nFar away:")
print(g.check_position("Si", p3))
