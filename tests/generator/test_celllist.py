import numpy as np

from celllist import NeighborGrid

box = 53.4
cell = 3.0

grid = NeighborGrid(box, cell)

positions = [
    np.array([1.0, 1.0, 1.0]),
    np.array([2.5, 1.2, 1.4]),
    np.array([10.0, 10.0, 10.0]),
    np.array([52.9, 1.1, 1.2]),   # near periodic boundary
]

for i, p in enumerate(positions):
    grid.add_atom(i, p)

print("Number of occupied cells:", len(grid.grid))

query = np.array([0.2, 1.0, 1.0])

neighbors = grid.nearby_atoms(query)

print("Neighbors:", sorted(neighbors))

expected = [0, 1, 3]
assert sorted(neighbors) == expected
print("NeighborGrid test PASSED")
