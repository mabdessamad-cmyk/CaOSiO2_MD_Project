import numpy as np

from topology import find_si_neighbors

box = 20.0

si_positions = [

    np.array([0.0, 0.0, 0.0]),

    np.array([3.0, 0.0, 0.0]),

    np.array([8.0, 0.0, 0.0]),

    np.array([0.0, 8.0, 0.0]),

]

neighbors = find_si_neighbors(
    si_positions,
    box,
    cutoff=4.0
)

for i, n in enumerate(neighbors):
    print(f"Si {i}: {n}")
