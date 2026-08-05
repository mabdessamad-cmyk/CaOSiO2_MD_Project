import numpy as np

from topology import connect_tetrahedra

box = 20.0

si_positions = [
    np.array([0., 0., 0.]),
    np.array([3., 0., 0.]),
    np.array([8., 0., 0.]),
]

neighbors = [
    [1],
    [0, 2],
    [1],
]

oxygen_positions, oxygen_owners = connect_tetrahedra(
    si_positions,
    neighbors,
    box
)

print("Number of BO:", len(oxygen_positions))

for p, o in zip(oxygen_positions, oxygen_owners):
    print(o, p)
