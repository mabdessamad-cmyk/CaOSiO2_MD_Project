import numpy as np

from topology import complete_tetrahedra, validate_coordination

si_positions = [
    np.array([0., 0., 0.]),
    np.array([5., 0., 0.]),
    np.array([10., 0., 0.]),
]

oxygen_positions = [
    np.array([2.5, 0., 0.]),
    np.array([7.5, 0., 0.]),
]

oxygen_owners = [
    [0, 1],
    [1, 2],
]

coordination = validate_coordination(oxygen_owners, 3)

print("Before:", coordination)

oxygen_positions, oxygen_owners = complete_tetrahedra(
    si_positions,
    oxygen_positions,
    oxygen_owners,
    coordination
)

coordination = validate_coordination(oxygen_owners, 3)

print("After :", coordination)

print("Total oxygens:", len(oxygen_positions))
