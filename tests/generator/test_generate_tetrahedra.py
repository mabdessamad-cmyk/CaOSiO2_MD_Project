import numpy as np
from topology import generate_tetrahedra

si_positions = [
    np.array([0.0, 0.0, 0.0]),
    np.array([5.0, 5.0, 5.0])
]

oxygen_positions, oxygen_owner = generate_tetrahedra(si_positions)

print("Number of oxygens :", len(oxygen_positions))
print("Number of owners  :", len(oxygen_owner))
print("Owners            :", oxygen_owner)

for p in oxygen_positions:
    print(p)
