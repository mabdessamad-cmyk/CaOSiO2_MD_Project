import numpy as np

from topology import merge_oxygens

oxygen_positions = [

    np.array([0.00, 0.00, 0.00]),
    np.array([0.05, 0.00, 0.00]),

    np.array([5.00, 5.00, 5.00]),
    np.array([5.10, 5.00, 5.00]),

    np.array([9.00, 9.00, 9.00])

]

merged, mapping = merge_oxygens(oxygen_positions)

print("Original :", len(oxygen_positions))
print("Merged   :", len(merged))
print("Mapping  :", mapping)

for p in merged:
    print(p)
