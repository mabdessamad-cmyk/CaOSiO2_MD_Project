import numpy as np

from topology import merge_oxygens

oxygen_positions = [

    np.array([0.00, 0.00, 0.00]),
    np.array([0.05, 0.00, 0.00]),

    np.array([5.00, 5.00, 5.00]),
    np.array([5.10, 5.00, 5.00]),

    np.array([9.00, 9.00, 9.00])

]

oxygen_owner = [
    0,  # O1 belongs to Si 0
    1,  # O2 belongs to Si 1 (will merge with O1)
    2,  # O3 belongs to Si 2
    3,  # O4 belongs to Si 3 (will merge with O3)
    4   # O5 belongs to Si 4
]

merged, mapping, owners = merge_oxygens(
    oxygen_positions,
    oxygen_owner
)

print("Original :", len(oxygen_positions))
print("Merged   :", len(merged))
print("Mapping  :", mapping)
print("Owners   :", owners)
for p in merged:
    print(p)
