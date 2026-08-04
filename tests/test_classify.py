from topology import classify_oxygens

owners = [
    [0, 1],
    [2, 3],
    [4],
    [5],
    [6, 7]
]

bo, nbo = classify_oxygens(owners)

print("BO :", bo)
print("NBO:", nbo)
