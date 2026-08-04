from topology import validate_coordination

oxygen_owners = [

    [0, 1],
    [1, 2],
    [2],
    [0]

]

coord = validate_coordination(oxygen_owners, 3)

print(coord)
