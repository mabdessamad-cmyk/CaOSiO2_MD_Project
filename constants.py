"""
constants.py

Physical constants and Pedone (2006) parameters
for CaO–SiO2 glass generation.
"""

# Atomic masses (g/mol)
MASS = {
    "Si": 28.0855,
    "O": 15.999,
    "Ca": 40.078
}

# Formal charges
CHARGE = {
    "Si": 2.4,
    "O": -1.2,
    "Ca": 1.2
}

# LAMMPS atom types
TYPE = {
    "Si": 1,
    "O": 2,
    "Ca": 3
}

# Minimum allowed distances (Å)
CUTOFF = {
    ("Si", "Si"): 3.00,
    ("Si", "O"): 1.50,
    ("Si", "Ca"): 2.80,
    ("O", "O"): 2.20,
    ("O", "Ca"): 2.00,
    ("Ca", "Ca"): 3.00
}


def cutoff(species1, species2):
    """Return the cutoff distance for any pair of species."""

    if (species1, species2) in CUTOFF:
        return CUTOFF[(species1, species2)]

    return CUTOFF[(species2, species1)]
