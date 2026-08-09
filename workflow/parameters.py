"""
Physical parameters used throughout the CaO-SiO2 analysis workflow.
"""

# ------------------------------------------------------------------
# Atom types
# ------------------------------------------------------------------

ATOM_TYPES = {
    1: "Si",
    2: "O",
    3: "Ca",
}


# ------------------------------------------------------------------
# Bond cutoffs (Angstrom)
# ------------------------------------------------------------------

BOND_CUTOFFS = {
    (1, 2): 2.30,   # Si-O
    (2, 1): 2.30,   # O-Si

    (3, 2): 3.10,   # Ca-O
    (2, 3): 3.10,   # O-Ca
}

# ------------------------------------------------------------------
# Cutoff validation
# ------------------------------------------------------------------

# Si-O cutoff:
# RDF first minimum = 2.325 Å.
# Sensitivity study showed stable Si coordination from 2.10–2.40 Å.
# Adopted value: 2.30 Å.

# Ca-O cutoff:
# RDF first minimum = 3.675 Å.
# Sensitivity study showed strong coordination growth above 3.10 Å.
# Adopted value: 3.10 Å.
