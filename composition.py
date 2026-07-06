"""
composition.py

Stoichiometry and box-size calculations for
x(CaO)–(1-x)SiO2 glasses.

Single source of truth for:
- atom counts
- molecular weight
- system density
- simulation box size
"""

import math
from constants import MASS

AVOGADRO = 6.02214076e23


def calculate_composition(x, formula_units, density):
    """
    Robust stoichiometric construction for CaO–SiO2 glasses.
    """

    if not (0.0 <= x <= 1.0):
        raise ValueError("x must be between 0 and 1.")

    if formula_units <= 0:
        raise ValueError("formula_units must be positive.")

    total_atoms = (3-x) * formula_units
    n_si = round((1.0 - x) * formula_units)
    n_ca = round(x * formula_units)
    n_o  = round((2.0 - x) * formula_units)

    si = int(n_si)
    ca = int(n_ca)
    o  = int(n_o)
    natoms = si + ca + o

    mw = (
        si * MASS["Si"] +
        ca * MASS["Ca"] +
        o * MASS["O"]
    )

    mass_g = mw / AVOGADRO
    volume_cm3 = mass_g / density
    volume_A3 = volume_cm3 * 1.0e24
    box = volume_A3 ** (1.0 / 3.0)

    return {
        "x": x,
        "formula_units": formula_units,
        "Si": si,
        "Ca": ca,
        "O": o,
        "natoms": natoms,
        "density": density,
        "molecular_weight": mw,
        "box_length": box,
    }

def total_charge(comp):
    """
    Compute total system charge.
    """

    return (
        comp["Si"] * 2.4 +
        comp["Ca"] * 1.2 -
        comp["O"] * 1.2
    )
