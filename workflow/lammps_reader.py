"""
workflow/lammps_reader.py

Read LAMMPS data files (atom_style charge)

Author:
CaOSiO2 Workflow
"""
import numpy as np
from workflow.system import AtomicSystem


def read_lammps(filename):
    """
    Read a LAMMPS data file (atom_style charge).
    """

    print(f"Reading {filename}")

    with open(filename, "r") as f:
        lines = f.readlines()

    print(f"Lines read : {len(lines)}")

    # --------------------------------------------------------
    # Number of atoms
    # --------------------------------------------------------

    natoms = None

    for line in lines:

        words = line.split()

        if len(words) == 2 and words[1] == "atoms":
            natoms = int(words[0])
            break

    if natoms is None:
        raise RuntimeError("Cannot find number of atoms.")

    print(f"Atoms : {natoms}")

    # --------------------------------------------------------
    # Simulation box
    # --------------------------------------------------------

    box = np.zeros((3, 2))

    for line in lines:

        words = line.split()

        if len(words) == 4:

            if words[2] == "xlo":
                box[0] = [float(words[0]), float(words[1])]

            elif words[2] == "ylo":
                box[1] = [float(words[0]), float(words[1])]

            elif words[2] == "zlo":
                box[2] = [float(words[0]), float(words[1])]

    print("\nSimulation box")
    print(box)

    # --------------------------------------------------------
    # Find Atoms section
    # --------------------------------------------------------

    atoms_start = None

    for i, line in enumerate(lines):

        if line.startswith("Atoms"):

            atoms_start = i + 2
            break

    if atoms_start is None:
        raise RuntimeError("Atoms section not found.")

    # --------------------------------------------------------
    # Allocate arrays
    # --------------------------------------------------------

    ids = np.zeros(natoms, dtype=int)
    types = np.zeros(natoms, dtype=int)
    charges = np.zeros(natoms)
    positions = np.zeros((natoms, 3))

    # --------------------------------------------------------
    # Read atoms
    # --------------------------------------------------------

    for i in range(natoms):

        words = lines[atoms_start + i].split()

        ids[i] = int(words[0])
        types[i] = int(words[1])
        charges[i] = float(words[2])

        positions[i, 0] = float(words[3])
        positions[i, 1] = float(words[4])
        positions[i, 2] = float(words[5])

    print("\nAtom table successfully read.")

    system = AtomicSystem(
        ids=ids,
        types=types,
        charges=charges,
        positions=positions,
        box=box,
    )

    return system
