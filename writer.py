"""
writer.py

LAMMPS data-file writer for CaO–SiO₂ glasses.
"""

from constants import MASS

# ------------------------------------------------------------
# LAMMPS atom types
# ------------------------------------------------------------

ATOM_TYPES = {
    "Si": 1,
    "O": 2,
    "Ca": 3,
}

# ------------------------------------------------------------
# Partial charges (Pedone model)
# ------------------------------------------------------------

CHARGES = {
    "Si": 2.4,
    "O": -1.2,
    "Ca": 1.2,
}
def write_header(file, cfg):
    """
    Write the LAMMPS data-file header.

    Parameters
    ----------
    file : file object
        Open output file.

    cfg : GlassConfiguration
        Glass configuration to write.
    """

    file.write("LAMMPS data file\n\n")

    file.write(f"{cfg.natoms} atoms\n")
    file.write(f"{len(ATOM_TYPES)} atom types\n\n")

    file.write(f"0.0 {cfg.box:.6f} xlo xhi\n")
    file.write(f"0.0 {cfg.box:.6f} ylo yhi\n")
    file.write(f"0.0 {cfg.box:.6f} zlo zhi\n\n")

def write_masses(file):
    """
    Write the Masses section.
    """

    file.write("Masses\n\n")

    for species in ("Si", "O", "Ca"):

        atom_type = ATOM_TYPES[species]
        mass = MASS[species]

        file.write(f"{atom_type} {mass:.6f}\n")

    file.write("\n")

def write_atoms(file, cfg):
    """
    Write the Atoms section for atom_style charge.

    Parameters
    ----------
    file : file object
        Open output file.

    cfg : GlassConfiguration
        Configuration to write.
    """

    file.write("Atoms\n\n")

    for i in range(cfg.natoms):

        species = cfg.get_species(i)

        position = cfg.get_position(i)

        atom_type = ATOM_TYPES[species]

        charge = CHARGES[species]

        x, y, z = position

        file.write(
            f"{i + 1:d} "
            f"{atom_type:d} "
            f"{charge:.4f} "
            f"{x:.6f} "
            f"{y:.6f} "
            f"{z:.6f}\n"
        )

    file.write("\n")

def write_lammps_data(cfg, filename):
    """
    Write a complete LAMMPS data file.

    Parameters
    ----------
    cfg : GlassConfiguration
        Configuration to write.

    filename : str
        Output filename.
    """

    with open(filename, "w") as file:

        write_header(file, cfg)

        write_masses(file)

        write_atoms(file, cfg)
