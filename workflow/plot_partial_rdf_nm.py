#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Input RDF files
# ----------------------------------------------------------------------

BASE = Path(__file__).resolve().parents[1] / "simulations" / "data"

RDF_FILES = {
    0.10: BASE / "x010_partial_rdf_corrected.dat",
    0.20: BASE / "x020_partial_rdf_corrected.dat",
    0.30: BASE / "x030_partial_rdf_corrected.dat",
    0.40: BASE / "x040_partial_rdf_corrected.dat",
}


# ----------------------------------------------------------------------
# LAMMPS RDF column mapping
#
# Column 1 in the file = bin index
# Column 2 in the file = r
#
# Then each pair contributes:
#     g(r), coordination(r)
#
# Therefore, in zero-based Python indexing:
#     Si-Si -> column 2
#     Si-O  -> column 4
#     Si-Ca -> column 6
#     O-O   -> column 8
#     O-Ca -> column 10
#     Ca-Ca -> column 12
# ----------------------------------------------------------------------

PAIR_G_COLUMNS = {
    "Si-Si": 2,
    "Si-O": 4,
    "Si-Ca": 6,
    "O-O": 8,
    "O-Ca": 10,
    "Ca-Ca": 12,
}


# ----------------------------------------------------------------------
# Read and average the two RDF blocks
# ----------------------------------------------------------------------

def read_rdf_blocks(filename):
    """
    Read a LAMMPS fix ave/time RDF file containing two 300-bin blocks.

    Returns
    -------
    r_angstrom : ndarray, shape (300,)
        Radial coordinate in Angstrom.
    rdf : ndarray, shape (300, 6)
        Pointwise average of the two RDF blocks.
    """

    blocks = []
    current_block = []

    with filename.open("r") as fh:
        for line in fh:
            fields = line.split()

            # Skip comment lines and blank lines.
            if not fields or fields[0].startswith("#"):
                continue

            # A block header has two fields:
            # timestep, number-of-rows
            if len(fields) == 2:
                if current_block:
                    blocks.append(np.asarray(current_block, dtype=float))
                    current_block = []
                continue

            # RDF data rows contain 14 fields.
            if len(fields) == 14:
                current_block.append([float(x) for x in fields])

    if current_block:
        blocks.append(np.asarray(current_block, dtype=float))

    if len(blocks) != 2:
        raise ValueError(
            f"{filename}: expected 2 RDF blocks, found {len(blocks)}"
        )

    for i, block in enumerate(blocks, start=1):
        if block.shape != (300, 14):
            raise ValueError(
                f"{filename}: block {i} has shape {block.shape}; "
                "expected (300, 14)"
            )

    # Check that both blocks use exactly the same radial grid.
    if not np.allclose(blocks[0][:, 1], blocks[1][:, 1]):
        raise ValueError(
            f"{filename}: RDF blocks do not share the same radial grid"
        )

    # Pointwise average of the two 10-ps blocks.
    averaged = 0.5 * (blocks[0] + blocks[1])

    r_angstrom = averaged[:, 1]

    return r_angstrom, averaged


# ----------------------------------------------------------------------
# Load all compositions
# ----------------------------------------------------------------------

rdf_data = {}

for x, filename in RDF_FILES.items():
    if not filename.exists():
        raise FileNotFoundError(f"Missing RDF file: {filename}")

    r_angstrom, averaged = read_rdf_blocks(filename)

    rdf_data[x] = {
        "r_nm": r_angstrom / 10.0,
        "data": averaged,
    }

    print(f"x = {x:.2f}")
    print(f"  file   : {filename}")
    print(f"  bins   : {len(r_angstrom)}")
    print(f"  r range: {r_angstrom[0]:.3f}–{r_angstrom[-1]:.3f} Å")
    print()


# ----------------------------------------------------------------------
# Common plotting settings
# ----------------------------------------------------------------------

composition_labels = {
    0.10: "x = 0.10",
    0.20: "x = 0.20",
    0.30: "x = 0.30",
    0.40: "x = 0.40",
}


def make_figure(pairs, output_name):
    fig, axes = plt.subplots(
        1,
        len(pairs),
        figsize=(12, 3.8),
        sharex=True,
    )

    for ax, pair in zip(axes, pairs):
        g_col = PAIR_G_COLUMNS[pair]

        for x in sorted(rdf_data):
            r_nm = rdf_data[x]["r_nm"]
            g_r = rdf_data[x]["data"][:, g_col]

            ax.plot(
                r_nm,
                g_r,
                linewidth=1.5,
                label=composition_labels[x],
            )

        ax.set_xlabel(r"$r$ (nm)")
        ax.set_ylabel(r"$g_{ij}(r)$")
        ax.set_xlim(0.0, 1.5)
        ax.grid(True, alpha=0.25)

        ax.set_title(pair)

    axes[-1].legend(
        loc="upper right",
        frameon=False,
    )

    fig.tight_layout()

    output_path = BASE.parent / "figures" / output_name
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_path,
        dpi=600,
        bbox_inches="tight",
    )

    print(f"Saved: {output_path}")

    plt.close(fig)


# ----------------------------------------------------------------------
# Figure 6
# ----------------------------------------------------------------------

make_figure(
    ["Si-Si", "Si-O", "O-O"],
    "Figure6_partial_RDF_network_pairs_nm.png",
)


# ----------------------------------------------------------------------
# Figure 7
# ----------------------------------------------------------------------

make_figure(
    ["Si-Ca", "O-Ca", "Ca-Ca"],
    "Figure7_partial_RDF_modifier_pairs_nm.png",
)
