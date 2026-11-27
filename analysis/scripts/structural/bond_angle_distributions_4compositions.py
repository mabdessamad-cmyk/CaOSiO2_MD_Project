"""
Composition-resolved bond-angle distributions for CaO-SiO2 glasses.

Uses:
    - Validated LAMMPS reader
    - Existing periodic NeighborList
    - Adopted Si-O and Ca-O cutoffs
    - Minimum-image angle calculation

Outputs:
    analysis/figures/osi_o_distribution_4compositions.png
    analysis/figures/si_o_si_BO_distribution_4compositions.png
    analysis/reports/bond_angle_distributions_4compositions.csv
    analysis/reports/bond_angle_histograms_4compositions.csv

Existing project files are never overwritten.
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.bond_angles import BondAngleAnalysis
from workflow.parameters import BOND_CUTOFFS


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

STRUCTURES = {
    0.10: "simulations/data/x010_glass_final_300K.data",
    0.20: "simulations/data/x020_glass_final_300K.data",
    0.30: "simulations/data/glass_final_300K.data",
    0.40: "simulations/data/x040_glass_final_300K.data",
}

OUTPUTS = [
    Path("analysis/figures/osi_o_distribution_4compositions.png"),
    Path("analysis/figures/si_o_si_BO_distribution_4compositions.png"),
    Path("analysis/reports/bond_angle_distributions_4compositions.csv"),
    Path("analysis/reports/bond_angle_histograms_4compositions.csv"),
]

N_BINS = 360
ANGLE_MIN = 0.0
ANGLE_MAX = 180.0
SMOOTHING_SIGMA = 2.0

COLORS = {
    0.10: "tab:blue",
    0.20: "tab:orange",
    0.30: "tab:green",
    0.40: "tab:red",
}


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------

def calculate_statistics(values):
    values = np.asarray(values, dtype=float)

    if len(values) == 0:
        return {
            "number": 0,
            "mean_deg": np.nan,
            "std_deg": np.nan,
            "min_deg": np.nan,
            "max_deg": np.nan,
        }

    return {
        "number": int(len(values)),
        "mean_deg": float(np.mean(values)),
        "std_deg": float(np.std(values)),
        "min_deg": float(np.min(values)),
        "max_deg": float(np.max(values)),
    }


def calculate_bo_only_siosi(system, neighbors, angle_analysis):
    """
    Calculate Si-O-Si angles only at bridging oxygen atoms.

    A bridging oxygen is defined here as an oxygen with exactly
    two Si neighbours under the adopted Si-O cutoff.
    """

    angles = []
    bo_count = 0

    oxygen_indices = np.where(system.types == 2)[0]

    for oxygen in oxygen_indices:
        si_neighbors = [
            atom
            for atom in neighbors.neighbors_of(oxygen)
            if system.types[atom] == 1
        ]

        if len(si_neighbors) != 2:
            continue

        bo_count += 1

        si1, si2 = si_neighbors

        angles.append(
            angle_analysis.angle(si1, oxygen, si2)
        )

    return np.asarray(angles, dtype=float), bo_count


def calculate_histogram(values):
    counts, edges = np.histogram(
        values,
        bins=N_BINS,
        range=(ANGLE_MIN, ANGLE_MAX),
        density=True,
    )

    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = gaussian_filter1d(counts, sigma=SMOOTHING_SIGMA)

    # Preserve unit area after smoothing.
    bin_width = edges[1] - edges[0]
    area = np.sum(smoothed) * bin_width

    if area > 0:
        smoothed /= area

    return centers, smoothed


# ------------------------------------------------------------
# Main analysis
# ------------------------------------------------------------

def main():

    # Create output directories if needed.
    Path("analysis/figures").mkdir(parents=True, exist_ok=True)
    Path("analysis/reports").mkdir(parents=True, exist_ok=True)

    # Protect all existing files.
    existing = [str(path) for path in OUTPUTS if path.exists()]

    if existing:
        raise FileExistsError(
            "The following output files already exist. "
            "No files were modified:\n" + "\n".join(existing)
        )

    summary_rows = []
    histogram_rows = []

    distributions = {
        "O-Si-O": {},
        "Si-O-Si (BO only)": {},
    }

    for composition, filename in STRUCTURES.items():

        print("\n" + "=" * 65)
        print(f"Composition: x(CaO) = {composition:.2f}")
        print(f"Structure: {filename}")

        if not Path(filename).is_file():
            raise FileNotFoundError(filename)

        system = read_lammps(filename)

        neighbors = NeighborList(system, BOND_CUTOFFS)
        neighbors.build()

        angle_analysis = BondAngleAnalysis(system, neighbors)

        # O-Si-O: all pairs of O neighbours around Si.
        osi_o = angle_analysis.angles_for_triplet(
            outer_type=2,
            center_type=1,
        )

        # Si-O-Si: unrestricted, for comparison with existing results.
        si_o_si_all = angle_analysis.angles_for_triplet(
            outer_type=1,
            center_type=2,
        )

        # Si-O-Si: bridging oxygen only, exactly two Si neighbours.
        si_o_si_bo, bo_count = calculate_bo_only_siosi(
            system,
            neighbors,
            angle_analysis,
        )

        print(f"O-Si-O angles: {len(osi_o)}")
        print(f"Unrestricted Si-O-Si angles: {len(si_o_si_all)}")
        print(f"BO atoms (exactly 2 Si): {bo_count}")
        print(f"BO-only Si-O-Si angles: {len(si_o_si_bo)}")

        datasets = {
            "O-Si-O": osi_o,
            "Si-O-Si (all O with >=2 Si)": si_o_si_all,
            "Si-O-Si (BO only)": si_o_si_bo,
        }

        for angle_name, values in datasets.items():

            stats = calculate_statistics(values)

            summary_rows.append({
                "x_CaO": composition,
                "angle_type": angle_name,
                **stats,
            })

            print(
                f"{angle_name}: "
                f"N={stats['number']}, "
                f"mean={stats['mean_deg']:.4f} deg, "
                f"SD={stats['std_deg']:.4f} deg"
            )

        # Store the two distributions intended for plotting.
        distributions["O-Si-O"][composition] = osi_o
        distributions["Si-O-Si (BO only)"][composition] = si_o_si_bo

        for angle_name, values in [
            ("O-Si-O", osi_o),
            ("Si-O-Si (BO only)", si_o_si_bo),
        ]:

            centers, density = calculate_histogram(values)

            for center, probability_density in zip(centers, density):
                histogram_rows.append({
                    "x_CaO": composition,
                    "angle_type": angle_name,
                    "angle_deg": float(center),
                    "probability_density_per_degree": float(
                        probability_density
                    ),
                })

    # --------------------------------------------------------
    # Plot O-Si-O distributions
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for composition, values in distributions["O-Si-O"].items():

        centers, density = calculate_histogram(values)

        ax.plot(
            centers,
            density,
            linewidth=2.2,
            color=COLORS[composition],
            label=fr"$x={composition:.2f}$",
        )

    ax.set_xlabel(r"O--Si--O angle ($^\circ$)", fontsize=14)
    ax.set_ylabel("Probability density (degree$^{-1}$)", fontsize=13)
    ax.set_xlim(85, 135)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=11)
    fig.tight_layout()

    fig.savefig(OUTPUTS[0], dpi=600)
    plt.close(fig)

    # --------------------------------------------------------
    # Plot BO-only Si-O-Si distributions
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for composition, values in distributions["Si-O-Si (BO only)"].items():

        centers, density = calculate_histogram(values)

        ax.plot(
            centers,
            density,
            linewidth=2.2,
            color=COLORS[composition],
            label=fr"$x={composition:.2f}$",
        )

    ax.set_xlabel(r"Si--O--Si angle ($^\circ$)", fontsize=14)
    ax.set_ylabel("Probability density (degree$^{-1}$)", fontsize=13)
    ax.set_xlim(100, 180)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=11)
    fig.tight_layout()

    fig.savefig(OUTPUTS[1], dpi=600)
    plt.close(fig)

    # --------------------------------------------------------
    # Save numerical summary
    # --------------------------------------------------------

    with open(OUTPUTS[2], "x", newline="") as f:
        fieldnames = [
            "x_CaO",
            "angle_type",
            "number",
            "mean_deg",
            "std_deg",
            "min_deg",
            "max_deg",
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    # --------------------------------------------------------
    # Save histogram data
    # --------------------------------------------------------

    with open(OUTPUTS[3], "x", newline="") as f:
        fieldnames = [
            "x_CaO",
            "angle_type",
            "angle_deg",
            "probability_density_per_degree",
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(histogram_rows)

    print("\n" + "=" * 65)
    print("COMPOSITION-RESOLVED ANGLE ANALYSIS COMPLETED")
    print("=" * 65)

    for path in OUTPUTS:
        print(f"Saved: {path}")

    print("\nNo existing project files were overwritten.")


if __name__ == "__main__":
    main()
