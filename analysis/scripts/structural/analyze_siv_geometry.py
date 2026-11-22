#!/usr/bin/env python3
"""
Site-resolved geometrical characterization of five-coordinated Si
in CaO-SiO2 glasses.

Uses the adopted Si-O cutoff, validated NeighborList, QnAnalysis
and minimum-image BondAngleAnalysis.

Outputs:
    analysis/reports/siv_geometry_sites.csv
    analysis/reports/siv_geometry_summary.csv

Existing files are never overwritten.
"""

import csv
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.bond_angles import BondAngleAnalysis
from workflow.qn import QnAnalysis
from workflow.parameters import BOND_CUTOFFS


STRUCTURES = {
    0.10: "simulations/data/x010_glass_final_300K.data",
    0.20: "simulations/data/x020_glass_final_300K.data",
    0.30: "simulations/data/glass_final_300K.data",
    0.40: "simulations/data/x040_glass_final_300K.data",
}

EXPECTED_SIV = {
    0.10: {"total": 8, "BO4_NBO1": 1, "BO5_NBO0": 7},
    0.20: {"total": 11, "BO4_NBO1": 1, "BO5_NBO0": 10},
    0.30: {"total": 5, "BO4_NBO1": 3, "BO5_NBO0": 2},
    0.40: {"total": 7, "BO4_NBO1": 2, "BO5_NBO0": 3},
}

SITE_OUTPUT = Path("analysis/reports/siv_geometry_sites.csv")
SUMMARY_OUTPUT = Path("analysis/reports/siv_geometry_summary.csv")

SI = 1
O = 2


def mean_std(values):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return np.nan, np.nan
    return float(np.mean(values)), float(np.std(values))


def analyze_site(system, neighbors, qn, angles, si_index):
    """Calculate the five bond lengths and ten angles for one SiV site."""

    oxygen_indices = sorted(
        j for j in neighbors.neighbors_of(si_index)
        if system.types[j] == O
    )

    if len(oxygen_indices) != 5:
        raise ValueError(
            f"Si atom index {si_index}: expected 5 O neighbours, "
            f"found {len(oxygen_indices)}"
        )

    distances = []
    ligand_classes = []

    for oxygen_index in oxygen_indices:
        vector = angles._vector(si_index, oxygen_index)
        distances.append(float(np.linalg.norm(vector)))
        ligand_classes.append(qn.oxygen_class(oxygen_index))

    angle_values = [
        angles.angle(i, si_index, j)
        for i, j in combinations(oxygen_indices, 2)
    ]

    if len(angle_values) != 10:
        raise ValueError("A five-coordinate site must have 10 angles.")

    n_bo = ligand_classes.count("BO")
    n_nbo = ligand_classes.count("NBO")
    n_other = ligand_classes.count("other")

    if n_bo == 4 and n_nbo == 1 and n_other == 0:
        subtype = "BO4_NBO1"
    elif n_bo == 5 and n_nbo == 0 and n_other == 0:
        subtype = "BO5_NBO0"
    else:
        subtype = "OTHER"

    return {
        "oxygen_indices": oxygen_indices,
        "distances": distances,
        "angles": angle_values,
        "n_bo": n_bo,
        "n_nbo": n_nbo,
        "n_other": n_other,
        "subtype": subtype,
    }


def write_csv(path, rows, fieldnames):
    with path.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    for path in (SITE_OUTPUT, SUMMARY_OUTPUT):
        if path.exists():
            raise FileExistsError(
                f"Output already exists; refusing to overwrite: {path}"
            )

    site_rows = []
    grouped = defaultdict(list)

    for composition, filename in STRUCTURES.items():

        print("\n" + "=" * 65)
        print(f"Composition x(CaO) = {composition:.2f}")
        print(f"Structure: {filename}")

        if not Path(filename).is_file():
            raise FileNotFoundError(filename)

        system = read_lammps(filename)

        neighbors = NeighborList(system, BOND_CUTOFFS)
        neighbors.build()

        qn = QnAnalysis(system, neighbors)
        angle_analysis = BondAngleAnalysis(system, neighbors)

        si_indices = np.where(system.types == SI)[0]
        siv_indices = []

        for i in si_indices:
            oxygen_neighbors = [
                j for j in neighbors.neighbors_of(i)
                if system.types[j] == O
            ]
            if len(oxygen_neighbors) == 5:
                siv_indices.append(i)

        print(f"Identified SiV sites: {len(siv_indices)}")
        expected = EXPECTED_SIV[composition]

        if len(siv_indices) != expected["total"]:
            raise ValueError(
                f"x={composition:.2f}: expected "
                f"{expected['total']} SiV sites, "
                f"found {len(siv_indices)}"
            )


        for si_index in siv_indices:

            result = analyze_site(
                system, neighbors, qn, angle_analysis, si_index
            )

            atom_id = int(system.ids[si_index])
            distances = result["distances"]
            angle_values = result["angles"]

            bond_mean, bond_std = mean_std(distances)
            angle_mean, angle_std = mean_std(angle_values)

            row = {
                "x_CaO": composition,
                "Si_atom_ID": atom_id,
                "Si_array_index": int(si_index),
                "SiV_subtype": result["subtype"],
                "BO_count": result["n_bo"],
                "NBO_count": result["n_nbo"],
                "other_O_count": result["n_other"],
                "O_atom_IDs": ";".join(
                    str(int(system.ids[j]))
                    for j in result["oxygen_indices"]
                ),
                "Si_O_distances_A": ";".join(
                    f"{v:.6f}" for v in distances
                ),
                "mean_Si_O_A": bond_mean,
                "std_Si_O_A": bond_std,
                "min_Si_O_A": min(distances),
                "max_Si_O_A": max(distances),
                "O_Si_O_angles_deg": ";".join(
                    f"{v:.6f}" for v in angle_values
                ),
                "mean_O_Si_O_deg": angle_mean,
                "std_O_Si_O_deg": angle_std,
                "min_O_Si_O_deg": min(angle_values),
                "max_O_Si_O_deg": max(angle_values),
            }

            site_rows.append(row)
            grouped[(composition, result["subtype"])].append(row)

        subtype_counts = {
            subtype: sum(
                row["SiV_subtype"] == subtype
                for row in site_rows
                if row["x_CaO"] == composition
            )
            for subtype in ("BO4_NBO1", "BO5_NBO0")
        }

        if subtype_counts != {
            "BO4_NBO1": expected["BO4_NBO1"],
            "BO5_NBO0": expected["BO5_NBO0"],
        }:
            raise ValueError(
                f"x={composition:.2f}: subtype mismatch: "
                f"{subtype_counts}"
            )

        print(
            f"Validated: {len(siv_indices)} sites, "
            f"{10 * len(siv_indices)} angles"
        )

    # Summary: independent units are SiV sites, not individual angles.
    summary_rows = []

    for (composition, subtype), rows in sorted(grouped.items()):

        bond_means = [r["mean_Si_O_A"] for r in rows]
        bond_stds = [r["std_Si_O_A"] for r in rows]
        angle_means = [r["mean_O_Si_O_deg"] for r in rows]
        angle_stds = [r["std_O_Si_O_deg"] for r in rows]

        bond_mean, bond_between_site_std = mean_std(bond_means)
        bond_distortion_mean, bond_distortion_std = mean_std(bond_stds)
        angle_mean, angle_between_site_std = mean_std(angle_means)
        angle_distortion_mean, angle_distortion_std = mean_std(angle_stds)

        summary_rows.append({
            "x_CaO": composition,
            "SiV_subtype": subtype,
            "n_SiV_sites": len(rows),
            "mean_site_Si_O_A": bond_mean,
            "between_site_SD_Si_O_A": bond_between_site_std,
            "mean_within_site_SD_Si_O_A": bond_distortion_mean,
            "SD_within_site_SD_Si_O_A": bond_distortion_std,
            "mean_site_O_Si_O_deg": angle_mean,
            "between_site_SD_O_Si_O_deg": angle_between_site_std,
            "mean_within_site_SD_O_Si_O_deg": angle_distortion_mean,
            "SD_within_site_SD_O_Si_O_deg": angle_distortion_std,
        })

    site_fields = list(site_rows[0].keys()) if site_rows else []
    summary_fields = [
        "x_CaO",
        "SiV_subtype",
        "n_SiV_sites",
        "mean_site_Si_O_A",
        "between_site_SD_Si_O_A",
        "mean_within_site_SD_Si_O_A",
        "SD_within_site_SD_Si_O_A",
        "mean_site_O_Si_O_deg",
        "between_site_SD_O_Si_O_deg",
        "mean_within_site_SD_O_Si_O_deg",
        "SD_within_site_SD_O_Si_O_deg",
    ]

    if not site_rows:
        raise RuntimeError("No SiV sites were identified.")

    SITE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    write_csv(SITE_OUTPUT, site_rows, site_fields)
    write_csv(SUMMARY_OUTPUT, summary_rows, summary_fields)

    print("\n" + "=" * 65)
    print("SiV GEOMETRY ANALYSIS COMPLETED")
    print(f"Total SiV sites: {len(site_rows)}")
    print(f"Site-level data: {SITE_OUTPUT}")
    print(f"Summary data:    {SUMMARY_OUTPUT}")


if __name__ == "__main__":
    main()
