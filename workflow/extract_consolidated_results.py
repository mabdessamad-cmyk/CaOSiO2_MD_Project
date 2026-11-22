#!/usr/bin/env python3

"""
Read-only consolidation of CaO-SiO2 structural results.

Inputs:
    Four final 300 K LAMMPS structures
    Existing periodic_ring_counts.csv

Outputs:
    analysis/reports/consolidated_structural_results.csv
    analysis/reports/consolidated_periodic_ring_details.csv
    analysis/reports/consolidated_structural_validation.txt

Existing outputs are never overwritten.
No simulations, ring enumeration, or figure generation are performed.
"""

import csv
import hashlib
import math
from pathlib import Path

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.report import StructuralReport
from workflow.modifiers import ModifierAnalysis


ROOT = Path(__file__).resolve().parents[1]

STRUCTURES = [
    (0.10, "simulations/data/x010_glass_final_300K.data"),
    (0.20, "simulations/data/x020_glass_final_300K.data"),
    (0.30, "simulations/data/glass_final_300K.data"),
    (0.40, "simulations/data/x040_glass_final_300K.data"),
]

RING_INPUT = ROOT / "analysis/reports/periodic_ring_counts.csv"

OUTPUTS = {
    "structural": ROOT / "analysis/reports/consolidated_structural_results.csv",
    "rings": ROOT / "analysis/reports/consolidated_periodic_ring_details.csv",
    "validation": ROOT / "analysis/reports/consolidated_structural_validation.txt",
}

RING_LIMITS = (6, 8, 10)


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def close(a, b, tolerance=1e-6):
    return math.isclose(
        float(a),
        float(b),
        rel_tol=tolerance,
        abs_tol=tolerance,
    )


def read_ring_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)

        required = {
            "x_CaO",
            "n_Si",
            "max_ring_size",
            "ring_size_Si",
            "ring_count",
            "rings_per_Si",
            "fraction_of_rings",
            "total_rings",
            "contractible_count",
            "winding_count",
            "invalid_count",
        }

        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                "Periodic-ring CSV is missing required columns."
            )

        return list(reader), reader.fieldnames


def calculate_composition(x, relative_path):
    path = ROOT / relative_path

    if not path.is_file():
        raise FileNotFoundError(path)

    system = read_lammps(str(path))

    neighbors = NeighborList(system, BOND_CUTOFFS)
    neighbors.build()

    report = StructuralReport(system, neighbors)
    result = report.results()

    modifier = ModifierAnalysis(system, neighbors)

    composition = result["composition"]
    density = result["density"]
    coordination = result["coordination"]
    qn = result["qn"]
    topology = result["topology"]
    angles = result["bond_angles"]

    n_si = composition["Si"]
    n_o = composition["O"]
    n_ca = composition["Ca"]
    n_atoms = n_si + n_o + n_ca

    # The Qn distribution is defined over SiIV sites only.
    qn_distribution = qn["distribution"]
    siiv_count = sum(qn_distribution.values())
    siv_count = qn["SiV"]["total"]

    other_si_count = n_si - siiv_count - siv_count

    oxygen = qn["oxygen"]

    # Independent Ca modifier statistics.
    ca_environment = modifier.global_environment()
    ca_statistics = modifier.coordination_statistics()

    row = {
        "x_CaO": x,
        "x_SiO2": 1.0 - x,
        "n_Si": n_si,
        "n_O": n_o,
        "n_Ca": n_ca,
        "n_atoms": n_atoms,
        "density_g_cm3": density["density_g_cm3"],

        "mean_Si_O_coord": coordination["Si-O"]["mean"],
        "std_Si_O_coord": coordination["Si-O"]["std"],
        "mean_O_Si_coord": coordination["O-Si"]["mean"],
        "std_O_Si_coord": coordination["O-Si"]["std"],
        "mean_Ca_O_coord": coordination["Ca-O"]["mean"],
        "std_Ca_O_coord": coordination["Ca-O"]["std"],
        "mean_O_Ca_coord": coordination["O-Ca"]["mean"],
        "std_O_Ca_coord": coordination["O-Ca"]["std"],

        "BO_count": oxygen["BO"],
        "NBO_count": oxygen["NBO"],
        "other_O_count": oxygen["other"],
        "BO_fraction": oxygen["BO_fraction"],
        "NBO_fraction": oxygen["NBO_fraction"],
        "other_O_fraction": oxygen["other_fraction"],

        "SiIV_count": siiv_count,
        "SiV_count": siv_count,
        "other_Si_coordination_count": other_si_count,

        "SiV_BO4_NBO1_count": qn["SiV"]["BO4_NBO1"],
        "SiV_BO5_NBO0_count": qn["SiV"]["BO5_NBO0"],

        "Ca_O_BO_contacts": ca_environment["BO"],
        "Ca_O_NBO_contacts": ca_environment["NBO"],
        "Ca_O_other_contacts": ca_environment["other"],

        "mean_Ca_O_contacts_per_Ca": ca_statistics["mean_ca_o"],
        "std_Ca_O_contacts_per_Ca": ca_statistics["std_ca_o"],
        "mean_BO_contacts_per_Ca": ca_statistics["mean_bo"],
        "mean_NBO_contacts_per_Ca": ca_statistics["mean_nbo"],
        "mean_other_contacts_per_Ca": ca_statistics["mean_other"],

        "network_mean_degree": topology["network"]["mean_degree"],
        "network_std_degree": topology["network"]["std_degree"],
        "network_min_degree": topology["network"]["minimum_degree"],
        "network_max_degree": topology["network"]["maximum_degree"],

        "number_components": topology["components"]["number_components"],
        "largest_component_size": topology["components"]["largest_component"],
        "smallest_component_size": topology["components"]["smallest_component"],
        "mean_component_size": topology["components"]["mean_component"],
        "isolated_Si": topology["components"]["isolated_si"],
        "percolating_fraction": topology["components"]["percolating_fraction"],
    }

    for q in range(5):
        row[f"Q{q}_count"] = qn_distribution.get(q, 0)

    for q in range(5):
        row[f"Q{q}_fraction"] = (
            qn_distribution.get(q, 0) / siiv_count
            if siiv_count else float("nan")
        )

    for name, stats in angles.items():
        key = name.replace("-", "_")
        row[f"{key}_number"] = stats["number"]
        row[f"{key}_mean_deg"] = stats["mean"]
        row[f"{key}_std_deg"] = stats["std"]
        row[f"{key}_min_deg"] = stats["min"]
        row[f"{key}_max_deg"] = stats["max"]

    return row, path, system


def validate_ring_data(ring_rows, structural_rows):
    errors = []
    ring_summary = {}

    expected_x = {round(row["x_CaO"], 2): row for row in structural_rows}

    for row in ring_rows:
        x = round(float(row["x_CaO"]), 2)
        limit = int(row["max_ring_size"])
        ring_size = int(row["ring_size_Si"])

        if x not in expected_x:
            errors.append(f"Unexpected ring composition: x={x}")
            continue

        structure = expected_x[x]

        if int(row["n_Si"]) != structure["n_Si"]:
            errors.append(f"x={x}: ring CSV Si count does not match structure")

        if limit not in RING_LIMITS:
            errors.append(f"x={x}: unexpected maximum ring size {limit}")

        if ring_size < 3 or ring_size > limit:
            errors.append(
                f"x={x}, max={limit}: invalid ring size {ring_size}"
            )

        key = (x, limit)

        if key not in ring_summary:
            ring_summary[key] = {
                "counts": [],
                "total": int(row["total_rings"]),
                "contractible": int(row["contractible_count"]),
                "winding": int(row["winding_count"]),
                "invalid": int(row["invalid_count"]),
            }

        entry = ring_summary[key]

        if int(row["total_rings"]) != entry["total"]:
            errors.append(f"x={x}, max={limit}: inconsistent total_rings")

        for field in ("contractible_count", "winding_count", "invalid_count"):
            if int(row[field]) != entry[field.replace("_count", "")]:
                errors.append(
                    f"x={x}, max={limit}: inconsistent {field}"
                )

        entry["counts"].append(int(row["ring_count"]))

    for x in expected_x:
        for limit in RING_LIMITS:
            key = (x, limit)

            if key not in ring_summary:
                errors.append(f"Missing ring data for x={x}, max={limit}")
                continue

            entry = ring_summary[key]

            if sum(entry["counts"]) != entry["total"]:
                errors.append(
                    f"x={x}, max={limit}: ring-size counts do not sum "
                    "to total_rings"
                )

            if (
                entry["contractible"]
                + entry["winding"]
                + entry["invalid"]
                != entry["total"]
            ):
                errors.append(
                    f"x={x}, max={limit}: ring classifications "
                    "do not sum to total_rings"
                )

    return ring_summary, errors


def validate_structural_rows(rows):
    errors = []

    for row in rows:
        x = row["x_CaO"]

        if row["BO_count"] + row["NBO_count"] + row["other_O_count"] != row["n_O"]:
            errors.append(f"x={x}: oxygen classifications do not sum to n_O")

        qn_total = sum(row[f"Q{q}_count"] for q in range(5))

        if qn_total != row["SiIV_count"]:
            errors.append(f"x={x}: Qn counts do not sum to SiIV_count")

        if (
            row["SiIV_count"]
            + row["SiV_count"]
            + row["other_Si_coordination_count"]
            != row["n_Si"]
        ):
            errors.append(f"x={x}: silicon coordination populations do not sum to n_Si")

        ca_contacts = (
            row["Ca_O_BO_contacts"]
            + row["Ca_O_NBO_contacts"]
            + row["Ca_O_other_contacts"]
        )

        if ca_contacts != row["mean_Ca_O_contacts_per_Ca"] * row["n_Ca"]:
            if not close(ca_contacts, row["mean_Ca_O_contacts_per_Ca"] * row["n_Ca"]):
                errors.append(f"x={x}: Ca-O contact totals are inconsistent")

        if row["n_Ca"] + row["n_Si"] == 0:
            errors.append(f"x={x}: invalid composition")

        calculated_x = row["n_Ca"] / (row["n_Ca"] + row["n_Si"])

        if not close(calculated_x, x):
            errors.append(f"x={x}: composition does not match atom counts")

    return errors


def main():
    # Check every destination before performing calculations.
    existing = [str(path) for path in OUTPUTS.values() if path.exists()]

    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing output(s):\n"
            + "\n".join(existing)
        )

    if not RING_INPUT.is_file():
        raise FileNotFoundError(RING_INPUT)

    ring_rows, ring_fields = read_ring_csv(RING_INPUT)

    structural_rows = []
    input_paths = []
    validation_lines = []

    for x, relative_path in STRUCTURES:
        print(f"Processing x={x:.2f}: {relative_path}")

        row, path, system = calculate_composition(x, relative_path)

        structural_rows.append(row)
        input_paths.append(path)

        validation_lines.append(
            f"x={x:.2f}: atoms={row['n_atoms']}, "
            f"Si={row['n_Si']}, O={row['n_O']}, Ca={row['n_Ca']}, "
            f"density={row['density_g_cm3']:.6f} g/cm3"
        )

    structural_errors = validate_structural_rows(structural_rows)

    ring_summary, ring_errors = validate_ring_data(
        ring_rows,
        structural_rows,
    )

    errors = structural_errors + ring_errors

    if errors:
        print("\nVALIDATION FAILED:")
        for error in errors:
            print(" -", error)

        raise RuntimeError(
            "No output files have been written. Resolve validation errors first."
        )

    # Add periodic-ring summary columns to the main structural table.
    for row in structural_rows:
        x = round(row["x_CaO"], 2)

        for limit in RING_LIMITS:
            entry = ring_summary[(x, limit)]

            row[f"ring_paths_max{limit}"] = entry["total"]
            row[f"ring_paths_per_Si_max{limit}"] = (
                entry["total"] / row["n_Si"]
            )
            row[f"ring_contractible_max{limit}"] = entry["contractible"]
            row[f"ring_winding_max{limit}"] = entry["winding"]
            row[f"ring_invalid_max{limit}"] = entry["invalid"]

    structural_fields = list(structural_rows[0].keys())

    validation_lines.extend([
        "",
        "STRUCTURAL VALIDATION",
        "All structural consistency checks passed.",
        "Qn distributions are normalized over SiIV sites.",
        "SiV sites are reported separately from SiIV/Qn.",
        "Oxygen populations are BO, NBO, and other.",
        "Ca-O environment values are contact counts.",
        "",
        "PERIODIC RING VALIDATION",
        "Existing periodic_ring_counts.csv was read without modification.",
        "All ring-size counts and classifications passed consistency checks.",
        "Ring statistics refer to oxygen-labelled primitive ring paths.",
        "No ring enumeration was performed by this script.",
        "",
        "INPUT FILE SHA256",
    ])

    for path in input_paths + [RING_INPUT]:
        validation_lines.append(f"{path.relative_to(ROOT)}: {sha256(path)}")

    validation_lines.extend([
        "",
        "COMPOSITION SUMMARY",
        *validation_lines[:4],
    ])

    # Prepare all outputs in memory before creating any destination.
    structural_output = []
    structural_output.append(",".join(structural_fields))

    for row in structural_rows:
        structural_output.append(
            ",".join(str(row[field]) for field in structural_fields)
        )

    ring_output = [",".join(ring_fields)]

    for row in ring_rows:
        ring_output.append(
            ",".join(row[field] for field in ring_fields)
        )

    validation_output = "\n".join(validation_lines) + "\n"

    # Recheck destinations immediately before writing.
    existing = [str(path) for path in OUTPUTS.values() if path.exists()]

    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing output(s):\n"
            + "\n".join(existing)
        )

    # Exclusive creation: 'x' fails if a file already exists.
    with OUTPUTS["structural"].open("x") as handle:
        handle.write("\n".join(structural_output) + "\n")

    with OUTPUTS["rings"].open("x") as handle:
        handle.write("\n".join(ring_output) + "\n")

    with OUTPUTS["validation"].open("x") as handle:
        handle.write(validation_output)

    print("\nSUCCESS: all validations passed.")
    print("\nCreated files:")

    for path in OUTPUTS.values():
        print(f"  {path.relative_to(ROOT)}")

    print("\nComposition summary:")
    for row in structural_rows:
        print(
            f"x={row['x_CaO']:.2f}: "
            f"density={row['density_g_cm3']:.4f} g/cm3, "
            f"SiIV={row['SiIV_count']}, "
            f"SiV={row['SiV_count']}, "
            f"BO={row['BO_count']}, "
            f"NBO={row['NBO_count']}, "
            f"components={row['number_components']}"
        )


if __name__ == "__main__":
    main()
