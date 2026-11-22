
#!/usr/bin/env python3
"""
Create a versioned extension of the consolidated CaO-SiO2 structural results.

Inputs (read-only):
    analysis/reports/consolidated_structural_results.csv
    analysis/reports/siv_geometry_sites.csv
    Four final 300 K structures
    analysis/reports/periodic_ring_counts.csv

Outputs (new files only):
    analysis/reports/consolidated_structural_results_v2.csv
    analysis/reports/consolidated_structural_validation_v2.txt

No simulations, structural recalculations, or ring enumeration are performed.
Existing files are never overwritten.
"""

from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import io
import platform
import sys


ROOT = Path(__file__).resolve().parents[1]

REPORT_DIR = ROOT / "analysis" / "reports"

STRUCTURAL_CSV = REPORT_DIR / "consolidated_structural_results.csv"
SIV_SITES_CSV = REPORT_DIR / "siv_geometry_sites.csv"

RING_CSV = REPORT_DIR / "periodic_ring_counts.csv"

OUTPUT_CSV = REPORT_DIR / "consolidated_structural_results_v2.csv"
OUTPUT_VALIDATION = (
    REPORT_DIR / "consolidated_structural_validation_v2.txt"
)

EXPECTED_X = (0.10, 0.20, 0.30, 0.40)

STRUCTURE_FILES = [
    ROOT / "simulations/data/x010_glass_final_300K.data",
    ROOT / "simulations/data/x020_glass_final_300K.data",
    ROOT / "simulations/data/glass_final_300K.data",
    ROOT / "simulations/data/x040_glass_final_300K.data",
]

PROVENANCE_SCRIPTS = [
    ROOT / "workflow/extract_consolidated_results.py",
    ROOT / "analysis/scripts/structural/analyze_siv_geometry.py",
]


def composition_key(value):
    """Normalize composition values to two decimal places."""
    return round(float(value), 2)


def sha256_file(path):
    """Calculate the SHA256 hash of a file."""
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def read_csv(path):
    """Read a CSV file without modifying it."""
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames:
            raise ValueError(f"Missing CSV header: {path}")

        return reader.fieldnames, list(reader)


def require_columns(path, fieldnames, required):
    """Check that all required CSV columns exist."""
    missing = sorted(set(required) - set(fieldnames))

    if missing:
        raise ValueError(
            f"{path}: missing required columns: {missing}"
        )


def classify_siv_site(row):
    """
    Classify a five-coordinated Si site from its actual ligand counts.

    The original OTHER label is retained as source information.
    BO3_NBO2 is assigned only when the measured counts are exactly
    BO=3, NBO=2, other=0.
    """

    bo = int(row["BO_count"])
    nbo = int(row["NBO_count"])
    other = int(row["other_O_count"])

    if bo + nbo + other != 5:
        raise ValueError(
            f"SiV site {row['Si_atom_ID']} does not have "
            f"five classified oxygen neighbors."
        )

    original_label = row["SiV_subtype"]

    if (bo, nbo, other) == (4, 1, 0):
        if original_label != "BO4_NBO1":
            raise ValueError(
                f"Inconsistent original subtype for Si "
                f"{row['Si_atom_ID']}: {original_label}"
            )
        return "BO4_NBO1"

    if (bo, nbo, other) == (5, 0, 0):
        if original_label != "BO5_NBO0":
            raise ValueError(
                f"Inconsistent original subtype for Si "
                f"{row['Si_atom_ID']}: {original_label}"
            )
        return "BO5_NBO0"

    if (bo, nbo, other) == (3, 2, 0):
        if original_label != "OTHER":
            raise ValueError(
                f"Expected original OTHER label for Si "
                f"{row['Si_atom_ID']}, got {original_label}"
            )
        return "BO3_NBO2"

    if original_label != "OTHER":
        raise ValueError(
            f"Unexpected subtype/count combination for Si "
            f"{row['Si_atom_ID']}: {original_label}, "
            f"BO={bo}, NBO={nbo}, other={other}"
        )

    return "UNCLASSIFIED"


def build_results():
    """Validate and construct the extended structural dataset."""

    structural_fields, structural_rows = read_csv(STRUCTURAL_CSV)
    site_fields, site_rows = read_csv(SIV_SITES_CSV)

    require_columns(
        STRUCTURAL_CSV,
        structural_fields,
        [
            "x_CaO",
            "n_Si",
            "SiV_count",
            "SiV_BO4_NBO1_count",
            "SiV_BO5_NBO0_count",
        ],
    )

    require_columns(
        SIV_SITES_CSV,
        site_fields,
        [
            "x_CaO",
            "Si_atom_ID",
            "SiV_subtype",
            "BO_count",
            "NBO_count",
            "other_O_count",
        ],
    )

    # Check the structural CSV composition coverage.
    structural_by_x = {}

    for row in structural_rows:
        x = composition_key(row["x_CaO"])

        if x in structural_by_x:
            raise ValueError(
                f"Duplicate structural CSV composition: {x}"
            )

        structural_by_x[x] = row

    if set(structural_by_x) != set(EXPECTED_X):
        raise ValueError(
            "Unexpected structural CSV compositions: "
            f"{sorted(structural_by_x)}"
        )

    # Group site-resolved SiV records by composition.
    sites_by_x = defaultdict(list)

    for row in site_rows:
        x = composition_key(row["x_CaO"])

        if x not in EXPECTED_X:
            raise ValueError(
                f"Unexpected SiV site composition: {x}"
            )

        sites_by_x[x].append(row)

    if set(sites_by_x) != set(EXPECTED_X):
        raise ValueError(
            "Missing or unexpected SiV site compositions."
        )

    new_rows = []

    summary = {}

    for x in EXPECTED_X:
        structural_row = structural_by_x[x]
        sites = sites_by_x[x]

        expected_siv = int(structural_row["SiV_count"])

        if len(sites) != expected_siv:
            raise ValueError(
                f"x={x:.2f}: site-level SiV records "
                f"({len(sites)}) do not match "
                f"structural SiV_count ({expected_siv})."
            )

        # Check uniqueness of Si atom IDs within each composition.
        site_ids = [row["Si_atom_ID"] for row in sites]

        if len(site_ids) != len(set(site_ids)):
            raise ValueError(
                f"x={x:.2f}: duplicate Si atom IDs in site CSV."
            )

        counts = {
            "BO4_NBO1": 0,
            "BO5_NBO0": 0,
            "BO3_NBO2": 0,
            "UNCLASSIFIED": 0,
        }

        for site in sites:
            subtype = classify_siv_site(site)
            counts[subtype] += 1

        # Preserve and independently validate original subtype counts.
        original_bo4 = int(
            structural_row["SiV_BO4_NBO1_count"]
        )
        original_bo5 = int(
            structural_row["SiV_BO5_NBO0_count"]
        )

        if counts["BO4_NBO1"] != original_bo4:
            raise ValueError(
                f"x={x:.2f}: BO4_NBO1 count mismatch."
            )

        if counts["BO5_NBO0"] != original_bo5:
            raise ValueError(
                f"x={x:.2f}: BO5_NBO0 count mismatch."
            )

        classified_total = sum(counts.values())

        if classified_total != expected_siv:
            raise ValueError(
                f"x={x:.2f}: SiV subtype counts do not "
                "sum to the total SiV population."
            )

        extended = dict(structural_row)

        extended["SiV_BO3_NBO2_count"] = counts["BO3_NBO2"]
        extended["SiV_other_unclassified_count"] = (
            counts["UNCLASSIFIED"]
        )

        new_rows.append(extended)

        summary[x] = {
            **counts,
            "total": expected_siv,
        }

    return structural_fields, new_rows, summary


def build_validation_report(summary):
    """Construct a provenance and validation report."""

    lines = [
        "CaO-SiO2 CONSOLIDATED STRUCTURAL RESULTS V2",
        "=" * 55,
        "",
        "PURPOSE",
        "Versioned extension of the existing consolidated "
        "structural results.",
        "",
        "No molecular dynamics simulations were performed.",
        "No coordination analysis was recalculated.",
        "No ring enumeration was performed.",
        "Existing source files were read without modification.",
        "",
        "SiV CLASSIFICATION",
        "-" * 55,
        "BO4_NBO1: BO=4, NBO=1, other=0",
        "BO5_NBO0: BO=5, NBO=0, other=0",
        "BO3_NBO2: BO=3, NBO=2, other=0",
        "UNCLASSIFIED: any remaining five-coordinate pattern",
        "",
        "VALIDATED SiV POPULATIONS",
        "-" * 55,
        "x_CaO  BO4_NBO1  BO5_NBO0  BO3_NBO2  "
        "Unclassified  Total",
    ]

    totals = {
        "BO4_NBO1": 0,
        "BO5_NBO0": 0,
        "BO3_NBO2": 0,
        "UNCLASSIFIED": 0,
        "total": 0,
    }

    for x in EXPECTED_X:
        result = summary[x]

        lines.append(
            f"{x:.2f}    "
            f"{result['BO4_NBO1']:>8}  "
            f"{result['BO5_NBO0']:>8}  "
            f"{result['BO3_NBO2']:>8}  "
            f"{result['UNCLASSIFIED']:>12}  "
            f"{result['total']:>5}"
        )

        for key in totals:
            totals[key] += result[key]

    lines.extend([
        "",
        f"TOTAL   "
        f"{totals['BO4_NBO1']:>8}  "
        f"{totals['BO5_NBO0']:>8}  "
        f"{totals['BO3_NBO2']:>8}  "
        f"{totals['UNCLASSIFIED']:>12}  "
        f"{totals['total']:>5}",
        "",
        "VALIDATION",
        "-" * 55,
        "All four compositions were checked.",
        "Every SiV site has five classified oxygen neighbors.",
        "Original BO4_NBO1 counts were preserved.",
        "Original BO5_NBO0 counts were preserved.",
        "BO3_NBO2 sites were identified from site-level ligand counts.",
        "Unclassified sites were retained as a separate category.",
        "All subtype populations sum to the structural SiV_count.",
        "",
        "SOFTWARE ENVIRONMENT",
        "-" * 55,
        f"Python: {sys.version.split()[0]}",
        f"Platform: {platform.platform()}",
        "",
        "INPUT FILE SHA256",
        "-" * 55,
    ])

    hash_files = [
        STRUCTURAL_CSV,
        SIV_SITES_CSV,
        RING_CSV,
        *STRUCTURE_FILES,
        *PROVENANCE_SCRIPTS,
        Path(__file__).resolve(),
    ]

    for path in hash_files:
        if path.is_file():
            relative = path.relative_to(ROOT)
            lines.append(f"{relative}: {sha256_file(path)}")
        else:
            lines.append(
                f"{path}: NOT FOUND (not hashed)"
            )

    lines.extend([
        "",
        "END OF VALIDATION REPORT",
        "",
    ])

    return "\n".join(lines)


def write_new_file(path, content):
    """Create a new file exclusively; never overwrite."""
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(content)


def main():
    print("CaO-SiO2 structural results V2")
    print("=" * 45)

    # Check all required source files before processing.
    required_inputs = [
        STRUCTURAL_CSV,
        SIV_SITES_CSV,
        RING_CSV,
        *STRUCTURE_FILES,
    ]

    missing = [path for path in required_inputs if not path.is_file()]

    if missing:
        raise FileNotFoundError(
            "Required input files are missing:\n"
            + "\n".join(str(path) for path in missing)
        )

    # Refuse to overwrite either output.
    existing = [
        path for path in (OUTPUT_CSV, OUTPUT_VALIDATION)
        if path.exists()
    ]

    if existing:
        raise FileExistsError(
            "Output file(s) already exist. Nothing was written:\n"
            + "\n".join(str(path) for path in existing)
        )

    fields, rows, summary = build_results()

    # Prepare both complete outputs in memory before writing.
    output_fields = fields + [
        "SiV_BO3_NBO2_count",
        "SiV_other_unclassified_count",
    ]

    csv_buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        csv_buffer,
        fieldnames=output_fields,
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)

    validation_text = build_validation_report(summary)

    # Create new outputs exclusively.
    write_new_file(OUTPUT_CSV, csv_buffer.getvalue())
    write_new_file(OUTPUT_VALIDATION, validation_text)

    print("\nValidation passed.")
    print(f"Created: {OUTPUT_CSV.relative_to(ROOT)}")
    print(f"Created: {OUTPUT_VALIDATION.relative_to(ROOT)}")

    print("\nSiV subtype summary:")
    for x in EXPECTED_X:
        print(f"x={x:.2f}: {summary[x]}")

    print("\nOriginal source files were not modified.")


if __name__ == "__main__":
    main()
