"""Generate bounded periodic-ring data, figures, and a methods note."""

import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.parameters import BOND_CUTOFFS
from workflow.periodic_rings import PeriodicRingAnalysis


ROOT = Path(__file__).resolve().parents[1]

STRUCTURES = [
    (0.10, "simulations/data/x010_glass_final_300K.data"),
    (0.20, "simulations/data/x020_glass_final_300K.data"),
    (0.30, "simulations/data/glass_final_300K.data"),
    (0.40, "simulations/data/x040_glass_final_300K.data"),
]

MAX_RING_SIZES = (6, 8, 10)
CLASSIFICATIONS = ("contractible", "winding", "invalid")

CSV_PATH = ROOT / "analysis/reports/periodic_ring_counts.csv"
FIGURE_PNG = ROOT / "analysis/figures/periodic_ring_analysis.png"
FIGURE_PDF = ROOT / "analysis/figures/periodic_ring_analysis.pdf"
METHODS_PATH = ROOT / "analysis/reports/periodic_ring_methods.txt"


def main():
    for directory in (
        CSV_PATH.parent,
        FIGURE_PNG.parent,
        FIGURE_PDF.parent,
        METHODS_PATH.parent,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    csv_rows = []
    summary = {}

    for x, relative_path in STRUCTURES:
        structure_path = ROOT / relative_path
        if not structure_path.is_file():
            raise FileNotFoundError(
                f"Required structure is missing: {structure_path}"
            )

        print(f"\n{'=' * 65}")
        print(f"CaO fraction x = {x:.2f}")
        print(f"Structure: {relative_path}")

        system = read_lammps(str(structure_path))
        n_si = int(np.count_nonzero(np.asarray(system.types) == 1))
        if n_si == 0:
            raise ValueError(f"No Si atoms found in {relative_path}")

        # Build the validated atomic neighbor list only once per structure.
        neighbors = NeighborList(system, BOND_CUTOFFS)
        neighbors.build()

        for max_size in MAX_RING_SIZES:
            print(f"Enumerating rings up to size {max_size}...")
            analysis = PeriodicRingAnalysis(
                system,
                neighbors,
                max_ring_size=max_size,
            )
            records = analysis.rings()

            size_counts = Counter(int(record["size"]) for record in records)
            class_counts = Counter(
                record["classification"] for record in records
            )
            total = len(records)

            summary[(x, max_size)] = {
                "n_si": n_si,
                "total": total,
                "size_counts": size_counts,
                "class_counts": class_counts,
            }

            for ring_size in range(3, max_size + 1):
                count = size_counts.get(ring_size, 0)
                csv_rows.append({
                    "x_CaO": f"{x:.2f}",
                    "n_Si": n_si,
                    "max_ring_size": max_size,
                    "ring_size_Si": ring_size,
                    "ring_count": count,
                    "rings_per_Si": f"{count / n_si:.8f}",
                    "fraction_of_rings": (
                        f"{count / total:.8f}" if total else "0.00000000"
                    ),
                    "total_rings": total,
                    "contractible_count": class_counts.get(
                        "contractible", 0
                    ),
                    "winding_count": class_counts.get("winding", 0),
                    "invalid_count": class_counts.get("invalid", 0),
                })

            print(
                f"  total={total}; "
                f"size distribution={dict(sorted(size_counts.items()))}; "
                f"contractible={class_counts.get('contractible', 0)}; "
                f"winding={class_counts.get('winding', 0)}; "
                f"invalid={class_counts.get('invalid', 0)}"
            )

    fieldnames = [
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
    ]

    with CSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

    # Figure: ring-size distribution at limit 10, plus limit sensitivity.
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12, 4.8), constrained_layout=True
    )

    max_size = 10
    ring_sizes = np.arange(3, max_size + 1)

    for x, _ in STRUCTURES:
        result = summary[(x, max_size)]
        values = [
            result["size_counts"].get(int(size), 0) / result["n_si"]
            for size in ring_sizes
        ]
        ax1.plot(
            ring_sizes,
            values,
            marker="o",
            markersize=4,
            linewidth=1.5,
            label=f"x = {x:.2f}",
        )

    ax1.set_xlabel("Ring size (number of Si atoms)")
    ax1.set_ylabel("Ring count per Si atom")
    ax1.set_title("(a) Ring-size distribution (limit = 10)")
    ax1.set_xticks(ring_sizes)
    ax1.grid(True, alpha=0.25)
    ax1.legend(frameon=False, ncol=2)

    for max_size in MAX_RING_SIZES:
        xs = [x for x, _ in STRUCTURES]
        values = [
            summary[(x, max_size)]["total"]
            / summary[(x, max_size)]["n_si"]
            for x in xs
        ]
        ax2.plot(
            xs,
            values,
            marker="o",
            markersize=4,
            linewidth=1.5,
            label=f"Maximum size = {max_size}",
        )

    ax2.set_xlabel("CaO mole fraction, x")
    ax2.set_ylabel("Total bounded ring count per Si atom")
    ax2.set_title("(b) Sensitivity to maximum ring size")
    ax2.set_xticks([x for x, _ in STRUCTURES])
    ax2.grid(True, alpha=0.25)
    ax2.legend(frameon=False)

    fig.savefig(FIGURE_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(FIGURE_PDF, bbox_inches="tight")
    plt.close(fig)

    methods = """Periodic Si-network ring analysis
=================================

Input structures
----------------
x = 0.10: simulations/data/x010_glass_final_300K.data
x = 0.20: simulations/data/x020_glass_final_300K.data
x = 0.30: simulations/data/glass_final_300K.data
x = 0.40: simulations/data/x040_glass_final_300K.data

Structural definition
---------------------
The Si network is constructed from bridging oxygens (BOs), each
connecting exactly two Si atoms under the existing pair-specific
neighbor cutoffs. The adopted Si-O and Ca-O cutoffs remain 2.30 A
and 3.10 A, respectively.

Candidate rings are chordless, primitive cycles in the underlying
simple Si graph, with 3 through the selected maximum number of Si
vertices. Each distinct BO identity on a cycle edge is retained;
therefore, distinct oxygen-labelled paths can generate distinct ring
records even when they share the same Si cycle.

Periodic classification
-----------------------
For each enumerated ring, the local minimum-image bond displacement
vectors are summed around the closed path. The resulting translation
is classified as contractible, winding, or invalid using the analyzer's
periodic-winding tolerance.

Reported limits and normalization
---------------------------------
Maximum ring sizes: 6, 8, and 10 Si atoms.
The CSV reports counts by ring size, ring count divided by the number
of Si atoms, fraction of all enumerated rings at that limit, and the
counts in each winding classification.

Interpretation and limitations
------------------------------
These are bounded, method-dependent primitive ring counts, not a
complete census of every possible ring in the periodic multigraph.
The analysis does not enumerate two-edge, two-Si cycles created by
parallel BO edges. Rings overlap, so rings per Si is a normalized
count, not the fraction of Si atoms belonging to rings.

A zero winding count means no winding ring was identified among the
rings enumerated at that size limit; it does not prove that no larger
winding rings exist. Comparisons across compositions should therefore
consider both the common size limits and the dependence on that limit.
"""

    METHODS_PATH.write_text(methods, encoding="utf-8")

    print(f"\n{'=' * 65}")
    print("Analysis completed.")
    print(f"CSV:         {CSV_PATH.relative_to(ROOT)}")
    print(f"Figure PNG:  {FIGURE_PNG.relative_to(ROOT)}")
    print(f"Figure PDF:  {FIGURE_PDF.relative_to(ROOT)}")
    print(f"Methods:     {METHODS_PATH.relative_to(ROOT)}")

    print("\nSummary: total bounded rings and rings per Si")
    print("x_CaO  max_size  total_rings  rings_per_Si")
    for x, _ in STRUCTURES:
        for max_size in MAX_RING_SIZES:
            result = summary[(x, max_size)]
            print(
                f"{x:5.2f}  {max_size:8d}  "
                f"{result['total']:12d}  "
                f"{result['total'] / result['n_si']:.8f}"
            )


if __name__ == "__main__":
    main()
