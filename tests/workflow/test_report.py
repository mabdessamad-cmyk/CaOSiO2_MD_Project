"""
Integration test for the master structural report.
"""

from workflow.lammps_reader import read_lammps
from workflow.parameters import BOND_CUTOFFS
from workflow.neighbors import NeighborList
from workflow.report import StructuralReport


DATA_FILE = "simulations/data/glass_final_300K.data"


print("===== Structural Report Validation =====")

# ---------------------------------------------------------
# Read structure
# ---------------------------------------------------------

system = read_lammps(DATA_FILE)

# ---------------------------------------------------------
# Build validated neighbor list
# ---------------------------------------------------------

neighbors = NeighborList(system, BOND_CUTOFFS)
neighbors.build()

# ---------------------------------------------------------
# Build master report
# ---------------------------------------------------------

report = StructuralReport(
    system,
    neighbors,
    rdf_file="analysis/tables/partial_rdf.dat",
)

results = report.results()

# ---------------------------------------------------------
# Validate required sections
# ---------------------------------------------------------

required = {
    "composition",
    "density",
    "coordination",
    "qn",
    "topology",
    "rings",
    "modifiers",
    "bond_angles",
    "rdf",
}

missing = required.difference(results.keys())

assert not missing, (
    f"Missing report sections: {sorted(missing)}"
)

# ---------------------------------------------------------
# Basic structural checks
# ---------------------------------------------------------

assert results["composition"]["Si"] == 1680
assert results["composition"]["O"] == 4080
assert results["composition"]["Ca"] == 720

assert results["density"]["density_g_cm3"] > 0.0

assert results["qn"]["oxygen"]["total"] == 4080

assert results["topology"]["network"]["number_si"] == 1680

assert results["modifiers"]["coordination"]["number_ca"] == 720

assert "Si-O" in results["rdf"]
assert "O-Ca" in results["rdf"]

assert sum(results["qn"]["distribution"].values()) == 1675

assert results["qn"]["distribution"] == {
    0: 4,
    1: 37,
    2: 280,
    3: 729,
    4: 625,
}

assert results["qn"]["SiV"]["total"] == 5
assert results["qn"]["SiV"]["BO4_NBO1"] == 3
assert results["qn"]["SiV"]["BO5_NBO0"] == 2
print("Report sections : PASS")
print("Composition      : PASS")
print("Density           : PASS")
print("Qn                : PASS")
print("Topology          : PASS")
print("Modifiers         : PASS")
print("Partial RDF       : PASS")
print()
print("Structural report validation: PASS")
