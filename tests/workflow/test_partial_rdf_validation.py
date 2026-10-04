from workflow.partial_rdf import PartialRDF


RDF_FILE = "analysis/tables/partial_rdf.dat"


rdf = PartialRDF(RDF_FILE).read()

# Current RDF representation:
# two 10-ps blocks, each containing 300 radial bins and 14 columns.
assert rdf.data.shape == (300, 14)
assert rdf.blocks.shape == (2, 300, 14)
assert rdf.n_blocks == 2
assert rdf.n_bins == 300
assert len(rdf.r) == 300

expected_pairs = {
    "Si-Si",
    "Si-O",
    "Si-Ca",
    "O-O",
    "O-Ca",
    "Ca-Ca",
}

report = rdf.report()

# _metadata is an auxiliary report entry, not an RDF pair.
assert set(report) - {"_metadata"} == expected_pairs

# Previously validated structural RDF minima for the current RDF dataset.
assert abs(report["Si-O"]["first_minimum"]["r"] - 2.425) < 0.051
assert abs(report["O-Ca"]["first_minimum"]["r"] - 3.725) < 0.051

# Basic physical consistency
for pair in expected_pairs:
    result = report[pair]

    peak_r = result["peak"]["r"]
    minimum_r = result["first_minimum"]["r"]

    assert peak_r > 0.0
    assert minimum_r > peak_r
    assert result["peak"]["g"] > 0.0

print("===== Partial RDF Validation =====")
print(f"RDF rows       : {len(rdf.r)}")
print(f"RDF columns    : {rdf.data.shape[1]}")
print(f"RDF blocks     : {rdf.n_blocks}")
print(f"Pair channels  : {len(expected_pairs)}")
print(
    f"Si-O minimum   : "
    f"{report['Si-O']['first_minimum']['r']:.3f} Å"
)
print(
    f"O-Ca minimum   : "
    f"{report['O-Ca']['first_minimum']['r']:.3f} Å"
)
print("Partial RDF validation: PASS")
