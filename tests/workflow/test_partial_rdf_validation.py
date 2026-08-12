from workflow.partial_rdf import PartialRDF


RDF_FILE = "analysis/tables/partial_rdf.dat"


rdf = PartialRDF(RDF_FILE).read()

assert rdf.data.shape == (600, 14)
assert len(rdf.r) == 600

expected_pairs = {
    "Si-Si",
    "Si-O",
    "Si-Ca",
    "O-O",
    "O-Ca",
    "Ca-Ca",
}

report = rdf.report()

assert set(report) == expected_pairs

# Previously validated structural cutoffs
assert abs(report["Si-O"]["first_minimum"]["r"] - 2.325) < 0.051
assert abs(report["O-Ca"]["first_minimum"]["r"] - 3.675) < 0.051

# Basic physical consistency
for pair, result in report.items():
    peak_r = result["peak"]["r"]
    minimum_r = result["first_minimum"]["r"]

    assert peak_r > 0.0
    assert minimum_r > peak_r
    assert result["peak"]["g"] > 0.0

print("===== Partial RDF Validation =====")
print(f"RDF rows       : {len(rdf.r)}")
print(f"RDF columns    : {rdf.data.shape[1]}")
print(f"Pair channels  : {len(report)}")
print(
    f"Si-O minimum   : "
    f"{report['Si-O']['first_minimum']['r']:.3f} Å"
)
print(
    f"O-Ca minimum   : "
    f"{report['O-Ca']['first_minimum']['r']:.3f} Å"
)
print("Partial RDF validation: PASS")
