from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.modifiers import ModifierAnalysis
from workflow.parameters import BOND_CUTOFFS


system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

neighbors = NeighborList(
    system,
    BOND_CUTOFFS
)

neighbors.build()

analysis = ModifierAnalysis(
    system,
    neighbors
)

stats = analysis.coordination_statistics()
totals = analysis.global_environment()

assert stats["number_ca"] == 720
assert abs(stats["mean_ca_o"] - 6.157) < 0.01

assert totals["BO"] == 1239
assert totals["NBO"] == 3160
assert totals["other"] == 34

print("===== Modifier Validation =====")
print(f"Ca atoms       : {stats['number_ca']}")
print(f"Mean Ca-O      : {stats['mean_ca_o']:.3f}")
print(f"BO contacts    : {totals['BO']}")
print(f"NBO contacts   : {totals['NBO']}")
print(f"Other contacts : {totals['other']}")
print("Modifier validation: PASS")
