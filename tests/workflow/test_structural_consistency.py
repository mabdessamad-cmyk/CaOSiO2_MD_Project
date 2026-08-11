from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.qn import QnAnalysis
from workflow.topology import TopologyAnalysis
from workflow.parameters import BOND_CUTOFFS


system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

neighbors = NeighborList(system, BOND_CUTOFFS)
neighbors.build()

qn = QnAnalysis(system, neighbors)
topology = TopologyAnalysis(system, neighbors)

oxygen = qn.oxygen_statistics()
qdist = qn.qn_distribution()
network = topology.network_statistics()

print("===== Structural Consistency =====")

# ---------------------------------------------------------
# Composition
# ---------------------------------------------------------

n_si = sum(system.types == 1)
n_o = sum(system.types == 2)
n_ca = sum(system.types == 3)

assert n_si == 1680
assert n_o == 4080
assert n_ca == 720

# ---------------------------------------------------------
# BO / NBO / Other oxygen consistency
# ---------------------------------------------------------

bo = oxygen["BO"]
nbo = oxygen["NBO"]
other = oxygen["other"]

assert bo + nbo + other == n_o

# ---------------------------------------------------------
# BO == Si-Si network edges
# ---------------------------------------------------------

edges = network["bridging_oxygen_count"]

assert bo == edges

# ---------------------------------------------------------
# Qn population
# ---------------------------------------------------------

qn_total = sum(
    qdist[q]
    for q in range(5)
)

invalid = qdist["invalid"]

assert qn_total + invalid == n_si

# ---------------------------------------------------------
# NBO predicted from Qn
# ---------------------------------------------------------

nbo_from_qn = sum(
    (4 - q) * qdist[q]
    for q in range(5)
)

print(
    f"NBO from ideal Qn relation : {nbo_from_qn}"
)

print(
    f"Actual NBO                 : {nbo}"
)

print(
    f"Difference                 : "
    f"{nbo - nbo_from_qn}"
)

assert nbo_from_qn <= nbo

print("Structural consistency: PASS")
print(f"Si atoms              : {n_si}")
print(f"O atoms               : {n_o}")
print(f"Ca atoms              : {n_ca}")
print(f"BO                    : {bo}")
print(f"NBO                   : {nbo}")
print(f"Other O               : {other}")
print(f"Si-Si edges           : {edges}")
print(f"NBO from Qn           : {nbo_from_qn}")
print("Structural consistency: PASS")
