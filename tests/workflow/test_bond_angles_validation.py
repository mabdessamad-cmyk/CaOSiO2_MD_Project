from workflow.lammps_reader import read_lammps
from workflow.neighbors import NeighborList
from workflow.bond_angles import BondAngleAnalysis
from workflow.parameters import BOND_CUTOFFS


FILENAME = "simulations/data/glass_final_300K.data"


def test_bond_angles_validation():
    system = read_lammps(FILENAME)

    neighbors = NeighborList(system, BOND_CUTOFFS)
    neighbors.build()

    analysis = BondAngleAnalysis(system, neighbors)

    for name, outer, center in [
        ("O-Si-O", 2, 1),
        ("Si-O-Si", 1, 2),
        ("O-Ca-O", 2, 3),
    ]:
        values = analysis.angles_for_triplet(outer, center)

        assert len(values) > 0, f"{name}: no angles found"
        assert values.min() >= 0.0, f"{name}: angle below 0"
        assert values.max() <= 180.0, f"{name}: angle above 180"
        assert values.shape[0] == len(values)
