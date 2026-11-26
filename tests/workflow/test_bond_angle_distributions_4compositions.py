from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy as np


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "analysis/scripts/structural/bond_angle_distributions_4compositions.py"
)

spec = spec_from_file_location("bond_angle_distributions_4compositions", SCRIPT)
module = module_from_spec(spec)
spec.loader.exec_module(module)


class FakeNeighbors:
    def __init__(self, mapping):
        self.mapping = mapping

    def neighbors_of(self, atom):
        return self.mapping.get(atom, [])


class FakeAngleAnalysis:
    def __init__(self):
        self.calls = []

    def angle(self, atom1, center, atom2):
        self.calls.append((atom1, center, atom2))
        return 123.45


def test_bo_only_siosi_requires_exactly_two_si_neighbors():
    # Atom types: Si, O, O, O, Si, Si
    # Oxygen 1 has one Si neighbour.
    # Oxygen 2 has exactly two Si neighbours.
    # Oxygen 3 has three Si neighbours.
    system = SimpleNamespace(
        types=np.array([1, 2, 2, 2, 1, 1])
    )

    neighbors = FakeNeighbors({
        1: [0],
        2: [0, 4],
        3: [0, 4, 5],
    })

    angle_analysis = FakeAngleAnalysis()

    angles, bo_count = module.calculate_bo_only_siosi(
        system,
        neighbors,
        angle_analysis,
    )

    assert bo_count == 1
    assert len(angles) == 1
    assert np.allclose(angles, [123.45])
    assert angle_analysis.calls == [(0, 2, 4)]
