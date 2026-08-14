"""
Master structural report for CaO-SiO2 glass analysis.

This module orchestrates the validated analysis modules without
duplicating their scientific logic.
"""

from workflow.density import DensityAnalysis
from workflow.coordination import CoordinationAnalysis
from workflow.partial_rdf import PartialRDF
from workflow.qn import QnAnalysis
from workflow.topology import TopologyAnalysis
from workflow.rings import RingAnalysis
from workflow.modifiers import ModifierAnalysis
from workflow.bond_angles import BondAngleAnalysis


class StructuralReport:
    """
    Collect all validated structural analyses for one configuration.
    """

    def __init__(
        self,
        system,
        neighbors,
        rdf_file=None,
    ):
        self.system = system
        self.neighbors = neighbors
        self.rdf_file = rdf_file

    # ---------------------------------------------------------
    # Complete machine-readable result
    # ---------------------------------------------------------

    def results(self):
        """
        Return all structural results as a nested dictionary.
        """

        density = DensityAnalysis(self.system)

        coordination = CoordinationAnalysis(
            self.system,
            self.neighbors,
        )

        qn = QnAnalysis(
            self.system,
            self.neighbors,
        )

        topology = TopologyAnalysis(
            self.system,
            self.neighbors,
        )

        rings = RingAnalysis(
            self.system,
            topology,
        )

        modifiers = ModifierAnalysis(
            self.system,
            self.neighbors,
        )

        angles = BondAngleAnalysis(
            self.system,
            self.neighbors,
        )

        si_indices = [
            i for i, atom_type in enumerate(self.system.types)
            if atom_type == 1
        ]

        si4_qn = {
            0: 0,
            1: 0,
            2: 0,
            3: 0,
            4: 0,
        }

        siv_count = 0
        siv_bo4_nbo1 = 0
        siv_bo5_nbo0 = 0

        for i in si_indices:

            oxygen_neighbors = [
                j
                for j in self.neighbors.neighbors_of(i)
                if self.system.types[j] == 2
            ]

            si_coordination = len(oxygen_neighbors)

            bo = sum(
                qn.oxygen_class(j) == "BO"
                for j in oxygen_neighbors
            )

            nbo = sum(
                qn.oxygen_class(j) == "NBO"
                for j in oxygen_neighbors
            )

            if si_coordination == 4:

                if bo in si4_qn:
                    si4_qn[bo] += 1

            elif si_coordination == 5:

                siv_count += 1

                if bo == 4 and nbo == 1:
                    siv_bo4_nbo1 += 1

                elif bo == 5 and nbo == 0:
                    siv_bo5_nbo0 += 1

        result = {
            "composition": density.composition(),
            "density": {
                "mass_amu": density.total_mass_amu(),
                "volume_angstrom3": density.volume_angstrom3(),
                "density_g_cm3": density.density_g_cm3(),
            },

            "coordination": {
                "Si-O": coordination.statistics(1, 2),
                "O-Si": coordination.statistics(2, 1),
                "Ca-O": coordination.statistics(3, 2),
                "O-Ca": coordination.statistics(2, 3),
            },

            "qn": {
                "oxygen": qn.oxygen_statistics(),
                "distribution": si4_qn,
                "SiV": {
                    "total": siv_count,
                    "BO4_NBO1": siv_bo4_nbo1,
                    "BO5_NBO0": siv_bo5_nbo0,
                },
            },
            "topology": {
                "network": topology.network_statistics(),
                "components": topology.component_statistics(),
            },

            "rings": {
                "statistics": rings.statistics(),
                "distribution": rings.distribution(),
            },

            "modifiers": {
                "environment": modifiers.global_environment(),
                "coordination": modifiers.coordination_statistics(),
            },

            "bond_angles": {
                "O-Si-O": angles.statistics(2, 1),
                "Si-O-Si": angles.statistics(1, 2),
                "O-Ca-O": angles.statistics(2, 3),
            },
        }

        # RDF is optional because it is an external LAMMPS output file.
        if self.rdf_file is not None:
            rdf = PartialRDF(self.rdf_file)
            result["rdf"] = rdf.report()

        return result

    # ---------------------------------------------------------
    # Human-readable report
    # ---------------------------------------------------------

    def print_report(self):
        """
        Print a compact complete structural report.
        """

        data = self.results()

        print()
        print("=" * 60)
        print("           CaO-SiO2 STRUCTURAL REPORT")
        print("=" * 60)

        print()
        print("===== Composition =====")
        print(data["composition"])

        print()
        print("===== Density =====")
        d = data["density"]
        print(f"Mass        : {d['mass_amu']:.4f} amu")
        print(f"Volume      : {d['volume_angstrom3']:.3f} Å^3")
        print(f"Density     : {d['density_g_cm3']:.4f} g/cm^3")

        print()
        print("===== Coordination =====")
        for pair, stats in data["coordination"].items():
            print(
                f"{pair:6s} : "
                f"mean={stats['mean']:.3f}, "
                f"std={stats['std']:.3f}, "
                f"range={stats['min']}–{stats['max']}"
            )

        print()
        print("===== Qn / Oxygen Network =====")
        oxygen = data["qn"]["oxygen"]
        print(
            f"BO={oxygen['BO']} "
            f"({100 * oxygen['BO_fraction']:.2f} %)  "
            f"NBO={oxygen['NBO']} "
            f"({100 * oxygen['NBO_fraction']:.2f} %)  "
            f"Other={oxygen['other']} "
            f"({100 * oxygen['other_fraction']:.2f} %)"
        )

        print(
            f"Qn distribution (SiIV only) : "
            f"{data['qn']['distribution']}"
        )

        qn_distribution = data["qn"]["distribution"]
        qn_total = sum(qn_distribution.values())

        qn_fractions = {
            q: qn_distribution[q] / qn_total
            for q in range(5)
        }

        print(
            "Qn fractions (SiIV only) : "
            f"{qn_fractions}"
        )

        print(
            "SiV total       : "
            f"{data['qn']['SiV']['total']}"
        )
        print(
            "SiV(BO4,NBO1)   : "
            f"{data['qn']['SiV']['BO4_NBO1']}"
        )
        print(
            "SiV(BO5,NBO0)   : "
            f"{data['qn']['SiV']['BO5_NBO0']}"
        )

        print()
        print("===== Silicate Network =====")
        print(data["topology"]["network"])
        print(data["topology"]["components"])

        print()
        print("===== Ring Topology =====")
        print(data["rings"]["statistics"])
        print(f"Distribution : {data['rings']['distribution']}")

        print()
        print("===== Ca Modifier Environment =====")
        print(data["modifiers"]["environment"])
        print(data["modifiers"]["coordination"])

        print()
        print("===== Bond Angles =====")
        for name, stats in data["bond_angles"].items():
            print(
                f"{name:8s} : "
                f"mean={stats['mean']:.3f} deg, "
                f"std={stats['std']:.3f} deg"
            )

        if "rdf" in data:
            print()
            print("===== Partial RDF =====")
            for pair, values in data["rdf"].items():
                peak = values["peak"]
                minimum = values["first_minimum"]

                print(
                    f"{pair:6s} : "
                    f"peak={peak['r']:.3f} Å "
                    f"(g={peak['g']:.3f}), "
                    f"minimum={minimum['r']:.3f} Å"
                )

        print()
        print("=" * 60)
        print("Structural report completed.")
        print("=" * 60)
