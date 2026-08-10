"""
Network-topology analysis for CaO-SiO2 glasses.

The topology is defined through Si-O-Si bridging oxygens.
"""

import numpy as np


SI = 1
O = 2


class TopologyAnalysis:
    """
    Analyze the silicate network formed by Si atoms connected
    through bridging oxygen atoms.
    """

    def __init__(self, system, neighbors):
        self.system = system
        self.neighbors = neighbors

    def si_neighbors_of_oxygen(self, oxygen_index):
        """Return Si neighbors of one oxygen atom."""
        return [
            i
            for i in self.neighbors.neighbors_of(oxygen_index)
            if self.system.types[i] == SI
        ]

    def bridging_oxygens(self):
        """
        Return oxygen indices connected to exactly two Si atoms.
        """

        result = []

        for oxygen in np.where(self.system.types == O)[0]:

            if len(self.si_neighbors_of_oxygen(oxygen)) == 2:
                result.append(oxygen)

        return np.asarray(result, dtype=int)

    def si_connections(self):
        """
        Build Si-Si connections through bridging oxygens.

        Returns
        -------
        dict
            Si index -> set of connected Si indices.
        """

        si_indices = np.where(
            self.system.types == SI
        )[0]

        connections = {
            int(i): set()
            for i in si_indices
        }

        for oxygen in self.bridging_oxygens():

            si = self.si_neighbors_of_oxygen(oxygen)

            i, j = si

            connections[i].add(j)
            connections[j].add(i)

        return connections

    def connected_components(self):
        """
        Identify connected Si-network components.

        Two Si atoms belong to the same component if they are
        connected through one or more bridging oxygens.

        Returns
        -------
        list of sets
            Each set contains the Si indices belonging to one
            connected component.
        """

        connections = self.si_connections()

        visited = set()
        components = []

        for start in connections:

            if start in visited:
                continue

            component = set()
            stack = [start]

            while stack:

                current = stack.pop()

                if current in visited:
                    continue

                visited.add(current)
                component.add(current)

                for neighbor in connections[current]:

                    if neighbor not in visited:
                        stack.append(neighbor)

            components.append(component)

        return components

    def component_statistics(self):
        """
        Calculate connected-component statistics.
        """

        components = self.connected_components()

        sizes = np.asarray(
            [len(component) for component in components],
            dtype=int
        )

        largest = int(np.max(sizes))

        percolating_fraction = (
            largest / len(self.system.types[self.system.types == 1])
        )

        return {
            "number_components": len(sizes),
            "largest_component": largest,
            "smallest_component": int(np.min(sizes)),
            "mean_component": float(np.mean(sizes)),
            "isolated_si": int(np.sum(sizes == 1)),
            "percolating_fraction": float(percolating_fraction),
        }


    def network_statistics(self):
        """
        Calculate basic Si-network statistics.
        """

        connections = self.si_connections()

        degrees = np.asarray(
            [len(v) for v in connections.values()],
            dtype=int
        )

        return {
            "number_si": len(degrees),
            "mean_degree": float(np.mean(degrees)),
            "std_degree": float(np.std(degrees)),
            "minimum_degree": int(np.min(degrees)),
            "maximum_degree": int(np.max(degrees)),
            "bridging_oxygen_count": len(
                self.bridging_oxygens()
            ),
        }

    def report(self):
        """Print network-topology statistics."""

        stats = self.network_statistics()

        print()
        print("===== Silicate Network Topology =====")
        print()
        print(
            f"Si atoms              : "
            f"{stats['number_si']}"
        )
        print(
            f"Bridging oxygens      : "
            f"{stats['bridging_oxygen_count']}"
        )
        print(
            f"Mean Si network degree: "
            f"{stats['mean_degree']:.3f}"
        )
        print(
            f"Std deviation         : "
            f"{stats['std_degree']:.3f}"
        )
        print(
            f"Minimum degree        : "
            f"{stats['minimum_degree']}"
        )
        print(
            f"Maximum degree        : "
            f"{stats['maximum_degree']}"
        )

        components = self.component_statistics()

        print()
        print("Connected components")
        print(
            f"Number of components : "
            f"{components['number_components']}"
        )
        print(
            f"Largest component     : "
            f"{components['largest_component']} Si"
        )
        print(
            f"Smallest component    : "
            f"{components['smallest_component']} Si"
        )
        print(
            f"Mean component size   : "
            f"{components['mean_component']:.3f} Si"
        )
        print(
            f"Isolated Si atoms     : "
            f"{components['isolated_si']}"
        )
        print(
            f"Percolating Si fraction : "
            f"{100.0 * components['percolating_fraction']:.2f} %"
        )
