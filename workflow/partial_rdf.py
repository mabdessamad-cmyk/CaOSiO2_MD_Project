"""
Partial RDF analysis for CaO-SiO2 glass.

LAMMPS ordering:

    Si-Si, Si-O, Si-Ca, O-O, O-Ca, Ca-Ca

Each pair contributes:
    g(r), coordination(r)

RDF files may contain multiple time-averaged blocks. The blocks are
averaged at identical radial positions before peak/minimum analysis.
"""

from pathlib import Path
import numpy as np


PAIR_CHANNELS = {
    "Si-Si": (2, 3),
    "Si-O":  (4, 5),
    "Si-Ca": (6, 7),
    "O-O":   (8, 9),
    "O-Ca":  (10, 11),
    "Ca-Ca": (12, 13),
}


class PartialRDF:
    """Analyze chemically resolved partial RDFs."""

    def __init__(self, filename):
        self.filename = Path(filename)

        self.r = None
        self.data = None

        self.blocks = None
        self.n_blocks = 0
        self.n_bins = 0

    def read(self):
        """Read and average all RDF blocks in the file."""

        blocks = []
        current_block = []

        with self.filename.open() as fh:
            for line in fh:
                fields = line.split()

                # RDF data rows contain 14 columns.
                if len(fields) == 14:
                    current_block.append([float(x) for x in fields])
                    continue

                # A new block header means the previous block is complete.
                if len(fields) == 2:
                    if current_block:
                        blocks.append(current_block)
                        current_block = []

            # Store final block.
            if current_block:
                blocks.append(current_block)

        if not blocks:
            raise ValueError("No RDF blocks found.")

        # Validate that all blocks have identical dimensions.
        n_bins = len(blocks[0])

        if n_bins == 0:
            raise ValueError("First RDF block is empty.")

        if any(len(block) != n_bins for block in blocks):
            raise ValueError(
                "RDF blocks have different numbers of rows."
            )

        data_blocks = np.asarray(blocks, dtype=float)

        if data_blocks.ndim != 3 or data_blocks.shape[2] != 14:
            raise ValueError(
                f"Expected RDF data with shape "
                f"(n_blocks, n_bins, 14), got {data_blocks.shape}."
            )

        # All blocks must use the same radial grid.
        reference_r = data_blocks[0, :, 1]

        for i in range(1, data_blocks.shape[0]):
            if not np.allclose(
                data_blocks[i, :, 1],
                reference_r,
                rtol=0.0,
                atol=1.0e-12,
            ):
                raise ValueError(
                    "RDF blocks use different radial grids."
                )

        # Equal-weight average of the RDF blocks.
        data = np.mean(data_blocks, axis=0)

        self.blocks = data_blocks
        self.n_blocks = data_blocks.shape[0]
        self.n_bins = data_blocks.shape[1]

        self.data = data
        self.r = data[:, 1]

        return self

    def pair_rdf(self, pair):
        """Return r and averaged g(r) for one chemical pair."""

        if pair not in PAIR_CHANNELS:
            raise KeyError(f"Unknown pair: {pair}")

        g_col, _ = PAIR_CHANNELS[pair]

        return self.r, self.data[:, g_col]

    def pair_coordination(self, pair):
        """Return r and averaged coordination(r) for one chemical pair."""

        if pair not in PAIR_CHANNELS:
            raise KeyError(f"Unknown pair: {pair}")

        _, cn_col = PAIR_CHANNELS[pair]

        return self.r, self.data[:, cn_col]

    def peak(self, pair):
        """Return principal RDF peak position and height."""

        r, g = self.pair_rdf(pair)

        i = int(np.argmax(g))

        return {
            "r": float(r[i]),
            "g": float(g[i]),
        }

    def first_minimum(self, pair):
        """
        Return the first local minimum after the principal RDF peak.
        """

        r, g = self.pair_rdf(pair)

        peak_index = int(np.argmax(g))

        for i in range(peak_index + 1, len(g) - 1):
            if g[i] <= g[i - 1] and g[i] <= g[i + 1]:
                return {
                    "r": float(r[i]),
                    "g": float(g[i]),
                }

        raise ValueError(f"No local minimum after peak for {pair}.")

    def report(self):
        """Return a compact structural RDF report."""

        if self.data is None:
            self.read()

        report = {
            "_metadata": {
                "n_blocks": self.n_blocks,
                "n_bins": self.n_bins,
            }
        }

        for pair in PAIR_CHANNELS:
            report[pair] = {
                "peak": self.peak(pair),
                "first_minimum": self.first_minimum(pair),
            }

        return report
