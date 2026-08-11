from workflow.lammps_reader import read_lammps
from workflow.density import DensityAnalysis


system = read_lammps(
    "simulations/data/glass_final_300K.data"
)

density = DensityAnalysis(system)

density.report()
