# CaOSiO2 Glass Molecular Dynamics Workflow

## Project

Molecular Dynamics simulation and structural analysis of
xCaO–(1−x)SiO₂ glasses using LAMMPS.

Current force field:
- Pedone et al. (2006)

Current composition:
- 30CaO–70SiO₂

Future compositions:
- 10CaO–90SiO₂
- 20CaO–80SiO₂
- 30CaO–70SiO₂
- 40CaO–60SiO₂

---

## Directory structure

CaOSiO2_generator/

generator/
    Glass generation scripts

simulations/
    inputs/
    data/
    restart/
    trajectories/
    logs/

analysis/
    scripts/
    figures/
    tables/

tests/
    Validation scripts

workflow/
    Modular version under development

documentation/
    Notes and documentation

---

## Simulation workflow

1. Generate initial structure
2. Energy minimization
3. Relaxation at 300 K
4. Heating
5. Melting
6. Quenching
7. Final relaxation
8. RDF
9. Coordination analysis
10. Bond-angle analysis
11. BO/NBO analysis
12. Qn analysis

---

## Current status

✓ Stable glass generation

✓ Stable melt-quench procedure

✓ Density validated

✓ RDF validated

✓ Coordination validated

✓ Bond angles validated

✓ Qn validated

---

## Next objectives

- Complete workflow modularization
- Analyze all four compositions
- Publish structural paper
- Compute dielectric properties
- Publish dielectric/conductivity paper
