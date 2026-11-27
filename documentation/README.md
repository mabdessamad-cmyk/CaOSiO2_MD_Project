# CaO–SiO₂ Glass Molecular Dynamics and Structural Analysis
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23260918.svg)](https://doi.org/10.5281/zenodo.23260918)

## Project

Molecular dynamics simulation and structural analysis of binary calcium silicate glasses of composition

**x(CaO)–(1−x)SiO₂**

with

**x = 0.10, 0.20, 0.30, and 0.40.**

The simulations were performed with **LAMMPS** using the **Pedone et al. (2006)** interaction potential. The project provides the simulation inputs, final and intermediate atomic configurations, structural-analysis scripts, partial radial distribution functions, figures, tables, and consolidated structural reports associated with the CaO–SiO₂ glass study.

## Scientific scope

The study examines the composition-driven structural evolution of CaO–SiO₂ glasses through:

* atomic density and structural relaxation;
* Si–O and Ca–O coordination;
* Si oxygen-coordination states;
* Qⁿ speciation of SiIV sites;
* bridging and non-bridging oxygen populations;
* Si–O–Si network connectivity;
* largest-network fraction and connected-component analysis;
* ring statistics;
* O–Si–O, Si–O–Si, and O–Ca–O bond angles;
* partial radial distribution functions.

The structural analysis is based on chemically specific distance cutoffs and a common analysis protocol applied across all four compositions.

## Simulation method

The production simulations use:

* **LAMMPS 29 Aug 2024**
* **Pedone et al. (2006)** potential
* periodic boundary conditions
* `pedone` short-range interactions combined with `coul/long`
* PPPM treatment of long-range Coulomb interactions
* classical molecular dynamics with a 0.5 fs timestep during the production melt–quench protocol
* melt–quench simulations followed by final 300 K NPT relaxation.

The final structures contain approximately 6,500–6,800 atoms per composition.

The exact simulation parameters are contained in the LAMMPS input files included in the repository.

## Compositions and final structures

| CaO fraction x | Composition  | Final atom count |
| -------------: | ------------ | ---------------: |
|           0.10 | 10CaO–90SiO₂ |             6496 |
|           0.20 | 20CaO–80SiO₂ |             6720 |
|           0.30 | 30CaO–70SiO₂ |             6480 |
|           0.40 | 40CaO–60SiO₂ |             6760 |

## Repository structure

```text
CaOSiO2_MD_Project/

├── lammps/
│   ├── data/
│   ├── forcefields/
│   ├── inputs/
│   └── logs/
│
├── simulations/
│   ├── data/
│   ├── inputs/
│   └── logs/
│
├── analysis/
│   ├── scripts/
│   ├── figures/
│   └── tables/
│
├── workflow/
│   └── analysis and workflow scripts
│
├── generator/
│   └── structure-generation and validation scripts
│
├── tests/
│   └── validation tests
│
├── output/
│   └── generated and consolidated project outputs
│
└── documentation/
    └── project documentation
```

Large trajectory and restart files are intentionally excluded from the public repository. The retained atomic configurations and analysis outputs are sufficient to document and reproduce the reported structural analysis.

## Structural-analysis workflow

The principal analysis sequence is:

1. read the final equilibrated glass configurations;
2. determine chemically specific coordination numbers;
3. classify Si oxygen coordination;
4. calculate Qⁿ distributions for SiIV sites;
5. classify oxygen atoms as bridging, non-bridging, or other according to Si connectivity;
6. construct the Si–O network through bridging oxygen;
7. determine connected components and the largest-network fraction;
8. calculate ring statistics;
9. calculate O–Si–O, Si–O–Si, and O–Ca–O angle distributions;
10. calculate partial radial distribution functions for the relevant network and modifier pairs;
11. generate consolidated structural reports, tables, and figures.

The analysis scripts used to generate the reported structural quantities are included in `analysis/` and `workflow/`.

## Reproducibility

The repository preserves the validated scientific state of the project together with the simulation inputs, final configurations, structural-analysis scripts, and derived structural results.

The scientific/data state used for the manuscript is identified by Git commit:

```text
a5c065e
```

The subsequent commits contain repository-cleanup operations only; they do not change the scientific workflow or reported structural results.

## Reference

The interaction potential used in the simulations is:

Pedone, A.; Malavasi, G.; Menziani, M. C.; Cormack, A. N.; Segre, U.
*Atomic-scale modeling of structural, vibrational and electronic properties of (CaO)x(SiO2)1−x glasses.*
**J. Phys. Chem. B 110 (2006) 11780–11795.**
DOI: 10.1021/jp0611018.

The manuscript associated with this repository provides the complete methodological and structural-analysis description.

## Data availability

The simulation inputs, atomic configurations, structural-analysis scripts, partial radial distribution functions, figures, tables, and consolidated structural reports are provided in this repository.

The repository is archived on Zenodo. The [concept DOI](https://doi.org/10.5281/zenodo.23260918) resolves to the record for all versions, while the [v1.1.0 version DOI](https://doi.org/10.5281/zenodo.24478109) identifies the specific archived release published on 23 November 2026. The version-specific DOI should be cited when referring to the exact release used.

## License

The source code and software components in this repository are licensed under the MIT License. See `LICENSE`.

Research data, figures, tables, and derived research materials are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), unless otherwise stated. See `DATA_LICENSE.md`.

