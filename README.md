# Trachoma Transmission Simulation

An integrated epidemiological and genomic simulation framework for *Chlamydia trachomatis* (trachoma), combining agent-based SIS epidemic dynamics with forward genetic simulation (SLiM + msprime) to evaluate transmission inference methods.

This is the companion code repository for the manuscript:

> **"A simulation-based evaluation of Chlamydia trachomatis transmission directionality."**
> Mouhamadou Fadel Diop et al. *(in preparation)*

---

## Overview

This framework simulates trachoma outbreaks in structured populations (households + community) across three endemicity levels (hypo-, meso-, hyperendemic). It then runs pathogen genomes through SLiM and msprime, computes pairwise genetic metrics (IBD, IBS, phylogenetic distance), and benchmarks their ability to recover the known transmission network.

**Key findings:**
- IBD-based methods fail in trachoma (AUPRC ≈ 0.12) due to small bottleneck sizes, rapid recombination, and superinfection
- IBS and phylogenetic approaches retain moderate utility within a **50-day critical window**
- Bottleneck size ≤ 3 strains is the critical threshold for any genetic inference to work

---

## Repository Structure

```
.
├── R/                    # Core reusable R functions
│   ├── trachoma_epidemic_core.R     # SIS epidemic engine
│   ├── trachoma_calibration.R       # Parameter calibration to target incidence
│   ├── utils.R                      # General utilities
│   ├── utils_environment.R          # Python/SLiM environment setup
│   ├── fn_diagnostics.R             # Diagnostic and summary functions
│   ├── load_correct_python.R        # reticulate environment loader
│   └── pretty2D.R                   # 2D visualisation helpers
│
├── scripts/              # Pipeline entry-points (run from project root)
│   ├── trachoma_pipeline_script.R       # Main end-to-end pipeline (latest)
│   ├── trachoma_pipeline_script_v1.0.R  # Version 1.0 pipeline
│   ├── runme_trachoma_flat_two_sis_v1.R # Full module: epidemic + genetics
│   ├── runme_trachoma_flat_two_sis.R    # Variant of above
│   ├── runme_trachoma_flat_three_sis.R  # Three-strain SIS variant
│   ├── runme_trachoma_flat_one_sis.R    # Single-strain SIS variant
│   ├── runme_trachoma_abm.R             # Agent-based model only
│   ├── run_trachoma_full_flat.R         # Full flat-structure run
│   ├── recap_and_mutate_patch.R         # Recapitation/mutation patch
│   └── GET_SIMULATION_STRUCTURE.R       # Inspect simulation output structure
│
├── slim/                 # SLiM simulation scripts and wrappers
│   ├── slim_probe_genomes.slim          # Core SLiM genome simulation
│   ├── slim_probe_genomes_v2.slim       # Updated SLiM script
│   ├── run_slim_haplo_dropin_v6.R       # SLiM haplotype wrapper (pipeline-required)
│   ├── run_slim_haplo_dropin_v10.R      # Latest SLiM haplotype wrapper
│   ├── run_slim_dropin_v11.R            # Latest SLiM run wrapper
│   ├── finish_with_msprime_v6.R         # msprime recapitation (pipeline-required)
│   ├── finish_with_msprime_v9.R         # Latest msprime recapitation
│   ├── msprime_finish_dropin.R          # msprime drop-in helper
│   └── run_slim_probe.R / _v2.R         # SLiM probe utilities
│
├── analysis/             # Genetic metric computation and evaluation
│   ├── pairwise_sanity_v6.R             # Pairwise IBD/IBS/tree comparisons (pipeline-required)
│   ├── roc_transmission_eval.R          # ROC/PR curve evaluation
│   ├── roc_transmission_eval_time_v2.R  # Time-stratified evaluation (pipeline-required)
│   └── trac.R                           # Miscellaneous analysis
│
├── figures/              # Figure generation scripts (one script per figure)
│   ├── figure1_framework.R
│   ├── figure2_epidemic_model.R
│   ├── figure3_methods_performance.R
│   ├── figure4_temporal_dynamics.R
│   ├── figure5_household_confounding.R
│   ├── figure6_bottleneck_sensitivity.R
│   ├── figure3_ibd_erosion.R            # IBD erosion analysis for Figure 3
│   ├── figure3_concrete_example.R       # Step-by-step Figure 3 walkthrough
│   └── export_figures_code.r            # Batch figure export
│
├── calibration/          # Model calibration scripts
│   ├── calibrate_K_for_incidence.R
│   └── calibrate_R0_given_K.R
│
├── docs/                 # Documentation and guides
│   ├── trachoma_quickstart.md           # Getting started guide
│   ├── EXECUTIVE_SUMMARY.md             # Project overview and manuscript strategy
│   ├── FIGURE3_IMPLEMENTATION_GUIDE.md
│   ├── FIGURE3_STRATEGY_SUMMARY.md
│   └── FIGURE3_VISUAL_REFERENCE.md
│
├── documentation/        # Reference documents and notes
│
├── submission/           # Manuscript files for journal submission
│   ├── trachoma_transmission.rmd        # R Markdown manuscript
│   └── main/                            # Final submission figures
│
├── archive/              # Old versioned scripts (kept for reference)
│
├── results/              # Simulation outputs — NOT tracked in git
│   └── .gitkeep
│
├── trachoma-slim-export.yml    # Conda environment specification
└── trachoma_transmission_simulation.Rproj
```

---

## Getting Started

### 1. Set up the environment

The pipeline requires R, Python (via `reticulate`), SLiM (v5.0+), and the conda environment:

```bash
conda env create -f trachoma-slim-export.yml
conda activate trachoma-slim
```

R packages required: `slendr`, `igraph`, `ape`, `reticulate`, `tidyverse`, `ggplot2`

### 2. Run the pipeline

Open `trachoma_transmission_simulation.Rproj` in RStudio (sets working directory to project root), then:

```r
source("scripts/trachoma_pipeline_script.R")

results <- run_complete_pipeline(
  scenario   = "hyperendemic",
  N_hosts    = 1000,
  days       = 180,
  output_dir = "results/run_001"
)
```

### 3. Calibrate to a target incidence

```r
source("R/trachoma_calibration.R")
cal <- calibrate_model(target_incidence = 0.025)  # ~hyperendemic
```

### 4. Generate figures

Each figure has its own self-contained script in `figures/`. Run from the project root:

```r
source("figures/figure2_epidemic_model.R")
source("figures/figure3_methods_performance.R")
# etc.
```

---

## Simulation Design

| Parameter | Hypoendemic | Mesoendemic | Hyperendemic |
|-----------|-------------|-------------|--------------|
| Target incidence (per 100 PD) | < 0.003 | 0.005–0.01 | 0.02–0.03 |
| Genome length | 1,040,000 bp | ← | ← |
| Mutation rate | 5×10⁻⁷ | ← | ← |
| Bottleneck size | 1–10 strains | ← | ← |
| Critical inference window | — | — | ≤ 50 days |

---

## Dependency Map

```
trachoma_pipeline_script.R
  └── scripts/runme_trachoma_flat_two_sis_v1.R
        ├── slim/run_slim_haplo_dropin_v6.R
        ├── slim/finish_with_msprime_v6.R
        ├── analysis/pairwise_sanity_v6.R
        ├── analysis/roc_transmission_eval_time_v2.R
        ├── R/fn_diagnostics.R
        ├── R/trachoma_calibration.R
        └── R/trachoma_epidemic_core.R
  └── R/utils_environment.R
```

---

## Citation

If you use this framework, please cite:

> Diop MF et al. *Genetic inference of trachoma transmission: fundamental limitations and a path forward.* (in preparation)

---

## Contact

Mouhamadou Fadel Diop — diop.mfadel@gmail.com
London School of Hygiene and Tropical Medicine
