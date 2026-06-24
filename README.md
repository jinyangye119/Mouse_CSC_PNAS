# Mouse_CSC_PNAS

Figure-generation code for:

> Singh, Jin, Hu, Palazzo, Cano, Hoang, Bhutto, Wang, Sinha, Blackshaw, Qian, and Handa. Molecular underpinnings of induced degenerative heterogeneity in the retinal pigment epithelium. PNAS, 2026.

This repository organizes the analysis code used to generate Figures 1-3 of the paper. The scripts were refactored from the working analysis file `All_figures_v12.R` into one script per published figure.

## Repository Layout

- `code/00_setup.R` - shared package loading, paths, annotation tables, and plotting helpers.
- `code/01_figure1_snRNA_seq.R` - Figure 1: snRNA-seq RPE/dedifferentiated RPE panels.
- `code/02_figure2_snATAC_seq.R` - Figure 2: snATAC-seq/accessibility panels.
- `code/03_figure3_transcriptional_phenotype.R` - Figure 3: DEG, KEGG, EMT/senescence/mitochondrial panels.
- `code/run_all_figures.R` - runs the three figure scripts in order.
- `data/README.md` - expected input files and where to place them.
- `CODE_REVIEW.md` - notes from refactoring and issues fixed from the original script.

## Quick Start

Run from the repository root:

```r
source("code/run_all_figures.R")
```

or run one figure at a time:

```bash
Rscript code/01_figure1_snRNA_seq.R
Rscript code/02_figure2_snATAC_seq.R
Rscript code/03_figure3_transcriptional_phenotype.R
```

By default, scripts look for inputs under `data/` and write outputs under `figures/`. You can keep large data outside git and point the scripts to it:

```bash
export MOUSE_CSC_DATA_DIR="/path/to/Mouse_CSC_PNAS_data"
export MOUSE_CSC_FIGURE_DIR="/path/to/output_figures"
Rscript code/run_all_figures.R
```

## Notes

The repository intentionally excludes large `.RData`, `.rds`, fragment, and generated figure files. See `data/README.md` for the expected local data layout.
