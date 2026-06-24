# Code Review Notes

These notes summarize the cleanup performed when refactoring `All_figures_v12.R` into the public figure scripts.

## Main Issues Found

- Published figure numbering and script comments were inconsistent. In the published paper, Figure 2 is the snATAC/accessibility figure and Figure 3 is the transcriptional DEG figure. The original script labeled the DEG block as Figure 2 and the snATAC block as Figure 3.
- The original script mixed all figures, supplementary analyses, exports, and exploratory commands in one 2,448-line file. The release now uses one script per published figure plus shared helpers.
- Several absolute local/server paths were embedded in the analysis, including fragment paths and manuscript export locations. The refactor replaces these with `MOUSE_CSC_DATA_DIR`, `MOUSE_CSC_FIGURE_DIR`, and documented local data layout.
- Interactive calls such as `View()` were removed because they break non-interactive runs.
- Some objects were used before being defined in the Figure 1-3 blocks, including `Marker_EMT`, `mouse_gene_SenMayo`, `DEG_sickRPE_remove_sig`, and `All_shared`. The refactor loads gene sets explicitly and computes required derived objects locally.
- A malformed duplicate block near the original Figure 3 KEGG export would not parse cleanly as a standalone script. The refactor rewrites the ATAC-explained DEG and KEGG steps.
- Coverage-plot code used `genes.use = "Rpe65"` while plotting `region = "Ctsl"` and saving to a `Basp1` filename. The refactor makes the coverage genes explicit and aligned with the published Figure 2 caption: `Abhd2` and `Col9a3`.
- The public repository excludes large/private inputs and generated artifacts via `.gitignore`.

## Remaining Assumptions

- The Seurat object `combined` is present in `Handa_snRNA_Final.RData`.
- The Signac/Seurat object `integrated` is present in `ATAC_integrated_final.RData`.
- The DEG workbook keeps the same sheet names used in the original analysis: `Young6d_CSEC9vsCSEC1`, `Young10d_CSEC9vsCSEC1`, and `Old10d_CONC9vsCONC1`.
- Fragment paths are supplied locally if fragment-dependent calculations are rebuilt.
