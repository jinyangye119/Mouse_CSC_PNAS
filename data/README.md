# Data Layout

The figure scripts do not commit large input data. Place files in this directory, or set `MOUSE_CSC_DATA_DIR` to another directory with the same layout.

## Required Objects

```text
data/
  Handa_snRNA_Final.RData
  ATAC_integrated_final.RData
  derived/
    DEG_sickRPE_remove.xlsx
    TotalFrag.RData
    DAPexpainDEG.xlsx
  gene_sets/
    New_rpe_marker.csv
    EMT_genes.xlsx
    125_SenMayo_mouse.xlsx
    Mouse.MitoCarta3.0.xls
  fragments/
    fragment_paths.csv
```

## Object Names Expected

- `Handa_snRNA_Final.RData` should contain a Seurat object named `combined`.
- `ATAC_integrated_final.RData` should contain a Signac/Seurat object named `integrated`.
- `TotalFrag.RData` should contain a data frame or list named `total_fragments` with cell barcodes and total fragments.

## Fragment Path CSV

`data/fragments/fragment_paths.csv` is only needed if `TotalFrag.RData` is not already available or if coverage tracks need to rebuild fragment objects. It should contain:

```text
group,path
con_3d,/path/to/con_3d/atac_fragments.tsv.gz
cse_3d,/path/to/cse_3d/atac_fragments.tsv.gz
con_6d,/path/to/con_6d/atac_fragments.tsv.gz
cse_6d,/path/to/cse_6d/atac_fragments.tsv.gz
con_10d,/path/to/con_10d/atac_fragments.tsv.gz
cse_10d,/path/to/cse_10d/atac_fragments.tsv.gz
con_10d_old,/path/to/con_10d_old/atac_fragments.tsv.gz
cse_10d_old,/path/to/cse_10d_old/atac_fragments.tsv.gz
```

## Generated Output

Scripts write figure panels to:

```text
figures/
  figure1/
  figure2/
  figure3/
```
