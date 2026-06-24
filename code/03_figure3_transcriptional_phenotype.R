# Figure 3: transcriptional phenotype of dedifferentiated RPE.
# Refactored from All_figures_v12.R, original DEG/KEGG/gene-set expression block.

source("code/00_setup.R")

out_dir <- ensure_dir(figure_path("figure3"))
derived_dir <- ensure_dir(results_path("figure3"))

combined <- load_rdata_object(data_path("Handa_snRNA_Final.RData"), "combined")
combined <- annotate_rna_clusters(combined)
Seurat::DefaultAssay(combined) <- "RNA"

comparison_names <- c("Young6d_CSEC9vsCSEC1", "Young10d_CSEC9vsCSEC1", "Old10d_CONC9vsCONC1")

compute_deg_list <- function(object) {
  object_deg <- Seurat::SCTransform(object)
  Seurat::Idents(object_deg) <- "seurat_clusters"
  Seurat::DefaultAssay(object_deg) <- "SCT"
  object_deg$celltype <- paste(object_deg$seurat_clusters, object_deg$group, sep = "_")
  Seurat::Idents(object_deg) <- "celltype"
  object_deg <- Seurat::PrepSCTFindMarkers(object_deg)

  comparisons <- list(
    Young6d_CSEC9vsCSEC1 = c("9_cse_6d", "1_cse_6d"),
    Young10d_CSEC9vsCSEC1 = c("9_cse_10d", "1_cse_10d"),
    Old10d_CONC9vsCONC1 = c("9_con_10d_old", "1_con_10d_old")
  )

  lapply(names(comparisons), function(name) {
    ids <- comparisons[[name]]
    Seurat::FindMarkers(object_deg, assay = "SCT", ident.1 = ids[[1]], ident.2 = ids[[2]], logfc.threshold = 0) %>%
      dplyr::mutate(
        group = name,
        dir = dplyr::if_else(avg_log2FC > 0, "Up-regulated", "Down-regulated")
      ) %>%
      tibble::rownames_to_column("Genename") %>%
      dplyr::mutate(sig = dplyr::if_else(abs(avg_log2FC) > 0.1 & p_val_adj < 0.05, "fdr<0.05", "non-sig")) %>%
      dplyr::filter(!Genename == "AY036118", !grepl("^mt-", Genename), !grepl("^Gm", Genename))
  }) %>%
    stats::setNames(names(comparisons))
}

deg_path <- first_existing_path(
  c(data_path("derived", "DEG_sickRPE_remove.xlsx"), data_path("DEG_sickRPE_remove.xlsx")),
  label = "DEG_sickRPE_remove.xlsx",
  must_work = FALSE
)
if (file.exists(deg_path)) {
  deg_sick_rpe <- read_excel_allsheets(deg_path)
} else {
  deg_sick_rpe <- compute_deg_list(combined)
  writexl::write_xlsx(deg_sick_rpe, file.path(derived_dir, "DEG_sickRPE_remove.xlsx"))
}

deg_sick_rpe <- deg_sick_rpe[comparison_names]
deg_sig <- lapply(deg_sick_rpe, function(x) {
  x %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.1)
})

# Figure 3A: DEG counts.
deg_counts <- dplyr::bind_rows(deg_sig, .id = "comparison") %>%
  dplyr::count(comparison, dir, name = "n") %>%
  dplyr::mutate(comparison = factor(comparison, levels = comparison_names))

fig3a <- ggplot(deg_counts, aes(x = comparison, y = n, fill = dir)) +
  geom_col(position = "dodge", width = 0.75, color = "black") +
  scale_fill_manual(values = c("Up-regulated" = "#e41a1c", "Down-regulated" = "#377eb8")) +
  labs(x = NULL, y = "Number of DEGs") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_panel(fig3a, file.path(out_dir, "Fig3A_DEG_counts.pdf"), width = 6, height = 4)
save_panel(fig3a, file.path(out_dir, "Fig3A_DEG_counts.png"), width = 6, height = 4)

# Figure 3B: volcano plots.
volcano_plots <- lapply(deg_sick_rpe, function(deg) {
  label_data <- deg %>%
    dplyr::filter(sig == "fdr<0.05", Genename %in% c("Ankfn1", "Ctsl"))

  ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
    geom_point(size = 0.5) +
    ggrepel::geom_text_repel(data = label_data, aes(label = Genename), max.overlaps = Inf) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("fdr<0.05" = "red", "non-sig" = "black")) +
    labs(x = "Log2 fold change in dedifferentiated RPE vs healthy RPE", y = "-log10 adjusted p value") +
    paper_theme() +
    theme(legend.position = "none")
})
fig3b <- cowplot::plot_grid(plotlist = volcano_plots, nrow = 1, labels = names(volcano_plots))
save_panel(fig3b, file.path(out_dir, "Fig3B_DEG_volcano.pdf"), width = 12, height = 4)
save_panel(fig3b, file.path(out_dir, "Fig3B_DEG_volcano.png"), width = 12, height = 4)

# Figure 3C: KEGG enrichment for DEGs.
universe_symbols <- rownames(combined@assays$RNA)
deg_for_kegg <- split(
  dplyr::bind_rows(deg_sig, .id = "comparison"),
  list(
    dplyr::bind_rows(deg_sig, .id = "comparison")$comparison,
    dplyr::bind_rows(deg_sig, .id = "comparison")$dir
  ),
  drop = TRUE
)
kegg <- run_kegg_by_group(deg_for_kegg, universe_symbols = universe_symbols)
kegg_for_export <- kegg
names(kegg_for_export) <- make.unique(substr(gsub("[^A-Za-z0-9_]", "_", names(kegg_for_export)), 1, 31))
openxlsx::write.xlsx(kegg_for_export, file.path(derived_dir, "KEGG_allGroups_Fig3C.xlsx"), overwrite = TRUE)

kegg_plot_data <- dplyr::bind_rows(lapply(names(kegg), function(name) {
  kegg[[name]] %>%
    dplyr::filter(!Description %in% c("Ribosome", "Coronavirus disease - COVID-19", "Metabolic pathways", "MicroRNAs in cancer", "Phototransduction")) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::mutate(
      group = gsub("\\..*", "", name),
      direction = gsub(".*\\.", "", name),
      group_direction = paste(direction, group, sep = "_"),
      Description = gsub(" - Mus.*", "", Description),
      logP = -log10(padj)
    )
}))

fig3c <- ggplot(kegg_plot_data, aes(x = Description, y = group_direction, color = group_direction)) +
  geom_point(aes(size = logP)) +
  coord_flip() +
  labs(x = NULL, y = NULL, size = "-log10(FDR)") +
  paper_theme() +
  theme(panel.background = element_blank())
save_panel(fig3c, file.path(out_dir, "Fig3C_DEG_KEGG.pdf"), width = 9, height = 7)
save_panel(fig3c, file.path(out_dir, "Fig3C_DEG_KEGG.png"), width = 9, height = 7)

# Figure 3D: EMT, SenMayo, and MitoCarta gene-set expression.
load_gene_set <- function(path, candidates) {
  if (!file.exists(path)) {
    stop("Gene-set file not found: ", path, call. = FALSE)
  }
  data <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    readxl::read_excel(path)
  }
  col <- first_existing_col(data, candidates)
  unique(stats::na.omit(data[[col]]))
}

emt_path <- first_existing_path(
  c(data_path("gene_sets", "EMT_genes.xlsx"), data_path("../../../raw/meta/EMT_genes.xlsx")),
  label = "EMT gene set"
)
senmayo_path <- first_existing_path(
  c(
    data_path("gene_sets", "125_SenMayo_mouse.xlsx"),
    data_path("gene_sets", "125_SenMayogenes.xlsx"),
    data_path("../../../raw/meta/125_SenMayogenes.xlsx")
  ),
  label = "SenMayo gene set"
)
mitocarta_path <- first_existing_path(
  c(data_path("gene_sets", "Mouse.MitoCarta3.0.xls"), data_path("../../../raw/meta/Mouse.MitoCarta3.0.xls")),
  label = "Mouse MitoCarta gene set"
)

emt_genes <- load_gene_set(emt_path, c("Gene_name", "Mouse", "mouse", "Symbol", "Gene"))
senmayo_genes <- load_gene_set(senmayo_path, c("Gene(murine)", "Mouse", "mouse", "Gene_name", "Symbol", "Gene"))
mitocarta_data <- readxl::read_excel(mitocarta_path, sheet = 2)
mitocarta_genes <- unique(stats::na.omit(mitocarta_data[[first_existing_col(mitocarta_data, c("Symbol", "Gene", "GeneSymbol"))]]))

gene_set_expression <- dplyr::bind_rows(
  summarize_gene_set_expression(combined, emt_genes) %>% dplyr::mutate(gene_set = "EMT genes"),
  summarize_gene_set_expression(combined, senmayo_genes) %>% dplyr::mutate(gene_set = "SenMayo genes"),
  summarize_gene_set_expression(combined, mitocarta_genes) %>% dplyr::mutate(gene_set = "MitoCarta genes")
) %>%
  dplyr::mutate(
    new_group = factor(
      new_group,
      levels = c(
        "healthy RPE",
        "cse_6d_dedifferentiated RPE",
        "cse_10d_dedifferentiated RPE",
        "con_10d_old_dedifferentiated RPE"
      )
    )
  )

fig3d <- ggplot(gene_set_expression, aes(x = new_group, y = total_expression, fill = new_group)) +
  geom_violin(trim = TRUE, scale = "width") +
  facet_wrap(~gene_set, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(scales::hue_pal()(19)[2], rep(scales::hue_pal()(19)[10], 3))) +
  labs(x = NULL, y = "Summed expression") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_panel(fig3d, file.path(out_dir, "Fig3D_gene_set_expression.pdf"), width = 12, height = 4)
save_panel(fig3d, file.path(out_dir, "Fig3D_gene_set_expression.png"), width = 12, height = 4)

message("Figure 3 panels written to: ", out_dir)
