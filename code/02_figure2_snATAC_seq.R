# Figure 2: chromatin accessibility changes following CSC and aging.
# Refactored from All_figures_v12.R, original snATAC block.

source("code/00_setup.R")

out_dir <- ensure_dir(figure_path("figure2"))
derived_dir <- ensure_dir(results_path("figure2"))

integrated <- load_rdata_object(data_path("ATAC_integrated_final.RData"), "integrated")
integrated <- annotate_atac_clusters(integrated)

# Figure 2A-B: snATAC UMAPs.
fig2a <- plot_rpe_umap(integrated, groups = young_groups, ncol = 4)
save_panel(fig2a, file.path(out_dir, "Fig2A_young_snATAC_umap.pdf"), width = 12, height = 4)
save_panel(fig2a, file.path(out_dir, "Fig2A_young_snATAC_umap.png"), width = 12, height = 4)

fig2b <- plot_rpe_umap(integrated, groups = old_groups, ncol = 2)
save_panel(fig2b, file.path(out_dir, "Fig2B_aged_snATAC_umap.pdf"), width = 7, height = 4)
save_panel(fig2b, file.path(out_dir, "Fig2B_aged_snATAC_umap.png"), width = 7, height = 4)

load_total_fragments <- function() {
  standardize_fragment_table <- function(x, group_name = NULL) {
    x <- as.data.frame(x)
    if (!"Cell" %in% colnames(x)) {
      if ("CB" %in% colnames(x)) {
        x <- dplyr::rename(x, Cell = CB)
      } else {
        x <- tibble::rownames_to_column(x, "Cell")
      }
    }
    if (!is.null(group_name) && !grepl(paste0("^", group_name, "_"), x$Cell[[1]])) {
      x$Cell <- paste(group_name, x$Cell, sep = "_")
    }
    if (!"fragments" %in% colnames(x) && "frequency_count" %in% colnames(x)) {
      x <- dplyr::rename(x, fragments = frequency_count)
    }
    x %>%
      dplyr::mutate(Cell = gsub(".*\\.", "", Cell)) %>%
      dplyr::transmute(Cell, total_fragments = fragments)
  }

  total_frag_path <- first_existing_path(
    c(data_path("derived", "TotalFrag.RData"), data_path("TotalFrag.RData")),
    label = "TotalFrag.RData",
    must_work = FALSE
  )
  if (file.exists(total_frag_path)) {
    env <- new.env(parent = emptyenv())
    load(total_frag_path, envir = env)
    total_fragments <- get("total_fragments", envir = env)
    if (is.list(total_fragments)) {
      return(dplyr::bind_rows(lapply(names(total_fragments), function(name) {
        standardize_fragment_table(total_fragments[[name]], group_name = name)
      })))
    }
    return(standardize_fragment_table(total_fragments))
  }

  fragment_csv <- data_path("fragments", "fragment_paths.csv")
  if (!file.exists(fragment_csv)) {
    stop(
      "FRiP analysis needs either data/derived/TotalFrag.RData or data/fragments/fragment_paths.csv.",
      call. = FALSE
    )
  }
  fragment_paths <- readr::read_csv(fragment_csv, show_col_types = FALSE)
  totals <- lapply(seq_len(nrow(fragment_paths)), function(i) {
    Signac::CountFragments(fragment_paths$path[[i]]) %>%
      dplyr::mutate(Cell = paste(fragment_paths$group[[i]], CB, sep = "_"))
  })
  dplyr::bind_rows(totals)
}

# Figure 2C-D: fraction of fragments in peaks.
total_fragments <- load_total_fragments()
if (!"total_fragments" %in% colnames(total_fragments)) {
  stop("Total fragment table must contain or be normalized to `total_fragments`.", call. = FALSE)
}

integrated$peak_region_fragments <- Matrix::colSums(integrated@assays$peaks@counts)

frip <- integrated@meta.data %>%
  dplyr::mutate(cluster_group = paste(RNA_cluster, group, sep = "_")) %>%
  tibble::rownames_to_column("Cell") %>%
  dplyr::left_join(total_fragments, by = "Cell") %>%
  dplyr::mutate(
    FRiP = peak_region_fragments / total_fragments,
    group = factor(group, levels = c("con_3d", "cse_3d", "con_6d", "cse_6d", "con_10d", "cse_10d", "con_10d_old", "cse_10d_old")),
    rpe_state = dplyr::case_when(
      new_cluster == "healthy RPE" ~ "healthy RPE",
      new_cluster == "dedifferentiated RPE" & group == "cse_6d" ~ "dedifferentiated RPE CSC 6 days",
      new_cluster == "dedifferentiated RPE" & group == "cse_10d" ~ "dedifferentiated RPE CSC 10 days",
      new_cluster == "dedifferentiated RPE" & group == "con_10d_old" ~ "dedifferentiated RPE aged control",
      TRUE ~ "Other clusters"
    )
  )

fig2c <- frip %>%
  dplyr::filter(group %in% c("con_6d", "cse_6d", "con_10d", "cse_10d")) %>%
  ggplot(aes(x = rpe_state, y = FRiP, fill = rpe_state)) +
  geom_violin(trim = TRUE, scale = "width") +
  facet_wrap(~group, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "Fraction of fragments in peaks") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_panel(fig2c, file.path(out_dir, "Fig2C_young_FRiP.pdf"), width = 12, height = 4)
save_panel(fig2c, file.path(out_dir, "Fig2C_young_FRiP.png"), width = 12, height = 4)

fig2d <- frip %>%
  dplyr::filter(group %in% old_groups) %>%
  ggplot(aes(x = rpe_state, y = FRiP, fill = rpe_state)) +
  geom_violin(trim = TRUE, scale = "width") +
  facet_wrap(~group, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "Fraction of fragments in peaks") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_panel(fig2d, file.path(out_dir, "Fig2D_aged_FRiP.pdf"), width = 7, height = 4)
save_panel(fig2d, file.path(out_dir, "Fig2D_aged_FRiP.png"), width = 7, height = 4)

# Figure 2E: coverage tracks. Requires local fragment paths.
fragment_csv <- data_path("fragments", "fragment_paths.csv")
if (file.exists(fragment_csv) && file.exists(data_path("Handa_snRNA_Final.RData"))) {
  combined <- load_rdata_object(data_path("Handa_snRNA_Final.RData"), "combined")
  combined <- annotate_rna_clusters(combined)

  fragment_paths <- readr::read_csv(fragment_csv, show_col_types = FALSE)
  fragment_paths <- stats::setNames(fragment_paths$path, fragment_paths$group)

  meta_rna <- combined@meta.data %>%
    dplyr::rename(RNA_cluster = seurat_clusters) %>%
    dplyr::select(RNA_cluster, nCount_RNA, nFeature_RNA)

  Seurat::DefaultAssay(integrated) <- "ATAC"
  integrated_sub <- subset(integrated, subset = group %in% c("con_10d", "cse_10d", "con_10d_old"))
  integrated_sub <- Seurat::AddMetaData(integrated_sub, metadata = meta_rna)
  Seurat::DefaultAssay(integrated_sub) <- "peaks"
  integrated_sub <- Signac::RegionStats(integrated_sub, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  GenomeInfoDb::seqlevels(annotations) <- paste0("chr", GenomeInfoDb::seqlevels(annotations))
  GenomeInfoDb::genome(annotations) <- "mm10"
  Signac::Annotation(integrated_sub) <- annotations

  combined_sub <- combined[, colnames(combined) %in% colnames(integrated_sub)]
  integrated_sub[["RNA"]] <- combined_sub[["RNA"]]
  integrated_sub$celltype <- paste(integrated_sub$group, integrated_sub$RNA_cluster, sep = "_")
  Seurat::Idents(integrated_sub) <- "celltype"
  integrated_cov <- subset(integrated_sub, idents = c("cse_10d_1", "cse_10d_9", "con_10d_old_9"))

  fragment_list <- list()
  for (grp in intersect(names(fragment_paths), unique(integrated_cov$group))) {
    cells_grp <- colnames(integrated_cov)[integrated_cov$group == grp]
    fragment_list[[grp]] <- Signac::CreateFragmentObject(path = fragment_paths[[grp]], cells = cells_grp)
  }
  integrated_cov@assays$ATAC@fragments <- fragment_list

  for (gene in c("Abhd2", "Col9a3")) {
    coverage_plot <- Signac::CoveragePlot(
      object = integrated_cov,
      group.by = "celltype",
      region = gene,
      annotation = TRUE,
      extend.upstream = 500,
      extend.downstream = 5000,
      peaks = FALSE
    ) & ggplot2::scale_fill_manual(values = c(scales::hue_pal()(20)[2], scales::hue_pal()(20)[10], scales::hue_pal()(20)[10]))
    save_panel(coverage_plot, file.path(out_dir, paste0("Fig2E_coverage_", gene, ".pdf")), width = 5, height = 3)
    save_panel(coverage_plot, file.path(out_dir, paste0("Fig2E_coverage_", gene, ".png")), width = 5, height = 3)
  }
} else {
  message("Skipping Figure 2E coverage plots because fragment paths and/or RNA object are not available.")
}

# Figure 2F-H: differential accessibility and ATAC-explained DEGs.
Seurat::DefaultAssay(integrated) <- "peaks"
integrated$group_cluster <- paste(integrated$group, integrated$RNA_cluster, sep = "_")
Seurat::Idents(integrated) <- "group_cluster"

idents1 <- c("cse_6d_9", "cse_10d_9", "con_10d_old_9")
idents2 <- c("cse_6d_1", "cse_10d_1", "con_10d_old_1")
comparison_names <- c("Young6d_CSEC9vsCSEC1", "Young10d_CSEC9vsCSEC1", "Old10d_CONC9vsCONC1")

dap_avg <- Seurat::AverageExpression(integrated, assays = "peaks", group.by = "group_cluster")$peaks
dap <- vector("list", length(idents1))
names(dap) <- comparison_names

for (i in seq_along(idents1)) {
  dap[[i]] <- Seurat::FindMarkers(
    integrated,
    ident.1 = idents1[[i]],
    ident.2 = idents2[[i]],
    min.pct = 0.1,
    logfc.threshold = 0,
    test.use = "LR"
  ) %>%
    dplyr::mutate(
      group = comparison_names[[i]],
      dir = dplyr::if_else(avg_log2FC > 0, "Up-regulated", "Down-regulated")
    ) %>%
    tibble::rownames_to_column("query_region")

  dap[[i]] <- Signac::ClosestFeature(integrated, regions = dap[[i]]$query_region) %>%
    dplyr::left_join(dap[[i]], by = "query_region") %>%
    dplyr::mutate(
      p_val_adj = p.adjust(p_val, method = "BH"),
      sig = dplyr::if_else(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1, "FDR<0.05", "No-sig")
    )
}

volcano_plots <- lapply(dap, function(x) {
  ggplot(x, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
    geom_point(size = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("FDR<0.05" = "red", "No-sig" = "black")) +
    labs(x = "log2 fold change", y = "-log10 adjusted p value") +
    paper_theme() +
    theme(legend.position = "none")
})
fig2f <- cowplot::plot_grid(plotlist = volcano_plots, nrow = 1, labels = names(volcano_plots))
save_panel(fig2f, file.path(out_dir, "Fig2F_DAP_volcano.pdf"), width = 12, height = 4)
save_panel(fig2f, file.path(out_dir, "Fig2F_DAP_volcano.png"), width = 12, height = 4)

gene_ranges <- ensembldb::genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
gene_ranges <- gene_ranges[gene_ranges$gene_biotype == "protein_coding"]
gene_ranges <- GenomeInfoDb::keepStandardChromosomes(gene_ranges, pruning.mode = "coarse")
GenomeInfoDb::seqlevelsStyle(gene_ranges) <- "UCSC"
gene_promoters <- GenomicRanges::promoters(gene_ranges)

dap_promoter <- lapply(dap, function(dap_i) {
  peaks <- Signac::StringToGRanges(regions = dap_i$query_region, sep = c(":", "-"))
  promoter_hits <- IRanges::findOverlaps(query = peaks, subject = gene_promoters)
  dap_i[S4Vectors::queryHits(promoter_hits), ] %>%
    dplyr::select(query_region, avg_log2FC, p_val, pct.1, pct.2, p_val_adj, group) %>%
    dplyr::bind_cols(as.data.frame(gene_promoters[S4Vectors::subjectHits(promoter_hits)]) %>% dplyr::select(gene_name)) %>%
    dplyr::mutate(sig = dplyr::if_else(p_val_adj < 0.05, "FDR<0.05", "No-sig"))
})

deg_path <- first_existing_path(
  c(data_path("derived", "DEG_sickRPE_remove.xlsx"), data_path("DEG_sickRPE_remove.xlsx")),
  label = "DEG_sickRPE_remove.xlsx",
  must_work = FALSE
)
if (file.exists(deg_path)) {
  deg_sick_rpe <- read_excel_allsheets(deg_path)
  deg_dap <- vector("list", length(dap_promoter))
  names(deg_dap) <- names(dap_promoter)
  for (i in seq_along(dap_promoter)) {
    deg_dap[[i]] <- dap_promoter[[i]] %>%
      dplyr::rename(Genename = gene_name, atac_log2fc = avg_log2FC, atac_p_val_adj = p_val_adj) %>%
      dplyr::inner_join(deg_sick_rpe[[comparison_names[[i]]]], by = "Genename") %>%
      dplyr::filter(
        atac_p_val_adj < 0.05,
        p_val_adj < 0.05,
        (atac_log2fc > 0 & avg_log2FC > 0) | (atac_log2fc < 0 & avg_log2FC < 0)
      ) %>%
      dplyr::distinct(Genename, .keep_all = TRUE)
  }

  deg_dap_summary <- tibble::tibble(
    group = comparison_names,
    atac_explained_deg = vapply(deg_dap, nrow, integer(1)),
    total_deg = vapply(deg_sick_rpe[comparison_names], function(x) sum(x$p_val_adj < 0.05 & abs(x$avg_log2FC) > 0.1), integer(1))
  ) %>%
    dplyr::mutate(percent = 100 * atac_explained_deg / total_deg)

  writexl::write_xlsx(deg_dap, file.path(derived_dir, "DAPexplainableDEG.xlsx"))

  fig2g <- ggplot(deg_dap_summary, aes(x = group, y = percent, fill = group)) +
    geom_col(width = 0.65, color = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(60, max(deg_dap_summary$percent, na.rm = TRUE)))) +
    labs(x = NULL, y = "DEGs with same-direction promoter DAP (%)") +
    paper_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  save_panel(fig2g, file.path(out_dir, "Fig2G_ATAC_explained_DEG_percent.pdf"), width = 5, height = 4)
  save_panel(fig2g, file.path(out_dir, "Fig2G_ATAC_explained_DEG_percent.png"), width = 5, height = 4)

  deg_dap_sig <- lapply(deg_dap, function(x) x %>% dplyr::select(Genename, avg_log2FC, p_val_adj, group, dir))
  universe_symbols <- if ("RNA" %in% names(integrated@assays)) {
    rownames(integrated@assays$RNA)
  } else {
    unique(unlist(lapply(deg_sick_rpe, `[[`, "Genename")))
  }
  kegg <- run_kegg_by_group(deg_dap_sig, universe_symbols = universe_symbols)
  openxlsx::write.xlsx(kegg, file.path(derived_dir, "KEGG_DAPexplainableDEGs.xlsx"), overwrite = TRUE)

  kegg_plot_data <- dplyr::bind_rows(lapply(names(kegg), function(name) {
    kegg[[name]] %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::mutate(group = name, Description = gsub(" - Mus.*", "", Description), logP = -log10(P.DE))
  }))

  if (nrow(kegg_plot_data) > 0) {
    fig2h <- ggplot(kegg_plot_data, aes(x = group, y = reorder(Description, logP), size = logP, fill = group)) +
      geom_point(shape = 21) +
      labs(x = NULL, y = NULL, size = "-log10(P)") +
      paper_theme()
    save_panel(fig2h, file.path(out_dir, "Fig2H_ATAC_explained_DEG_KEGG.pdf"), width = 8, height = 6)
    save_panel(fig2h, file.path(out_dir, "Fig2H_ATAC_explained_DEG_KEGG.png"), width = 8, height = 6)
  }
} else {
  message("Skipping Figure 2G-H because data/derived/DEG_sickRPE_remove.xlsx is not available.")
}

message("Figure 2 panels written to: ", out_dir)
