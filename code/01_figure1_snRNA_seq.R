# Figure 1: snRNA-seq of control and CSC-treated mice.
# Refactored from All_figures_v12.R, original Figure 1 block.

source("code/00_setup.R")

out_dir <- ensure_dir(figure_path("figure1"))
combined <- load_rdata_object(data_path("Handa_snRNA_Final.RData"), "combined")
combined <- annotate_rna_clusters(combined)
Seurat::DefaultAssay(combined) <- "RNA"

# Figure 1A-B: RPE and dedifferentiated RPE UMAPs.
fig1a <- plot_rpe_umap(combined, groups = young_groups, ncol = 4)
save_panel(fig1a, file.path(out_dir, "Fig1A_young_snRNA_umap.pdf"), width = 12, height = 4)
save_panel(fig1a, file.path(out_dir, "Fig1A_young_snRNA_umap.png"), width = 12, height = 4)

fig1b <- plot_rpe_umap(combined, groups = old_groups, ncol = 2)
save_panel(fig1b, file.path(out_dir, "Fig1B_aged_snRNA_umap.pdf"), width = 7, height = 4)
save_panel(fig1b, file.path(out_dir, "Fig1B_aged_snRNA_umap.png"), width = 7, height = 4)

# Figure 1C: expression of RPE marker genes highlighted in the paper.
marker_genes <- c("Slc6a20a", "Myrip")
marker_expression <- fetch_expression_long(combined, marker_genes) %>%
  dplyr::filter(new_cluster != "Other clusters") %>%
  dplyr::mutate(
    new_group = dplyr::if_else(
      new_cluster == "dedifferentiated RPE",
      paste(group, new_cluster, sep = "_"),
      as.character(new_cluster)
    ),
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

fig1c <- ggplot(marker_expression, aes(x = new_group, y = expression, fill = new_group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.75, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.45) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(scales::hue_pal()(19)[2], rep(scales::hue_pal()(19)[10], 3))) +
  labs(x = NULL, y = "Average expression level") +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_panel(fig1c, file.path(out_dir, "Fig1C_RPE_marker_expression.pdf"), width = 5, height = 7)
save_panel(fig1c, file.path(out_dir, "Fig1C_RPE_marker_expression.png"), width = 5, height = 7)

# Figure 1D: percentage of dedifferentiated RPE among total RPE cells.
cluster_sum_rpe <- combined@meta.data %>%
  dplyr::filter(as.character(seurat_clusters) %in% c("1", "9")) %>%
  dplyr::group_by(seurat_clusters, group) %>%
  dplyr::mutate(cell_no = dplyr::n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(total_cell = dplyr::n(), percentage = cell_no / total_cell) %>%
  dplyr::ungroup() %>%
  dplyr::select(seurat_clusters, group, cell_no, total_cell, percentage) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    rpe_state = dplyr::if_else(as.character(seurat_clusters) == "1", "healthy RPE", "dedifferentiated RPE"),
    group = factor(group, levels = c("cse_6d", "cse_10d", "con_10d_old"))
  ) %>%
  dplyr::filter(group %in% c("cse_6d", "cse_10d", "con_10d_old"))

fig1d <- ggplot(cluster_sum_rpe, aes(x = group, y = percentage, fill = rpe_state)) +
  geom_col(width = 0.75, color = "black") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c("healthy RPE" = scales::hue_pal()(19)[2], "dedifferentiated RPE" = scales::hue_pal()(19)[10])) +
  labs(x = NULL, y = "Percentage of total RPE cells") +
  paper_theme()
save_panel(fig1d, file.path(out_dir, "Fig1D_dedifferentiated_RPE_percentage.pdf"), width = 5, height = 4)
save_panel(fig1d, file.path(out_dir, "Fig1D_dedifferentiated_RPE_percentage.png"), width = 5, height = 4)

message("Figure 1 panels written to: ", out_dir)
