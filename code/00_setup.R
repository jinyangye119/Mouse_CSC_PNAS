# Shared setup for Mouse CSC PNAS Figure 1-3 code.

required_packages <- c(
  "Seurat", "Signac", "tidyverse", "cowplot", "scales", "readxl",
  "writexl", "openxlsx", "ggrepel", "limma", "AnnotationDbi",
  "org.Mm.eg.db", "EnsDb.Mmusculus.v79", "BSgenome.Mmusculus.UCSC.mm10",
  "Matrix", "ensembldb", "GenomeInfoDb", "GenomicRanges", "IRanges", "S4Vectors"
)

load_required_packages <- function(packages = required_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "),
      "\nInstall CRAN packages with install.packages() and Bioconductor packages with BiocManager::install().",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    invisible(lapply(packages, library, character.only = TRUE))
  })
}

load_required_packages()

project_root <- normalizePath(getwd(), mustWork = TRUE)
data_dir <- normalizePath(Sys.getenv("MOUSE_CSC_DATA_DIR", file.path(project_root, "data")), mustWork = FALSE)
figure_dir <- normalizePath(Sys.getenv("MOUSE_CSC_FIGURE_DIR", file.path(project_root, "figures")), mustWork = FALSE)
results_dir <- normalizePath(Sys.getenv("MOUSE_CSC_RESULTS_DIR", file.path(project_root, "results")), mustWork = FALSE)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

data_path <- function(...) file.path(data_dir, ...)
figure_path <- function(...) file.path(figure_dir, ...)
results_path <- function(...) file.path(results_dir, ...)

first_existing_path <- function(candidates, label = "file", must_work = TRUE) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  if (must_work) {
    stop("Could not find ", label, ". Checked:\n", paste(candidates, collapse = "\n"), call. = FALSE)
  }
  candidates[[1]]
}

load_rdata_object <- function(path, object_name) {
  if (!file.exists(path)) {
    stop("Required data file not found: ", path, call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object `", object_name, "` was not found in ", path, call. = FALSE)
  }
  get(object_name, envir = env, inherits = FALSE)
}

read_excel_allsheets <- function(path) {
  if (!file.exists(path)) {
    stop("Workbook not found: ", path, call. = FALSE)
  }
  sheets <- readxl::excel_sheets(path)
  out <- lapply(sheets, function(sheet) as.data.frame(readxl::read_excel(path, sheet = sheet)))
  names(out) <- sheets
  out
}

first_existing_col <- function(data, candidates) {
  hit <- candidates[candidates %in% colnames(data)]
  if (length(hit) == 0) {
    stop("None of these columns were found: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[[1]]
}

save_panel <- function(plot, filename, width = 6, height = 4, dpi = 300) {
  ensure_dir(dirname(filename))
  if (grepl("[.]pdf$", filename, ignore.case = TRUE)) {
    ggplot2::ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, useDingbats = FALSE)
  } else {
    ggplot2::ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
  }
}

cell_annotation_table <- function(cluster_col = "seurat_clusters") {
  tibble::tibble(
    "{cluster_col}" := as.character(0:19),
    cell_anno = c(
      "0: Stromal cells -1",
      "1: RPE",
      "2: Stromal cells -2",
      "3: Melanocytes",
      "4: Rod -1",
      "5: Smooth Muscle",
      "6: Microglia",
      "7: V/E cells",
      "8: Stromal cells -3",
      "9: Dediff RPE",
      "10: Fibroblast",
      "11: non-pigmented ciliary epithelial cells",
      "12: Rod -2",
      "13: Stromal cells -4",
      "14: Amacrine Cells",
      "15: Schwann cell",
      "16: Cone",
      "17: Muller glia",
      "18: Cillary body cell",
      "19: Bipolar Cells"
    )
  )
}

annotate_rna_clusters <- function(object) {
  meta <- cell_annotation_table("seurat_clusters")
  object@meta.data <- object@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    left_join(meta, by = "seurat_clusters") %>%
    tibble::column_to_rownames("cell") %>%
    mutate(
      seurat_clusters = factor(seurat_clusters, levels = as.character(0:19)),
      cell_anno = factor(cell_anno, levels = meta$cell_anno),
      new_cluster = dplyr::case_when(
        as.character(seurat_clusters) == "1" ~ "healthy RPE",
        as.character(seurat_clusters) == "9" ~ "dedifferentiated RPE",
        TRUE ~ "Other clusters"
      ),
      new_cluster = factor(new_cluster, levels = c("healthy RPE", "dedifferentiated RPE", "Other clusters"))
    )
  object
}

annotate_atac_clusters <- function(object) {
  meta <- cell_annotation_table("RNA_cluster")
  object@meta.data <- object@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(RNA_cluster = as.character(RNA_cluster)) %>%
    left_join(meta, by = "RNA_cluster") %>%
    tibble::column_to_rownames("cell") %>%
    mutate(
      RNA_cluster = factor(RNA_cluster, levels = as.character(0:19)),
      cell_anno = factor(cell_anno, levels = meta$cell_anno),
      new_cluster = dplyr::case_when(
        as.character(RNA_cluster) == "1" ~ "healthy RPE",
        as.character(RNA_cluster) == "9" ~ "dedifferentiated RPE",
        TRUE ~ "Other clusters"
      ),
      new_cluster = factor(new_cluster, levels = c("healthy RPE", "dedifferentiated RPE", "Other clusters"))
    )
  object
}

group_labels <- c(
  con_3d = "Control 3 days",
  cse_3d = "CSC 3 days",
  con_6d = "Control 6 days",
  cse_6d = "CSC 6 days",
  con_10d = "Control",
  cse_10d = "CSC 10 days",
  con_10d_old = "Aged control",
  cse_10d_old = "Aged CSC 10 days"
)

young_groups <- c("con_10d", "cse_3d", "cse_6d", "cse_10d")
old_groups <- c("con_10d_old", "cse_10d_old")
rpe_colors <- c("healthy RPE" = scales::hue_pal()(19)[2], "dedifferentiated RPE" = scales::hue_pal()(19)[10], "Other clusters" = "grey70")

paper_theme <- function(base_size = 14) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(colour = "black", linewidth = 0.6),
      strip.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_blank()
    )
}

plot_rpe_umap <- function(object, groups, ncol = NULL) {
  object_sub <- subset(object, subset = group %in% groups)
  object_sub$group_label <- factor(unname(group_labels[as.character(object_sub$group)]), levels = unname(group_labels[groups]))
  Seurat::DimPlot(
    object = object_sub,
    group.by = "new_cluster",
    reduction = "umap",
    split.by = "group_label",
    ncol = ncol
  ) &
    paper_theme() &
    ggplot2::scale_color_manual(values = rpe_colors) &
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 19, size = 4))) &
    ggplot2::ggtitle("")
}

fetch_expression_long <- function(object, genes, metadata_cols = c("group", "new_cluster"), assay = "RNA") {
  Seurat::DefaultAssay(object) <- assay
  genes <- intersect(genes, rownames(object[[assay]]))
  if (length(genes) == 0) {
    stop("None of the requested genes were found in assay `", assay, "`.", call. = FALSE)
  }
  Seurat::FetchData(object, vars = c(genes, metadata_cols)) %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(genes), names_to = "gene", values_to = "expression")
}

summarize_gene_set_expression <- function(object, genes, assay = "RNA") {
  fetch_expression_long(object, genes, assay = assay) %>%
    dplyr::filter(new_cluster != "Other clusters") %>%
    dplyr::mutate(
      new_group = dplyr::if_else(
        new_cluster == "dedifferentiated RPE",
        paste(group, new_cluster, sep = "_"),
        as.character(new_cluster)
      )
    ) %>%
    dplyr::group_by(cell, group, new_cluster, new_group) %>%
    dplyr::summarise(total_expression = sum(expression), .groups = "drop")
}

map_symbols_to_entrez <- function(symbols) {
  AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
}

get_kegg_reference <- function() {
  list(
    gene_links = limma::getGeneKEGGLinks(species.KEGG = "mmu") %>%
      dplyr::mutate(
        Gene_Des = AnnotationDbi::mapIds(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = GeneID,
          column = "GENENAME",
          keytype = "ENTREZID",
          multiVals = "first"
        ),
        Genename = AnnotationDbi::mapIds(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = GeneID,
          column = "SYMBOL",
          keytype = "ENTREZID",
          multiVals = "first"
        )
      ),
    pathway_names = limma::getKEGGPathwayNames(species.KEGG = "mmu")
  )
}

run_kegg_by_group <- function(deg_list, universe_symbols, p_cutoff = 0.05) {
  kegg_ref <- get_kegg_reference()
  universe <- tibble::tibble(Genename = universe_symbols, entrez = map_symbols_to_entrez(universe_symbols)) %>%
    dplyr::filter(!is.na(entrez))

  out <- lapply(deg_list, function(deg) {
    deg <- deg %>%
      dplyr::mutate(entrez = map_symbols_to_entrez(Genename)) %>%
      dplyr::filter(!is.na(entrez))

    limma::kegga(unique(deg$entrez), species.KEGG = "mmu", universe = universe$entrez) %>%
      tibble::rownames_to_column("PathwayID") %>%
      dplyr::mutate(
        padj = p.adjust(P.DE, method = "BH"),
        PathwayID = gsub("path:", "", PathwayID)
      ) %>%
      dplyr::left_join(kegg_ref$pathway_names, by = "PathwayID") %>%
      dplyr::arrange(P.DE) %>%
      dplyr::filter(padj < p_cutoff)
  })
  names(out) <- names(deg_list)
  out
}
