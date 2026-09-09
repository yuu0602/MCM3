#!/usr/bin/env Rscript

# Final RNA-seq analysis.  Inputs are Salmon quantifications generated in
# reanalysis(final); it intentionally never reads prior DEG result tables.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(tximport); library(readr); library(dplyr); library(tibble)
  library(edgeR); library(limma); library(ggplot2); library(pheatmap); library(grid)
  library(AnnotationDbi); library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: _rnaseq_deg.R <deg_work> <reference>")
ROOT <- normalizePath(args[[1]], mustWork = TRUE)
REF <- normalizePath(args[[2]], mustWork = TRUE)
MANIFEST <- file.path(ROOT, "metadata", "Salmon_quantifications.tsv")
OUT_DATA <- file.path(ROOT, "data")
OUT_VIS <- file.path(ROOT, "visuals")
COMPAT_ROOT <- file.path(ROOT, "limma_outputs")
COMPAT_DEG <- file.path(COMPAT_ROOT, "DEGs")
COMPAT_VOLCANO <- file.path(COMPAT_ROOT, "volcanos")
COMPAT_HEATMAP <- file.path(COMPAT_ROOT, "heatmaps")
dir.create(OUT_DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_VIS, recursive = TRUE, showWarnings = FALSE)
dir.create(COMPAT_DEG, recursive = TRUE, showWarnings = FALSE)
dir.create(COMPAT_VOLCANO, recursive = TRUE, showWarnings = FALSE)
dir.create(COMPAT_HEATMAP, recursive = TRUE, showWarnings = FALSE)

FDR_MAX <- 0.05
LFC_MIN <- 0.28

clean_id <- function(x) sub("\\..*$", "", sub("\\|.*$", "", as.character(x)))

meta <- read_tsv(MANIFEST, show_col_types = FALSE)
stopifnot(all(file.exists(meta$quant_sf)))
# Salmon transcript identifiers are mapped using the generated tx2gene file.
tx2_path <- file.path(REF, "tx2gene.tsv")
if (!file.exists(tx2_path)) stop("Missing reference/tx2gene.tsv")
tx2gene <- read_tsv(tx2_path, show_col_types = FALSE)
names(tx2gene)[1:2] <- c("TXNAME", "GENEID")
tx2gene <- tx2gene |> mutate(TXNAME = clean_id(TXNAME), GENEID = clean_id(GENEID)) |> distinct()
map_gene_symbols <- function(ids) {
  ids <- unique(clean_id(ids))
  mapped <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db, keys = ids, keytype = "ENSEMBL", columns = "SYMBOL"
  )) |> as_tibble() |> transmute(gene_id = clean_id(ENSEMBL), gene_symbol = as.character(SYMBOL)) |>
    filter(!is.na(gene_id), !is.na(gene_symbol), gene_symbol != "") |> group_by(gene_id) |>
    summarise(gene_symbol = dplyr::first(gene_symbol), .groups = "drop")
  tibble(gene_id = ids) |> left_join(mapped, by = "gene_id")
}

files <- setNames(meta$quant_sf, meta$sample_id)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE,
                ignoreAfterBar = TRUE, countsFromAbundance = "no")
counts <- txi$counts

write_tsv(meta, file.path(ROOT, "metadata", "Analysis_samples.tsv"))

# Re-execute the retained reanalysis11 NONO/PSPC1 workflow. Its TMM factors
# and expression universe are established on this fixed 22-library cohort
# before either 20180718 target contrast is subset and fitted.
if (!"legacy_nono_pspc1_context" %in% names(meta)) {
  stop("Salmon_quantifications.tsv lacks legacy_nono_pspc1_context")
}
legacy_meta <- meta |>
  filter(as.logical(legacy_nono_pspc1_context)) |>
  arrange(as.integer(legacy_context_order))
if (nrow(legacy_meta) != 22L || anyDuplicated(legacy_meta$sample_id)) {
  stop("Expected 22 unique libraries in the accepted NONO/PSPC1 normalization cohort")
}
legacy_dge <- DGEList(counts[, legacy_meta$sample_id, drop = FALSE]) |>
  calcNormFactors()
legacy_filter_meta <- legacy_meta |>
  transmute(
    condition = base::factor(condition, levels = c("NT", "KD")),
    experiment = base::factor(experiment)
  )
legacy_filter_design <- model.matrix(~0 + condition + experiment, data = legacy_filter_meta)
legacy_keep <- filterByExpr(legacy_dge, design = legacy_filter_design)
legacy_dge <- legacy_dge[legacy_keep, , keep.lib.sizes = FALSE] |>
  calcNormFactors()
write_tsv(legacy_meta, file.path(ROOT, "metadata", "Normalization_NONO_PSPC1.tsv"))
write_tsv(tibble(gene_id = rownames(legacy_dge)), file.path(OUT_DATA, "DEGs_NONO_PSPC1_tested_genes.tsv"))

plot_volcano <- function(tab, factor) {
  df <- tab |> mutate(
    y = -log10(pmax(adj.P.Val, 1e-300)),
    significance = case_when(adj.P.Val <= FDR_MAX & logFC >= LFC_MIN ~ "Up",
                             adj.P.Val <= FDR_MAX & logFC <= -LFC_MIN ~ "Down",
                             TRUE ~ "NS")
  )
  p <- ggplot(df, aes(logFC, y, colour = significance)) +
    geom_point(shape = 16, size = 1.6, stroke = 0) +
    geom_vline(xintercept = c(-LFC_MIN, LFC_MIN), linetype = "dotted", linewidth = 0.7) +
    geom_hline(yintercept = -log10(FDR_MAX), linetype = "dotted", linewidth = 0.7) +
    scale_colour_manual(values = c(Down = "#3B73B9", NS = "#B3B3B3", Up = "#C73A3A"), guide = "none") +
    labs(x = paste0("Log\u2082 Fold Change (", factor, "-KD/NT)"), y = "-Log\u2081\u2080 BH-Adjusted p-value") +
    theme_classic(base_size = 13) +
    theme(axis.title.x = element_text(face = "bold", colour = "black"),
          axis.title.y = element_text(face = "bold", colour = "black"),
          axis.text = element_text(colour = "black"), plot.margin = margin(8, 10, 8, 10))
  ggsave(file.path(OUT_VIS, paste0("Volcano_", factor, ".png")), p, width = 7.2, height = 5.6, dpi = 300, bg = "white")
  ggsave(file.path(COMPAT_VOLCANO, paste0("Volcano_", factor, "_KD_vs_NT_paper.png")),
         p, width = 7.2, height = 5.6, dpi = 300, bg = "white")
}

plot_heatmap <- function(logcpm, sig, subset_meta, factor) {
  genes <- intersect(sig$gene_id, rownames(logcpm))
  matrix <- logcpm[genes, subset_meta$sample_id, drop = FALSE]
  z <- t(scale(t(matrix))); z[is.na(z)] <- 0; z[z > 2] <- 2; z[z < -2] <- -2
  labels <- paste0(ifelse(subset_meta$condition == "NT", "NT", paste0("KD_", subset_meta$shRNA)), "_", subset_meta$bioreplicate)
  colnames(z) <- labels
  annotation <- data.frame(
    Replicate = base::factor(subset_meta$bioreplicate, levels = c("L001", "L002", "L003", "L004")),
    shRNA = base::factor(ifelse(subset_meta$condition == "NT", "NT", subset_meta$shRNA)),
    Condition = base::factor(subset_meta$condition, levels = c("NT", "KD")),
    row.names = labels
  )
  ph <- pheatmap(z, color = colorRampPalette(c("#2B6CB0", "#FFFFFF", "#C53030"))(101),
                 breaks = seq(-2, 2, length.out = 101), cluster_rows = TRUE, cluster_cols = FALSE,
                 show_rownames = FALSE, border_color = NA, annotation_col = annotation,
                 annotation_colors = list(Replicate = c(L001 = "#FB8C82", L002 = "#E874D7", L003 = "#A7C800", L004 = "#D879EA"),
                                          Condition = c(NT = "#808080", KD = "#C73A3A"),
                                          shRNA = c(NT = "#808080", sh1 = "#A80F14", sh2 = "#F39C12", sh455 = "#A80F14", sh777 = "#F39C12", sh197 = "#A80F14", sh873 = "#F39C12")),
                 legend_breaks = c(-2, -1, 0, 1, 2), legend_labels = c("-2", "-1", "0", "+1", "+2"),
                 fontsize = 10, fontsize_col = 9, fontsize_annotation = 9, fontsize_legend = 9, angle_col = 45, silent = TRUE)
  # Preserve the established heatmap styling while using a portrait canvas so
  # row-level up/down patterns are readable at publication size.
  width_in <- max(8.5, grid::convertWidth(sum(ph$gtable$widths), "in", valueOnly = TRUE) + 1.0)
  height_in <- max(width_in * 1.45, grid::convertHeight(sum(ph$gtable$heights), "in", valueOnly = TRUE) + .2)
  png_type <- if (capabilities("aqua")) "quartz" else if (capabilities("cairo")) "cairo" else "Xlib"
  png(file.path(OUT_VIS, paste0("Heatmap_", factor, ".png")), width = width_in, height = height_in, units = "in", res = 300, type = png_type)
  grid.newpage(); grid.draw(ph$gtable); dev.off()
  write_tsv(tibble(gene_id = rownames(z)), file.path(OUT_DATA, paste0("Heatmap_", factor, "_genes.tsv")))
  legacy_name <- if (factor == "MCM3") {
    paste0("Heatmap_MCM3_BIOREP_DEGcount", nrow(z), "_paper_NT20161213_sh1sh2.png")
  } else {
    paste0("Heatmap_", factor, "_BIOREP_DEGcount", nrow(z), "_paper_NT20180718.png")
  }
  png(file.path(COMPAT_HEATMAP, legacy_name), width = width_in, height = height_in, units = "in", res = 300, type = png_type)
  grid.newpage(); grid.draw(ph$gtable); dev.off()
  write_tsv(tibble(gene_id = rownames(z)),
            file.path(COMPAT_HEATMAP, sub("\\.png$", "_genes.tsv", legacy_name)))
}

write_compatible_tables <- function(tab, factor) {
  sig <- tab |> filter(adj.P.Val <= FDR_MAX, abs(logFC) >= LFC_MIN)
  paths <- list(
    ALL = tab,
    SIG_adj0.05_lfc0.28 = sig,
    UP = sig |> filter(logFC > 0),
    DOWN = sig |> filter(logFC < 0)
  )
  for (suffix in names(paths)) {
    name <- if (suffix == "ALL") {
      paste0("DEG_", factor, "_KD_vs_NT_ALL.tsv")
    } else if (startsWith(suffix, "SIG")) {
      paste0("DEG_", factor, "_KD_vs_NT_", suffix, ".tsv")
    } else {
      paste0("DEG_", factor, "_", suffix, ".tsv")
    }
    write_tsv(paths[[suffix]], file.path(COMPAT_DEG, name))
  }
}

run_factor <- function(factor) {
  if (factor == "MCM3") {
    chosen <- meta |> filter(experiment == "20161213", target == "MCM3",
                             (condition == "NT" | shRNA == "sh1" | (shRNA == "sh2" & bioreplicate %in% c("L002", "L003"))))
    chosen <- chosen |> mutate(model_condition = ifelse(condition == "NT", "NT", "KD"))
    filter_kind <- "CPM >=1 in >=2 selected samples"
  } else {
    chosen <- meta |> filter(experiment == "20180718", (target == "NT" & condition == "NT") | (target == factor & condition == "KD")) |>
      mutate(model_condition = ifelse(condition == "NT", "NT", "KD"))
    filter_kind <- "accepted cohort-level edgeR filterByExpr universe"
  }
  chosen <- chosen |> mutate(display_condition_order = ifelse(condition == "NT", 0L, 1L)) |>
    arrange(display_condition_order, shRNA, bioreplicate) |> dplyr::select(-display_condition_order)
  stopifnot(nrow(chosen) > 3, all(chosen$sample_id %in% colnames(counts)))
  design <- model.matrix(~0 + base::factor(chosen$model_condition, levels = c("NT", "KD")))
  colnames(design) <- c("NT", "KD")
  if (factor == "MCM3") {
    dge <- DGEList(counts[, chosen$sample_id, drop = FALSE]) |> calcNormFactors()
    keep <- rowSums(cpm(dge) >= 1) >= 2
    dge <- dge[keep, , keep.lib.sizes = FALSE] |> calcNormFactors()
  } else {
    dge <- legacy_dge[, chosen$sample_id, keep.lib.sizes = FALSE]
  }
  v <- if (factor == "MCM3") voom(dge, design, plot = FALSE) else voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(contrasts.fit(lmFit(v, design), makeContrasts(KD_vs_NT = KD - NT, levels = design)))
  gene_map <- if (factor == "MCM3") {
    tibble(gene_id = clean_id(rownames(dge)), gene_symbol = NA_character_)
  } else {
    map_gene_symbols(rownames(dge))
  }
  tab <- topTable(fit, coef = "KD_vs_NT", number = Inf, sort.by = "P") |>
    rownames_to_column("gene_id") |>
    mutate(gene_id = clean_id(gene_id)) |>
    left_join(gene_map, by = "gene_id") |>
    mutate(target = factor, gene = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id)) |>
    dplyr::select(gene_id, logFC, AveExpr, t, P.Value, adj.P.Val, B, gene_symbol, target, gene)
  sig <- tab |> filter(adj.P.Val <= FDR_MAX, abs(logFC) >= LFC_MIN)
  write_tsv(tab, file.path(OUT_DATA, paste0("Volcano_", factor, "_data.tsv")))
  write_tsv(sig, file.path(OUT_DATA, paste0("DEGs_", factor, ".tsv")))
  write_tsv(sig |> filter(logFC > 0), file.path(OUT_DATA, paste0("VennDiagram_UP_", factor, "_genes.tsv")))
  write_tsv(sig |> filter(logFC < 0), file.path(OUT_DATA, paste0("VennDiagram_DOWN_", factor, "_genes.tsv")))
  write_tsv(chosen, file.path(ROOT, "metadata", paste0("Samples_", factor, ".tsv")))
  write_tsv(tibble(target = factor, filter = filter_kind, model = ifelse(factor == "MCM3", "limma-voom", "limma-voomWithQualityWeights"),
                   tested_genes = nrow(tab), significant_genes = nrow(sig), threshold = "BH-FDR <=0.05; |log2FC| >=0.28"),
            file.path(OUT_DATA, paste0("DEGs_", factor, "_summary.tsv")))
  write_compatible_tables(tab, factor)
  plot_volcano(tab, factor); plot_heatmap(cpm(dge, log = TRUE, prior.count = 1), sig, chosen, factor)
  list(all = tab, up = sig |> filter(logFC > 0), down = sig |> filter(logFC < 0))
}

results <- setNames(lapply(c("MCM3", "NONO", "PSPC1"), run_factor), c("MCM3", "NONO", "PSPC1"))

# Venn geometry uses an internal three-circle renderer.  The helper is part of
# this final pipeline and receives only newly generated gene sets and counts.
for (direction in c("UP", "DOWN")) {
  sets <- lapply(results, function(x) if (direction == "UP") x$up$gene_id else x$down$gene_id)
  saveRDS(sets, file.path(OUT_DATA, paste0("VennDiagram_", direction, "_sets.rds")))
  write_tsv(bind_rows(lapply(names(sets), function(name) tibble(factor = name, gene_id = sort(unique(sets[[name]]))))),
            file.path(OUT_DATA, paste0("VennDiagram_", direction, "_data.tsv")))
}
message("[DONE] RNA-seq DEG analysis: ", ROOT)
