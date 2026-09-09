#!/usr/bin/env Rscript

# Intersect promoter-bound genes with directional DEGs and render summary figures.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(data.table)
  library(grid)
})

FACTORS <- c("MCM3", "NONO", "PSPC1")
MOUSE_ENSEMBL_KEYS <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "ENSEMBL")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: render_direction.R [run_root]", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[1]
GTF_M25 <- file.path(RUN_ROOT, "reference", "gencode.vM25.annotation.gtf")
PROMOTER_DIR <- file.path(RUN_ROOT, "cutrun_work", "05_promoters")
PROMOTER_SENTINEL <- "Venn_PromoterGenes_MCM3.tsv"
if (!file.exists(file.path(PROMOTER_DIR, PROMOTER_SENTINEL))) {
  stop("Could not find final promoter-gene tables in: ", PROMOTER_DIR, call. = FALSE)
}
DEG_DIR <- file.path(RUN_ROOT, "deg_work", "data")
OUTDIR_BASE <- file.path(RUN_ROOT, "regulatory_work")
VISUAL_DIR <- file.path(OUTDIR_BASE, "visuals")
OUTDIR <- file.path(OUTDIR_BASE, "data")
VENN_HELPER <- file.path(dirname(script_file), "venn.py")
CHAIN_VENN_HELPER <- file.path(dirname(script_file), "direction_panel.py")

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

read_gene_list <- function(path) {
  stop_if_missing(path)
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"gene" %in% names(tab)) stop("Missing gene column in: ", path, call. = FALSE)
  sort(unique(trimws(as.character(tab$gene))[nzchar(trimws(as.character(tab$gene))) & !is.na(tab$gene)]))
}

write_gene_list <- function(path, genes) {
  write.table(data.frame(gene = sort(unique(genes))), path, sep = "\t", quote = FALSE, row.names = FALSE)
}

load_gencode_m25_map <- function(path) {
  stop_if_missing(path)
  gene_rows <- data.table::fread(
    cmd = sprintf("awk '$3 == \"gene\" {print $0}' %s", shQuote(path)),
    sep = "\t", header = FALSE, select = 9, col.names = "attribute", showProgress = FALSE
  )
  mapping <- data.frame(
    gene_id = sub(".*gene_id \\\"([^\\\"]+)\\\".*", "\\1", gene_rows$attribute),
    gene_symbol = sub(".*gene_name \\\"([^\\\"]+)\\\".*", "\\1", gene_rows$attribute),
    stringsAsFactors = FALSE
  )
  mapping$gene_id <- sub("\\..*$", "", mapping$gene_id)
  mapping$gene_symbol[!nzchar(mapping$gene_symbol)] <- NA_character_
  mapping[!duplicated(mapping$gene_id), , drop = FALSE]
}

GENCODE_M25_MAP <- load_gencode_m25_map(GTF_M25)

read_deg_symbols <- function(path, factor, direction, recover_missing_symbols = FALSE) {
  stop_if_missing(path)
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"gene_id" %in% names(tab)) stop("Missing gene_id column in: ", path, call. = FALSE)
  ids <- sub("\\..*$", "", trimws(as.character(tab$gene_id)))
  existing <- if ("gene_symbol" %in% names(tab)) trimws(as.character(tab$gene_symbol)) else rep(NA_character_, length(ids))
  existing[!nzchar(existing) | existing %in% c("NA", "<NA>")] <- NA_character_
  final <- existing
  gencode_symbols <- rep(NA_character_, length(ids))
  mapped <- setNames(character(), character())
  if (recover_missing_symbols) {
    gencode_symbols <- unname(GENCODE_M25_MAP$gene_symbol[match(ids, GENCODE_M25_MAP$gene_id)])
    need_gencode <- is.na(final) & !is.na(ids) & nzchar(ids)
    final[need_gencode] <- gencode_symbols[need_gencode]
    missing_ids <- unique(ids[is.na(final) & !is.na(ids) & nzchar(ids)])
    missing_ensembl_ids <- intersect(missing_ids[grepl("^ENSMUSG", missing_ids)], MOUSE_ENSEMBL_KEYS)
    mapped <- if (length(missing_ensembl_ids)) {
      suppressWarnings(AnnotationDbi::mapIds(
        org.Mm.eg.db, keys = missing_ensembl_ids, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"
      ))
    } else {
      setNames(character(), character())
    }
    need_map <- is.na(final) & !is.na(ids) & nzchar(ids)
    final[need_map] <- unname(mapped[ids[need_map]])
  }
  final[!nzchar(final) | final %in% c("NA", "<NA>")] <- NA_character_
  audit <- data.frame(
    factor = factor, direction = direction, gene_id = ids,
    original_gene_symbol = existing, final_gene_symbol = final,
    symbol_source = ifelse(!is.na(existing), "DEG_table",
                           ifelse(!is.na(gencode_symbols), "GENCODE_M25",
                                  ifelse(!is.na(final), "org.Mm.eg.db",
                                         ifelse(recover_missing_symbols, "unmapped", "historical_NA_excluded")))),
    stringsAsFactors = FALSE
  )
  list(genes = sort(unique(final[!is.na(final)])), audit = audit, n_input = length(ids))
}

assign_classes <- function(bound, up, down) {
  deg <- union(up, down)
  list(
    A_bound_UP = intersect(bound, up),
    B_bound_DOWN = intersect(bound, down),
    C_bound_NOCHANGE = setdiff(bound, deg),
    D_DEG_only_indirect = setdiff(deg, bound)
  )
}

load_panel_helpers <- function() {
  helper_env <- new.env(parent = globalenv())
  helper_env$COL_MCM3 <- "#A80F14"
  helper_env$COL_NONO <- "#F39C12"
  helper_env$COL_PSPC1 <- "#1F4AA8"
  helper_env$stack_three_pngs <- function(png_paths, out_png, out_pdf = NULL) {
    suppressPackageStartupMessages({ library(png); library(grid); library(gridExtra) })
    grobs <- lapply(png_paths, function(p) rasterGrob(readPNG(p), interpolate = TRUE))
    g <- gridExtra::arrangeGrob(grobs[[1]], grobs[[2]], grobs[[3]], ncol = 1)
    png(out_png, width = 2200, height = 5200, res = 300)
    grid.newpage(); grid.draw(g); dev.off()
    if (!is.null(out_pdf)) {
      pdf(out_pdf, width = 7.5, height = 17.0, useDingbats = FALSE)
      grid.newpage(); grid.draw(g); dev.off()
    }
  }
  helper_env
}

draw_direction_panel <- function(panel_letter, factor, n_center, n_up, n_down,
                                 n_only_center, n_only_up, n_only_down,
                                 n_overlap_up, n_overlap_down, center_col,
                                 out_png) {
  stop_if_missing(CHAIN_VENN_HELPER)
  args <- c(
    CHAIN_VENN_HELPER, "--out", out_png,
    "--panel-letter", panel_letter, "--factor", factor,
    "--n-center", n_center, "--n-up", n_up, "--n-down", n_down,
    "--n-only-center", n_only_center, "--n-only-up", n_only_up, "--n-only-down", n_only_down,
    "--n-overlap-up", n_overlap_up, "--n-overlap-down", n_overlap_down,
    "--center-col", center_col, "--col-up", "#D87070", "--col-down", "#6FA37A",
    "--col-overlap", "#BFBFBF", "--stroke-col", "#222222", "--stroke-lwd", 2,
    "--r-side", 0.170, "--r-center", 0.205,
    "--cx", 0.50, "--cy", 0.53, "--lx", 0.305, "--ly", 0.53, "--rx", 0.695, "--ry", 0.53,
    "--width-px", 3200, "--height-px", 2400, "--dpi", 300,
    "--fs-panel-letter", 20, "--fs-title", 24, "--fs-subtitle", 20,
    "--fs-num-center", 22, "--fs-num-overlap", 17,
    "--fs-bottom-lab", 20, "--fs-bottom-n", 18
  )
  status <- system2(find_python(), args = shQuote(as.character(args)))
  if (status != 0L) stop("Chain-Venn helper failed with status ", status, call. = FALSE)
}

find_python <- function() {
  candidates <- c(Sys.which("python"), Sys.which("python3"), "/opt/anaconda3/envs/cutrun_env/bin/python")
  for (candidate in candidates) {
    if (file.exists(candidate) || nzchar(Sys.which(candidate))) return(candidate)
  }
  stop("No usable Python interpreter was found.", call. = FALSE)
}

draw_venn <- function(sets, out_png, hide_numbers = FALSE) {
  stop_if_missing(VENN_HELPER)
  a <- sets$MCM3; b <- sets$NONO; c <- sets$PSPC1
  regions <- c(
    n100 = length(setdiff(a, union(b, c))),
    n010 = length(setdiff(b, union(a, c))),
    n001 = length(setdiff(c, union(a, b))),
    n110 = length(setdiff(intersect(a, b), c)),
    n101 = length(setdiff(intersect(a, c), b)),
    n011 = length(setdiff(intersect(b, c), a)),
    n111 = length(Reduce(intersect, sets))
  )
  args <- c(
    VENN_HELPER, "--out", out_png,
    "--title", "Regulatory targets overlap (promoter-bound AND DEG)",
    "--a-name", "MCM3", "--b-name", "NONO", "--c-name", "PSPC1",
    "--a-total", length(a), "--b-total", length(b), "--c-total", length(c),
    "--n100", regions[["n100"]], "--n010", regions[["n010"]], "--n001", regions[["n001"]],
    "--n110", regions[["n110"]], "--n101", regions[["n101"]], "--n011", regions[["n011"]],
    "--n111", regions[["n111"]]
  )
  if (hide_numbers) args <- c(args, "--hide-numbers")
  status <- system2(find_python(), args = shQuote(as.character(args)))
  if (status != 0L) stop("Venn helper failed with status ", status, call. = FALSE)
  regions
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VISUAL_DIR, recursive = TRUE, showWarnings = FALSE)
panel <- load_panel_helpers()
summary_rows <- list()
count_rows <- list()
direct_sets <- list()
mapping_audits <- list()
panel_pngs <- character()
panel_letters <- c("A", "B", "C")

for (i in seq_along(FACTORS)) {
  factor <- FACTORS[[i]]
  bound <- read_gene_list(file.path(PROMOTER_DIR, paste0("Venn_PromoterGenes_", factor, ".tsv")))
  up_result <- read_deg_symbols(file.path(DEG_DIR, paste0("VennDiagram_UP_", factor, "_genes.tsv")), factor, "Up",
                                recover_missing_symbols = TRUE)
  down_result <- read_deg_symbols(file.path(DEG_DIR, paste0("VennDiagram_DOWN_", factor, "_genes.tsv")), factor, "Down",
                                  recover_missing_symbols = TRUE)
  up <- up_result$genes
  down <- down_result$genes
  mapping_audits[[factor]] <- rbind(up_result$audit, down_result$audit)
  classes <- assign_classes(bound, up, down)
  direct <- union(classes$A_bound_UP, classes$B_bound_DOWN)
  direct_sets[[factor]] <- direct

  write_gene_list(file.path(OUTDIR, paste0("promoter_vs_DEG_direction_", factor, "_promoter_genes.tsv")), bound)
  write_gene_list(file.path(OUTDIR, paste0("Venn_target_", factor, "_genes.tsv")), direct)
  universe <- sort(unique(c(bound, up, down)))
  write.table(data.frame(
    gene = universe,
    promoter_bound = universe %in% bound,
    upregulated = universe %in% up,
    downregulated = universe %in% down
  ), file.path(OUTDIR, paste0("promoter_vs_DEG_direction_", factor, "_data.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE)
  for (class_name in names(classes)) {
    write_gene_list(file.path(OUTDIR, paste0("promoter_vs_DEG_direction_", factor, "_", class_name, ".tsv")), classes[[class_name]])
  }
  summary_rows[[factor]] <- data.frame(
    factor = factor, n_promoter_bound = length(bound), n_deg_up_raw = up_result$n_input, n_deg_down_raw = down_result$n_input,
    n_deg_up = length(up), n_deg_down = length(down),
    n_deg_total = length(union(up, down)), n_A_bound_UP = length(classes$A_bound_UP),
    n_B_bound_DOWN = length(classes$B_bound_DOWN), n_C_bound_NOCHANGE = length(classes$C_bound_NOCHANGE),
    n_D_deg_only_indirect = length(classes$D_DEG_only_indirect), n_direct_targets = length(direct),
    stringsAsFactors = FALSE
  )
  count_rows[[factor]] <- data.frame(factor = factor, class = names(classes), n = lengths(classes))

  center_col <- switch(factor, MCM3 = panel$COL_MCM3, NONO = panel$COL_NONO, PSPC1 = panel$COL_PSPC1)
  panel_png <- file.path(VISUAL_DIR, paste0("promoter_vs_DEG_direction_", factor, ".png"))
  draw_direction_panel(
    panel_letter = panel_letters[[i]], factor = factor, n_center = length(bound), n_up = length(up), n_down = length(down),
    n_only_center = length(classes$C_bound_NOCHANGE), n_only_up = length(up) - length(classes$A_bound_UP),
    n_only_down = length(down) - length(classes$B_bound_DOWN), n_overlap_up = length(classes$A_bound_UP),
    n_overlap_down = length(classes$B_bound_DOWN), center_col = center_col,
    out_png = panel_png
  )
  panel_pngs <- c(panel_pngs, panel_png)
}

shared <- Reduce(intersect, direct_sets)
regions <- draw_venn(
  direct_sets,
  file.path(VISUAL_DIR, "venn_regulatory_targets_MCM3_NONO_PSPC1.png")
)
draw_venn(
  direct_sets,
  file.path(VISUAL_DIR, "venn_regulatory_targets_MCM3_NONO_PSPC1_no_numbers.png"),
  hide_numbers = TRUE
)
write.table(do.call(rbind, summary_rows), file.path(OUTDIR, "promoter_vs_DEG_direction_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, count_rows), file.path(OUTDIR, "promoter_vs_DEG_direction_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, mapping_audits), file.path(OUTDIR, "GeneSymbolMappingAudit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_gene_list(file.path(OUTDIR, "Venn_target_shared_genes.tsv"), shared)
write.table(data.frame(region = names(regions), n = as.integer(regions)),
            file.path(OUTDIR, "Venn_target_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

region_genes <- list(
  MCM3_only = setdiff(direct_sets$MCM3, union(direct_sets$NONO, direct_sets$PSPC1)),
  NONO_only = setdiff(direct_sets$NONO, union(direct_sets$MCM3, direct_sets$PSPC1)),
  PSPC1_only = setdiff(direct_sets$PSPC1, union(direct_sets$MCM3, direct_sets$NONO)),
  MCM3_NONO_only = setdiff(intersect(direct_sets$MCM3, direct_sets$NONO), direct_sets$PSPC1),
  MCM3_PSPC1_only = setdiff(intersect(direct_sets$MCM3, direct_sets$PSPC1), direct_sets$NONO),
  NONO_PSPC1_only = setdiff(intersect(direct_sets$NONO, direct_sets$PSPC1), direct_sets$MCM3),
  MCM3_NONO_PSPC1 = shared
)
write.table(do.call(rbind, lapply(names(region_genes), function(region) {
  data.frame(region = region, gene = sort(region_genes[[region]]), stringsAsFactors = FALSE)
})), file.path(OUTDIR, "Venn_target_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

panel$stack_three_pngs(panel_pngs,
  out_png = file.path(VISUAL_DIR, "Fig_promoter_vs_DEG_direction.png")
)

writeLines(c(
  "direct_target_rule\tpromoter-bound AND DEG",
  "CUTRUN_promoter_source\tcutrun_work/05_promoters/Venn_PromoterGenes_<FACTOR>.tsv",
  "CUTRUN_threshold\tPooled six-library support >=4/6; 2019 p <=1e-4, q <=0.05, FE >=3; 2020 pooled matched-IgG p <=1e-4, q <=0.01, FE >=3; TSS +/-1,000 bp; promoter overlap >=250 bp",
  "analysis_role\tPooled 2019/2020 promoter evidence; display coverage uses batch-balanced mean tracks",
  "RNAseq_source\tfinal-local DEG tables",
  "RNAseq_threshold\tBH FDR <= 0.05 and absolute log2FC >= 0.28",
  "gene_symbol_mapping\tMissing symbols mapped from Ensembl IDs through GENCODE M25 then org.Mm.eg.db"
), file.path(OUTDIR, "AnalysisParameters.tsv"))

message("[DONE] Wrote promoter-versus-DEG and regulatory-target outputs to: ", OUTDIR)
