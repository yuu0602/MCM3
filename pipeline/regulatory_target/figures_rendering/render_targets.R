#!/usr/bin/env Rscript

# Render regulatory-target pies and gene-body profiles.

options(stringsAsFactors = FALSE)
options(scipen = 999)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[1]
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
if (is.na(script_file)) stop("Cannot resolve pipeline script path", call. = FALSE)
DEFAULT_RUN_ROOT <- normalizePath(file.path(dirname(script_file), "..", "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: render_targets.R [run_root]", call. = FALSE)
RUN_ROOT <- normalizePath(if (length(args) >= 1L) args[[1]] else DEFAULT_RUN_ROOT, mustWork = TRUE)
COVERAGE_AXIS_LABEL <- "Coverage"
GROUP_BASE_DIR <- file.path(RUN_ROOT, "regulatory_work")
REGTARGET_BASE_DIR <- GROUP_BASE_DIR
GROUP_VISUAL_DIR <- file.path(GROUP_BASE_DIR, "visuals")
REGTARGET_VISUAL_DIR <- file.path(REGTARGET_BASE_DIR, "visuals")
GROUP_DIR <- file.path(GROUP_BASE_DIR, "data")
REGTARGET_DIR <- file.path(REGTARGET_BASE_DIR, "data")
GTF <- file.path(RUN_ROOT, "reference", "gencode.vM25.annotation.gtf")
COMPUTE_MATRIX <- Sys.which("computeMatrix")
PLOT_PROFILE <- Sys.which("plotProfile")
if (!nzchar(COMPUTE_MATRIX)) COMPUTE_MATRIX <- "/opt/anaconda3/envs/cutrun_env/bin/computeMatrix"
if (!nzchar(PLOT_PROFILE)) PLOT_PROFILE <- "/opt/anaconda3/envs/cutrun_env/bin/plotProfile"
BIGWIGS <- c(
  MCM3 = file.path(RUN_ROOT, "cutrun_work", "03_bigwig", "IGV_representation", "MCM3_mean.bw"),
  NONO = file.path(RUN_ROOT, "cutrun_work", "03_bigwig", "IGV_representation", "NONO_mean.bw"),
  PSPC1 = file.path(RUN_ROOT, "cutrun_work", "03_bigwig", "IGV_representation", "PSPC1_mean.bw")
)

GROUPS <- list(
  A = list(label = "MCM3 + NONO + PSPC1", region = "MCM3_NONO_PSPC1", tracks = c("MCM3", "NONO", "PSPC1"), title = "MCM3-PSPC1-NONO direct targets (n=%s)"),
  B = list(label = "MCM3 + PSPC1 only", region = "MCM3_PSPC1_only", tracks = c("MCM3", "PSPC1"), title = "MCM3-PSPC1-only direct targets (n=%s)"),
  C = list(label = "NONO + PSPC1 only", region = "NONO_PSPC1_only", tracks = c("NONO", "PSPC1"), title = "NONO-PSPC1-only direct targets (n=%s)"),
  D = list(label = "MCM3 + NONO only", region = "MCM3_NONO_only", tracks = c("MCM3", "NONO"), title = "MCM3-NONO-only direct targets (n=%s)")
)
COLORS <- c(MCM3 = "#cf111f", NONO = "#f2a010", PSPC1 = "#4f7db1")
DIR_COLORS <- c(Up = "#D87070", Down = "#6FA37A", Mix = "#F3D37A")
BIN_SIZE <- 25L
BEFORE_BP <- 3000L
BODY_BP <- 3000L
AFTER_BP <- 3000L
DISPLAY_SCALE <- 1

for (pkg in c("data.table", "ggplot2", "grid", "png", "rtracklayer", "GenomicRanges", "GenomeInfoDb", "S4Vectors")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg, call. = FALSE)
}
dir.create(GROUP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REGTARGET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GROUP_VISUAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REGTARGET_VISUAL_DIR, recursive = TRUE, showWarnings = FALSE)
for (path in c(GROUP_DIR, GTF, COMPUTE_MATRIX, PLOT_PROFILE, unname(BIGWIGS))) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}
Sys.setenv(MPLCONFIGDIR = file.path(tempdir(), "mcm3_final_mplconfig"))

read_genes <- function(path) {
  tab <- data.table::fread(path, data.table = FALSE)
  if (!"gene" %in% names(tab)) stop("Missing gene column: ", path, call. = FALSE)
  sort(unique(trimws(as.character(tab$gene))[nzchar(trimws(as.character(tab$gene))) & !is.na(tab$gene)]))
}

write_genes <- function(path, genes) {
  data.table::fwrite(data.frame(gene = sort(unique(genes))), path, sep = "\t", quote = FALSE)
}

class_genes <- function(factor, class_name) {
  read_genes(file.path(GROUP_DIR, sprintf("promoter_vs_DEG_direction_%s_%s.tsv", factor, class_name)))
}

direction <- function(gene, factor) {
  if (gene %in% class_sets[[factor]]$up) return("Up")
  if (gene %in% class_sets[[factor]]$down) return("Down")
  stop("Direct target has no assigned direction: ", gene, " / ", factor, call. = FALSE)
}

summarize_group <- function(genes, factors) {
  rows <- lapply(sort(unique(genes)), function(gene) {
    dirs <- vapply(factors, function(factor) direction(gene, factor), character(1))
    cls <- if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mix"
    out <- data.frame(gene = gene, class = cls, stringsAsFactors = FALSE)
    for (factor in factors) out[[factor]] <- dirs[[factor]]
    out
  })
  if (!length(rows)) return(data.frame(gene = character(), class = character(), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

pie_plot <- function(title, tab, show_numbers = TRUE, show_labels = TRUE) {
  dat <- as.data.frame(table(factor(tab$class, levels = c("Up", "Down", "Mix"))), stringsAsFactors = FALSE)
  names(dat) <- c("class", "n")
  dat <- dat[dat$n > 0, , drop = FALSE]
  dat$frac <- dat$n / sum(dat$n)
  dat$label <- if (!show_labels) "" else if (show_numbers) paste0(dat$class, "\n(n=", dat$n, ")") else as.character(dat$class)
  ggplot2::ggplot(dat, ggplot2::aes(x = "", y = frac, fill = class)) +
    ggplot2::geom_col(width = 1, color = "#222222", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 4, fontface = "bold") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = DIR_COLORS) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8))
}

binding_pie <- function(factor, all_targets, direct_targets, out, show_numbers = TRUE) {
  dat <- data.frame(
    class = c("Targets w/ promoter binding", "Targets w/o promoter binding"),
    n = c(length(direct_targets), length(all_targets) - length(direct_targets)),
    stringsAsFactors = FALSE
  )
  dat$frac <- dat$n / sum(dat$n)
  dat$label <- if (show_numbers) paste0(sprintf("%.1f", 100 * dat$frac), "%\n(n=", dat$n, ")") else ""
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = "", y = frac, fill = class)) +
    ggplot2::geom_col(width = 1, color = "#222222", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 4, fontface = "bold") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = c("Targets w/ promoter binding" = "#BFD9F2", "Targets w/o promoter binding" = "#9AA3C7")) +
    ggplot2::labs(
      title = if (show_numbers) sprintf("%s targets (n=%s)", factor, format(length(all_targets), big.mark = ",")) else sprintf("%s targets", factor),
      subtitle = "Promoter binding from CUT&RUN (TSS +/- 1 kb)"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 11), plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5), plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5))
  ggplot2::ggsave(out, p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
}

read_gtf_genes <- function(path) {
  gr <- rtracklayer::import(path)
  gr <- gr[S4Vectors::mcols(gr)$type == "gene"]
  data.frame(
    chrom = as.character(GenomeInfoDb::seqnames(gr)), start = as.integer(GenomicRanges::start(gr)) - 1L,
    end = as.integer(GenomicRanges::end(gr)), strand = as.character(GenomicRanges::strand(gr)),
    gene = as.character(S4Vectors::mcols(gr)$gene_name), stringsAsFactors = FALSE
  )
}

write_gene_body_bed <- function(gtf_genes, genes, path) {
  dat <- gtf_genes[gtf_genes$gene %in% genes, , drop = FALSE]
  missing <- setdiff(genes, unique(dat$gene))
  if (length(missing)) stop("Genes absent from GTF: ", paste(missing, collapse = ", "), call. = FALSE)
  dat <- dat[dat$end > dat$start & dat$strand %in% c("+", "-"), c("chrom", "start", "end", "gene", "strand")]
  dat$score <- 0L
  dat <- dat[, c("chrom", "start", "end", "gene", "score", "strand")]
  data.table::fwrite(dat, path, sep = "\t", quote = FALSE, col.names = FALSE)
  nrow(dat)
}

run_command <- function(exe, args) {
  status <- system2(exe, args = vapply(as.character(args), shQuote, character(1)))
  if (!identical(status, 0L)) stop("Command failed: ", exe, call. = FALSE)
}

read_profile <- function(path, sample_order) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  lines <- lines[!(tolower(trimws(lines)) %in% c("bin labels", "bins"))]
  rows <- list()
  for (line in lines) {
    parts <- strsplit(trimws(line), "[[:space:]]+", perl = TRUE)[[1]]
    if (length(parts) < 3 || tolower(parts[[1]]) %in% c("bin", "bins")) next
    nums <- suppressWarnings(as.numeric(parts[-c(1, 2)]))
    nums <- nums[is.finite(nums)]
    if (!length(nums)) next
    rows[[length(rows) + 1L]] <- list(sample = parts[[1]], region = parts[[2]], values = nums)
  }
  if (!length(rows)) stop("Cannot parse plotProfile table: ", path, call. = FALSE)
  max_bins <- max(vapply(rows, function(row) length(row$values), integer(1)))
  parsed <- do.call(rbind, lapply(rows, function(row) {
    values <- c(row$values, rep(NA_real_, max_bins - length(row$values)))
    data.frame(sample = row$sample, region = row$region, t(values), check.names = FALSE)
  }))
  bin_cols <- names(parsed)[-(1:2)]
  out <- data.table::melt(data.table::as.data.table(parsed), id.vars = c("sample", "region"), measure.vars = bin_cols, variable.name = "bin", value.name = "signal")
  out[, bin_index := match(bin, bin_cols)]
  # Display the full -3 kb / scaled gene body / +3 kb span on a -3 to +3 axis.
  # TSS and TES are consequently internal landmarks at -1 and +1, respectively.
  out[, pos_kb := -3 + (bin_index - 1L) * (6 / (length(bin_cols) - 1L))]
  as.data.frame(out[, .(signal = mean(signal, na.rm = TRUE)), by = .(sample, region, pos_kb)])
}

make_group_profile <- function(group_name, genes, tracks, gtf_genes) {
  prefix <- file.path(REGTARGET_DIR, sprintf("Metaprofile_Group%s", group_name))
  bed <- paste0(prefix, "_regions.bed")
  n_bed <- write_gene_body_bed(gtf_genes, genes, bed)
  matrix <- paste0(prefix, ".matrix.gz")
  profile <- paste0(prefix, "_data.tsv")
  raw_plot <- tempfile(pattern = paste0("metaprofile_", group_name, "_"), fileext = ".png")
  on.exit(unlink(raw_plot), add = TRUE)
  run_command(COMPUTE_MATRIX, c("scale-regions", "-S", unname(BIGWIGS[tracks]), "-R", bed, "--beforeRegionStartLength", BEFORE_BP, "--regionBodyLength", BODY_BP, "--afterRegionStartLength", AFTER_BP, "--binSize", BIN_SIZE, "--missingDataAsZero", "--samplesLabel", tracks, "-o", matrix))
  run_command(PLOT_PROFILE, c("-m", matrix, "--plotType", "lines", "--legendLocation", "upper-right", "--plotTitle", sprintf(GROUPS[[group_name]]$title, n_bed), "--plotFileFormat", "png", "--plotHeight", "6", "--plotWidth", "9", "--outFileNameData", profile, "-out", raw_plot))
  long <- read_profile(profile, tracks)
  long$signal <- long$signal * DISPLAY_SCALE
  long$sample <- factor(long$sample, levels = tracks)
  yy <- long$signal[is.finite(long$signal)]
  ymax <- max(yy, na.rm = TRUE) * 1.05
  ymin <- 0
  p <- ggplot2::ggplot(long, ggplot2::aes(pos_kb, signal, color = sample)) +
    ggplot2::geom_line(linewidth = 1.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 1, color = "#222222") +
    ggplot2::scale_color_manual(values = COLORS[tracks], breaks = tracks, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -1, 1, 3), labels = c("-3", "TSS", "TES", "3")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax)) +
    ggplot2::labs(title = sprintf(GROUPS[[group_name]]$title, n_bed), x = "Relative distance (kb)", y = COVERAGE_AXIS_LABEL) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5), axis.line = ggplot2::element_line(linewidth = 1.0), panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0), axis.ticks = ggplot2::element_line(linewidth = 1.0), axis.text = ggplot2::element_text(size = 14, face = "bold", color = "black"), axis.title = ggplot2::element_text(size = 16, face = "bold"), legend.position = "right", legend.text = ggplot2::element_text(size = 14, face = "bold"), plot.margin = ggplot2::margin(10, 12, 12, 12))
  out_png <- file.path(REGTARGET_VISUAL_DIR, paste0("metaprofile_geneBody_", group_name, ".single.paper.png"))
  # Match the single-panel Profile_MCM3_NONO aspect ratio (7.4:5.6).
  ggplot2::ggsave(out_png, p, width = 10.5, height = 10.5 * 5.6 / 7.4, dpi = 300, bg = "white")
  for (path in c(out_png, profile, bed)) {
    target <- file.path(if (grepl("\\.png$", path)) GROUP_VISUAL_DIR else GROUP_DIR, basename(path))
    if (!identical(normalizePath(path, mustWork = FALSE), normalizePath(target, mustWork = FALSE))) file.copy(path, target, overwrite = TRUE)
  }
  data.frame(group = group_name, definition = GROUPS[[group_name]]$label, genes = length(genes), bed_rows = n_bed, tracks = paste(tracks, collapse = ","), stringsAsFactors = FALSE)
}

make_mcm3_target_profile <- function(gtf_genes) {
  up <- class_sets$MCM3$up
  down <- class_sets$MCM3$down
  up_bed <- file.path(REGTARGET_DIR, "Metaprofile_MCM3_target_UP_regions.bed")
  down_bed <- file.path(REGTARGET_DIR, "Metaprofile_MCM3_target_DOWN_regions.bed")
  n_up <- write_gene_body_bed(gtf_genes, up, up_bed)
  n_down <- write_gene_body_bed(gtf_genes, down, down_bed)
  prefix <- file.path(REGTARGET_DIR, "Metaprofile_MCM3_target")
  matrix <- paste0(prefix, ".matrix.gz")
  profile <- paste0(prefix, "_data.tsv")
  raw_plot <- tempfile(pattern = "PanelC_metaprofile_", fileext = ".png")
  on.exit(unlink(raw_plot), add = TRUE)
  # This deepTools build has no computeMatrix --regionsLabel option. The two
  # region files retain their order and are relabelled explicitly below.
  run_command(COMPUTE_MATRIX, c("scale-regions", "-S", BIGWIGS[["MCM3"]], "-R", up_bed, down_bed, "-b", BEFORE_BP, "-a", AFTER_BP, "--binSize", BIN_SIZE, "--regionBodyLength", BODY_BP, "--missingDataAsZero", "--startLabel", "TSS", "--endLabel", "TES", "--samplesLabel", "MCM3", "-o", matrix))
  run_command(PLOT_PROFILE, c("-m", matrix, "--plotType", "lines", "--legendLocation", "upper-right", "--plotTitle", "MCM3 CUT&RUN across gene bodies (TSS to TES; scaled)", "--plotFileFormat", "png", "--plotHeight", "6", "--plotWidth", "9", "--outFileNameData", profile, "-out", raw_plot))
  long <- read_profile(profile, "MCM3")
  long$signal <- long$signal * DISPLAY_SCALE
  region_levels <- unique(long$region)
  if (length(region_levels) != 2L) stop("Expected two MCM3 target region classes in Panel C profile.", call. = FALSE)
  long$direction <- ifelse(long$region == region_levels[[1]], "MCM3 up-regulated targets", "MCM3 down-regulated targets")
  yy <- long$signal[is.finite(long$signal)]
  ymax <- max(yy, na.rm = TRUE) * 1.05
  ymin <- 0
  p <- ggplot2::ggplot(long, ggplot2::aes(pos_kb, signal, color = direction)) +
    ggplot2::geom_rect(
      data = data.frame(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#F2F2F2", color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.7, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.9, color = "#444444") +
    ggplot2::scale_color_manual(
      values = c("MCM3 up-regulated targets" = "#2855B2", "MCM3 down-regulated targets" = "#74E37D"),
      breaks = c("MCM3 up-regulated targets", "MCM3 down-regulated targets"), name = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -1, 1, 3), labels = c("-3.0", "TSS", "TES", "3.0")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax)) +
    ggplot2::labs(title = "MCM3 CUT&RUN across gene bodies (TSS\u2192TES; scaled)", x = "Relative distance (kb)", y = COVERAGE_AXIS_LABEL) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 21, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 10)), axis.line = ggplot2::element_line(linewidth = 1), panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0), axis.ticks = ggplot2::element_line(linewidth = 1), axis.text = ggplot2::element_text(size = 15, face = "bold", color = "black"), axis.title = ggplot2::element_text(size = 18, face = "bold"), legend.position = "top", legend.direction = "horizontal", legend.text = ggplot2::element_text(size = 15, face = "bold"), legend.key.width = grid::unit(1.1, "cm"), plot.margin = ggplot2::margin(18, 18, 18, 24))
  out_png <- file.path(REGTARGET_VISUAL_DIR, "PanelC_metaprofile_MCM3_targets_geneBody_TSS_TES_scaled_pub.png")
  ggplot2::ggsave(out_png, p, width = 14.6667, height = 14.6667 * 5.6 / 7.4, dpi = 300, bg = "white")
  target <- file.path(GROUP_VISUAL_DIR, basename(out_png))
  if (!identical(normalizePath(out_png, mustWork = FALSE), normalizePath(target, mustWork = FALSE))) file.copy(out_png, target, overwrite = TRUE)
  write.table(data.frame(direction = c("Up", "Down"), n_genes = c(n_up, n_down)), paste0(prefix, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

direct_sets <- setNames(lapply(c("MCM3", "NONO", "PSPC1"), function(factor) read_genes(file.path(GROUP_DIR, paste0("Venn_target_", factor, "_genes.tsv")))), c("MCM3", "NONO", "PSPC1"))
class_sets <- setNames(lapply(c("MCM3", "NONO", "PSPC1"), function(factor) list(up = class_genes(factor, "A_bound_UP"), down = class_genes(factor, "B_bound_DOWN"), indirect = class_genes(factor, "D_DEG_only_indirect"))), c("MCM3", "NONO", "PSPC1"))

groups <- list(
  A = Reduce(intersect, direct_sets),
  B = setdiff(intersect(direct_sets$MCM3, direct_sets$PSPC1), direct_sets$NONO),
  C = setdiff(intersect(direct_sets$NONO, direct_sets$PSPC1), direct_sets$MCM3),
  D = setdiff(intersect(direct_sets$MCM3, direct_sets$NONO), direct_sets$PSPC1)
)
regions <- data.table::fread(file.path(GROUP_DIR, "Venn_target_genes.tsv"), data.table = FALSE)
summary_rows <- list()
pie_paths <- character()
for (name in names(GROUPS)) {
  expected <- sort(unique(regions$gene[regions$region == GROUPS[[name]]$region]))
  if (!setequal(groups[[name]], expected)) stop("Venn mismatch for group ", name, call. = FALSE)
  tab <- summarize_group(groups[[name]], GROUPS[[name]]$tracks)
  write_genes(file.path(GROUP_DIR, paste0("Group", name, "_genes.tsv")), groups[[name]])
  write.table(tab, file.path(GROUP_DIR, paste0("Group", name, "_directions.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  counts <- table(factor(tab$class, levels = c("Up", "Down", "Mix")))
  summary_rows[[name]] <- data.frame(group = name, group_definition = GROUPS[[name]]$label, venn_region = GROUPS[[name]]$region, n_genes = nrow(tab), n_up = counts[["Up"]], n_down = counts[["Down"]], n_mix = counts[["Mix"]], stringsAsFactors = FALSE)
  plot <- pie_plot(sprintf("Group %s: %s\nn = %s", name, GROUPS[[name]]$label, nrow(tab)), tab)
  out <- file.path(GROUP_VISUAL_DIR, paste0("group_", name, "_pie.png"))
  ggplot2::ggsave(out, plot, width = 5.2, height = 5.2, dpi = 300, bg = "white")
  pie_paths <- c(pie_paths, out)
  plot_no_labels <- pie_plot(NULL, tab, show_numbers = FALSE, show_labels = FALSE)
  ggplot2::ggsave(file.path(GROUP_VISUAL_DIR, paste0("group_", name, "_pie_no_labels.png")), plot_no_labels, width = 5.2, height = 5.2, dpi = 300, bg = "white")
}
write.table(do.call(rbind, summary_rows), file.path(GROUP_DIR, "Pie_Groups_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(group = names(GROUPS), group_definition = vapply(GROUPS, `[[`, character(1), "label"), venn_region = vapply(GROUPS, `[[`, character(1), "region"), n_genes = vapply(groups, length, integer(1))), file.path(GROUP_DIR, "Pie_Groups_definitions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

grDevices::png(file.path(GROUP_VISUAL_DIR, "groups_A_B_C_D_pies.png"), width = 3150, height = 3300, res = 300, bg = "white")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 2, heights = grid::unit(c(0.6, 5.2, 5.2), "in"), widths = grid::unit(c(5.25, 5.25), "in"))))
grid::grid.text("Regulatory-target overlap groups (Up / Down / Mix)", vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = grid::gpar(fontsize = 18, fontface = "bold"))
for (i in seq_along(pie_paths)) {
  image <- png::readPNG(pie_paths[[i]])
  grid::pushViewport(grid::viewport(layout.pos.row = if (i <= 2) 2 else 3, layout.pos.col = if (i %% 2 == 1) 1 else 2))
  grid::grid.raster(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
  grid::popViewport()
}
grid::popViewport()
grDevices::dev.off()

input_parameters <- readLines(file.path(GROUP_DIR, "AnalysisParameters.tsv"), warn = FALSE)
writeLines(c(input_parameters, "mix_color\t#F3D37A (yellow)"), file.path(GROUP_DIR, "Pie_Groups_parameters.tsv"))

for (factor in c("NONO", "PSPC1")) {
  all_targets <- union(union(class_sets[[factor]]$up, class_sets[[factor]]$down), class_sets[[factor]]$indirect)
  binding_pie(factor, all_targets, direct_sets[[factor]], file.path(GROUP_VISUAL_DIR, paste0("pie_targets_promoter_binding_", factor, ".png")))
  binding_pie(factor, all_targets, direct_sets[[factor]], file.path(GROUP_VISUAL_DIR, paste0("pie_targets_promoter_binding_", factor, "_no_numbers.png")), show_numbers = FALSE)
}

gtf_genes <- read_gtf_genes(GTF)
profile_summary <- do.call(rbind, lapply(names(GROUPS), function(name) make_group_profile(name, groups[[name]], GROUPS[[name]]$tracks, gtf_genes)))
write.table(profile_summary, file.path(GROUP_DIR, "Metaprofile_Groups_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
make_mcm3_target_profile(gtf_genes)
message("[DONE] Regulatory-target figures written to: ", GROUP_DIR)
