#!/usr/bin/env Rscript
## ==================================================================
## 05_figures.R
##
## An R alternative to 05_figures.py, reading the same output files and
## producing the same six figures. Use whichever you prefer; they are
## not meant to be run together.
##
## Packages (install once):
##   install.packages(c("ggplot2","patchwork","ggdist","scico","dplyr",
##                      "readr","tidyr","forcats","ggtext","scales"))
##
##   ggplot2   grammar of graphics
##   patchwork multipanel assembly with automatic panel tagging, which is
##             far less fiddly than gridExtra or cowplot for this
##   ggdist    slab/interval geoms; shows a whole distribution plus its
##             summary in one layer, better than boxplot-plus-jitter
##   scico     perceptually uniform, colourblind-safe palettes (Crameri)
##   ggtext    markdown in titles, so a coloured group name can be inlined
##
## Usage:
##   Rscript 05_figures.R --focal 12
##   Rscript 05_figures.R --focal 12 --isout is_out --meta metadata.tsv \
##       --out figures_R --only fig6
##
## Figures: fig1 IS load, fig2 recency, fig3 site frequency spectrum,
## fig3c relative spectrum (focal vs other groups), fig4 recombination,
## fig5 accessory gene distance, fig6 shell genome IS frequency,
## fig7 insertion spectrum against the neutral synonymous reference.
## The supplementary assembly-balance figure is produced separately by
## 06_assembly_qc_by_group.py.
## ==================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(forcats); library(patchwork); library(scales)
})
has_ggdist <- requireNamespace("ggdist", quietly = TRUE)
## interval_size_range was renamed across ggdist versions; check once and
## fall back to defaults rather than failing inside a layer
gd_args <- list()
if (has_ggdist) {
  fm <- names(formals(ggdist::stat_pointinterval))
  if ("interval_size_range" %in% fm ||
      "..." %in% fm) gd_args <- list(interval_size_range = c(0.5, 1.4))
}
has_scico  <- requireNamespace("scico",  quietly = TRUE)
has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
if (!has_ggdist) message("note: ggdist not installed; using boxplots instead")

## ---------- arguments ---------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
ISOUT <- getarg("--isout", "is_out")
META  <- getarg("--meta",  "metadata.tsv")
OUT   <- getarg("--out",   "figures_R")
FOCAL <- trimws(as.character(getarg("--focal", NA)))
ONLY  <- if ("--only" %in% args) args[(match("--only", args) + 1):length(args)] else NULL
if (is.na(FOCAL)) stop("--focal is required, e.g. --focal 12")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---------- shared design -----------------------------------------
## One accent for the focal group, neutral grey for everything else.
## Colour encodes the comparison and nothing decorative.
C_FOCAL   <- "#B0353C"
C_FOCAL_L <- "#E8B4B7"
C_NEUT    <- "#7A7A7A"
C_NEUT_L  <- "#D6D6D6"
INK       <- "#1A1A1A"
MUTED     <- "#6E6E6E"
PAL <- c(focal = C_FOCAL, other = C_NEUT)

theme_pub <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.3),
      axis.line          = element_line(colour = "#4D4D4D", linewidth = 0.3),
      axis.ticks         = element_line(colour = "#4D4D4D", linewidth = 0.3),
      axis.text          = element_text(colour = MUTED, size = base_size - 1.5),
      axis.title         = element_text(colour = INK),
      plot.title         = element_text(colour = INK, size = base_size + 0.5,
                                        hjust = 0, margin = margin(b = 4)),
      plot.title.position = "plot",
      legend.position    = "top",
      legend.justification = "left",
      legend.key.size    = unit(8, "pt"),
      legend.margin      = margin(0, 0, 0, 0),
      legend.title       = element_blank(),
      plot.margin        = margin(4, 6, 4, 4)
    )
}
theme_set(theme_pub())

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         dpi = 400, bg = "white")
  message("  wrote ", name, ".pdf / .png")
}

want <- function(key) is.null(ONLY) || key %in% ONLY
exists_all <- function(...) all(file.exists(c(...)))

## ---------- fig1: IS load ------------------------------------------
fig1 <- function() {
  calls_p <- file.path(ISOUT, "IS_calls.tsv")
  if (!exists_all(calls_p, META)) { message("  skip fig1: missing inputs"); return() }
  calls <- read_tsv(calls_p, show_col_types = FALSE, progress = FALSE)
  meta  <- read_tsv(META, col_types = cols(fastbaps_group = col_character()),
                    progress = FALSE) |>
    mutate(genome = as.character(genome),
           fastbaps_group = trimws(fastbaps_group))
  load <- calls |> mutate(genome = as.character(genome)) |>
    count(genome, name = "n_IS")
  d <- meta |> left_join(load, by = "genome") |>
    mutate(n_IS = coalesce(n_IS, 0L),
           y = if ("total_len" %in% names(meta)) n_IS / total_len * 1e6 else n_IS,
           grp = fct_reorder(fastbaps_group, y, .fun = median, .desc = TRUE),
           arm = if_else(fastbaps_group == FOCAL, "focal", "other"))
  ylab <- if ("total_len" %in% names(meta)) "IS copies per Mb" else "IS copies per genome"
  ns <- d |> count(grp)

  ## ggdist's half-eye shows the density AND the interval, which conveys
  ## more than a boxplot without the clutter of overlaid points
  base <- ggplot(d, aes(grp, y, colour = arm, fill = arm))
  pA <- if (has_ggdist) {
    ## linewidth must be length 1 or nrow(data); to get thick 50% and thin
    ## 95% intervals use interval_size_range, which ggdist maps internally
    base +
      ggdist::stat_slab(aes(fill = arm), alpha = 0.35,
                        normalize = "groups", show.legend = FALSE) +
      do.call(ggdist::stat_pointinterval,
              c(list(.width = c(0.5, 0.95), point_size = 1.6,
                     show.legend = FALSE), gd_args))
  } else {
    base + geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.4,
                        show.legend = FALSE)
  }
  pA <- pA +
    geom_text(data = ns, aes(grp, Inf, label = n), inherit.aes = FALSE,
              vjust = 1.6, size = 2.1, colour = MUTED) +
    scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
    labs(x = "fastbaps group", y = ylab,
         title = "IS content across genetic groups") +
    expand_limits(y = 0)

  fe <- file.path(ISOUT, "02_focal_effect_summary.tsv")
  if (file.exists(fe)) {
    f <- read_tsv(fe, show_col_types = FALSE, progress = FALSE) |>
      mutate(model = fct_rev(fct_inorder(sub(":.*", "", model))))
    pB <- ggplot(f, aes(fold_change, model)) +
      geom_vline(xintercept = 1, linetype = "22", colour = MUTED,
                 linewidth = 0.3) +
      geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), colour = "#4D4D4D",
                     linewidth = 0.6) +
      geom_point(size = 2.2, colour = C_FOCAL) +
      labs(x = "fold change, focal vs other", y = NULL,
           title = "Effect across models")
    p <- pA + pB + plot_layout(widths = c(1.75, 1)) +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 11))
    save_fig(p, "fig1_IS_load", 7.4, 3.2)
  } else {
    save_fig(pA, "fig1_IS_load", 4.8, 3.2)
  }
}

## ---------- fig2: recency ------------------------------------------
fig2 <- function() {
  p0 <- file.path(ISOUT, "IS_recency_per_genome.tsv")
  if (!file.exists(p0)) { message("  skip fig2: missing ", p0); return() }
  r <- read_tsv(p0, show_col_types = FALSE, progress = FALSE)
  if (!"qfocal" %in% names(r)) { message("  skip fig2: no qfocal column"); return() }
  r <- r |> mutate(arm = as.character(qfocal)) |> filter(!is.na(p_recent))
  lab <- r |> count(arm) |> mutate(lab = paste0(arm, " (n=", n, ")"))
  r <- left_join(r, lab, by = "arm")

  ## histogram rather than a density: p_recent is a bounded proportion
  ## computed from a modest number of pairs per genome, so it is lumpy at
  ## simple fractions (1/2, 1/3, 2/3). A KDE smooths that structure away
  ## and invents support outside [0, 1].
  pA <- ggplot(r, aes(p_recent, fill = arm, colour = arm)) +
    geom_histogram(aes(y = after_stat(density)), bins = 33, alpha = 0.35,
                   position = "identity", linewidth = 0.4) +
    scale_fill_manual(values = PAL, labels = lab$lab, breaks = lab$arm) +
    scale_colour_manual(values = PAL, labels = lab$lab, breaks = lab$arm) +
    ## coord_cartesian zooms; scale_x_continuous(limits=) DROPS bars whose
    ## edges fall outside, which is what produced the "Removed 4 rows"
    ## warning and silently deleted the extreme bins
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "within-genome IS pairs >99% identical (proportion)",
         y = "density", title = "Recency of IS copies")

  pB <- ggplot(r, aes(arm, p_recent, colour = arm, fill = arm))
  pB <- if (has_ggdist) {
    pB + ggdist::stat_halfeye(adjust = 1, width = 0.5, .width = 0,
                              justification = -0.25, alpha = 0.35,
                              show.legend = FALSE) +
      do.call(ggdist::stat_pointinterval,
              c(list(.width = c(0.5, 0.95), point_size = 1.8,
                     show.legend = FALSE), gd_args)) +
      geom_jitter(width = 0.06, size = 0.35, alpha = 0.2,
                  show.legend = FALSE)
  } else {
    pB + geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.3,
                      show.legend = FALSE)
  }
  pB <- pB + scale_colour_manual(values = PAL) +
    scale_fill_manual(values = PAL) +
    labs(x = NULL, y = "proportion >99% identical",
         title = "Per-genome summary")

  p <- pA + pB + plot_layout(widths = c(1.5, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  save_fig(p, "fig2_recency", 7.2, 3.1)
}

## ---------- fig3: site frequency spectrum ---------------------------
fig3 <- function() {
  p0 <- file.path(ISOUT, "insertion_frequencies.tsv")
  if (!file.exists(p0)) { message("  skip fig3: missing ", p0); return() }
  ## python writes the site index as an unnamed first column, which
  ## readr renames to ...1 with a message; name it explicitly instead
  f <- read_tsv(p0, show_col_types = FALSE, progress = FALSE,
                name_repair = function(x) { x[x == ""] <- "site"; x })
  if (!nrow(f)) { message("  skip fig3: empty"); return() }
  if (!"freq" %in% names(f)) { message("  skip fig3: no freq column"); return() }

  pA <- ggplot(f, aes(freq)) +
    geom_histogram(bins = 25, fill = C_FOCAL, colour = "white",
                   linewidth = 0.25) +
    labs(x = "insertion frequency", y = "insertion sites",
         title = paste0("Site frequency spectrum (n = ", nrow(f), ")"))

  hi <- mean(f$freq >= 0.5)
  pB <- ggplot(f, aes(freq)) +
    stat_ecdf(geom = "step", colour = C_FOCAL, linewidth = 0.9) +
    geom_vline(xintercept = 0.5, linetype = "22", colour = MUTED,
               linewidth = 0.3) +
    annotate("text", x = 0.52, y = 0.25,
             label = sprintf("%.1f%% of sites\nat frequency >= 0.5", 100 * hi),
             hjust = 0, size = 2.6, colour = INK) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "insertion frequency", y = "cumulative fraction of sites",
         title = "Cumulative distribution")

  p <- pA + pB + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  save_fig(p, "fig3_sfs", 7.2, 3.1)
}


## ---------- fig3c: relative SFS, focal vs other groups ---------------
## Absolute counts cannot be compared between groups: they differ in both
## the number of segregating sites and the number of genomes. Recomputing
## frequency WITHIN each group and plotting each group's own proportion
## per bin puts the two spectra on a common scale.
fig3c <- function() {
  mp <- file.path(ISOUT, "insertion_matrix.tsv")
  if (!exists_all(mp, META)) {
    message("  skip fig3c: needs insertion_matrix.tsv and metadata"); return()
  }
  mat <- suppressMessages(read_tsv(mp, show_col_types = FALSE,
                                   progress = FALSE))
  names(mat)[1] <- "site"
  meta <- read_tsv(META, col_types = cols(fastbaps_group = col_character()),
                   progress = FALSE) |>
    mutate(genome = as.character(genome),
           fastbaps_group = trimws(fastbaps_group))
  gmap <- setNames(meta$fastbaps_group, meta$genome)

  gcols <- setdiff(names(mat), "site")
  focal_cols <- gcols[which(gmap[gcols] == FOCAL)]
  other_cols <- gcols[which(!is.na(gmap[gcols]) & gmap[gcols] != FOCAL)]
  if (length(focal_cols) < 5 || length(other_cols) < 5) {
    message("  skip fig3c: too few genomes matched to groups"); return()
  }

  ## a site scored in too few genomes gives an unreliable frequency; and a
  ## site absent from the group entirely is not segregating there
  spectrum <- function(cols, min_frac = 0.5) {
    sub <- as.matrix(mat[, cols, drop = FALSE])
    scored <- rowSums(!is.na(sub))
    keep <- scored >= min_frac * length(cols)
    fr <- rowMeans(sub, na.rm = TRUE)[keep]
    fr <- fr[is.finite(fr) & fr > 0]
    fr
  }
  fF <- spectrum(focal_cols); fO <- spectrum(other_cols)
  if (length(fF) < 20 || length(fO) < 20) {
    message("  skip fig3c: too few scorable sites per group"); return()
  }

  brks <- seq(0, 1, length.out = 21)
  ctr  <- head(brks, -1) + diff(brks) / 2
  binit <- function(x) {
    h <- hist(x, breaks = brks, plot = FALSE)$counts
    h / sum(h)
  }
  pF <- binit(fF); pO <- binit(fO)

  labF <- sprintf("group %s (%d sites, %d genomes)", FOCAL, length(fF),
                  length(focal_cols))
  labO <- sprintf("other groups (%d sites, %d genomes)", length(fO),
                  length(other_cols))
  dd <- tibble(mid = rep(ctr, 2),
               prop = c(pO, pF),
               set  = factor(rep(c(labO, labF), each = length(ctr)),
                             levels = c(labO, labF)))

  pA <- ggplot(dd, aes(mid, prop, fill = set)) +
    geom_col(position = position_dodge(width = 0.045), width = 0.04) +
    scale_fill_manual(values = setNames(c(C_NEUT, C_FOCAL), c(labO, labF))) +
    labs(x = "insertion frequency within group",
         y = "proportion of that group's sites",
         title = "Relative site frequency spectrum") +
    theme(legend.position = "top")

  rt <- tibble(mid = ctr, ratio = ifelse(pO > 0, pF / pO, NA_real_)) |>
    filter(is.finite(ratio), ratio > 0)
  pB <- ggplot(rt, aes(mid, ratio)) +
    geom_hline(yintercept = 1, linetype = "22", colour = MUTED,
               linewidth = 0.3) +
    geom_line(colour = C_FOCAL, linewidth = 0.8) +
    geom_point(colour = C_FOCAL, size = 1.6) +
    scale_y_log10() +
    labs(x = "insertion frequency within group",
         y = paste0("group ", FOCAL, " / other groups"),
         title = "Enrichment across the spectrum",
         caption = "above 1: relatively more sites at that frequency") +
    theme(plot.caption = element_text(hjust = 0, size = 6.5, colour = MUTED))

  p <- pA + pB + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  save_fig(p, "fig3c_relative_sfs", 7.4, 3.2)
}

## ---------- fig4: recombination -------------------------------------
fig4 <- function() {
  dens <- file.path(ISOUT, "recombination_density_by_group.tsv")
  corr <- file.path(ISOUT, "IS_vs_recombination_by_group.tsv")
  parts <- list()

  if (file.exists(dens)) {
    D <- read_tsv(dens, show_col_types = FALSE, progress = FALSE)
    pos <- D[[1]]
    grp_cols <- setdiff(names(D), names(D)[1])
    others <- setdiff(grp_cols, FOCAL)
    dd <- tibble(pos = pos,
                 other = if (length(others))
                   rowMeans(D[, others, drop = FALSE], na.rm = TRUE) else NA_real_,
                 focal = if (FOCAL %in% grp_cols) D[[FOCAL]] else NA_real_)
    parts$A <- ggplot(dd, aes(pos / 1e6)) +
      geom_area(aes(y = other), fill = C_NEUT_L, colour = NA) +
      geom_line(aes(y = other, colour = "other groups (mean)"),
                linewidth = 0.4) +
      geom_line(aes(y = focal, colour = paste("group", FOCAL)),
                linewidth = 0.6) +
      scale_colour_manual(values = setNames(c(C_NEUT, C_FOCAL),
                                            c("other groups (mean)",
                                              paste("group", FOCAL)))) +
      labs(x = "core-genome position (Mb)",
           y = "recombination density per genome",
           title = "Recombination along the core")
  }

  if (file.exists(corr)) {
    C <- read_tsv(corr, show_col_types = FALSE, progress = FALSE) |>
      mutate(group = as.character(group),
             arm = if_else(group == FOCAL, "focal", "other"),
             group = fct_reorder(group, spearman_rho))
    parts$B <- ggplot(C, aes(spearman_rho, group, colour = arm)) +
      geom_vline(xintercept = 0, linetype = "22", colour = MUTED,
                 linewidth = 0.3) +
      geom_segment(aes(x = 0, xend = spearman_rho, yend = group),
                   linewidth = 0.5, alpha = 0.5, show.legend = FALSE) +
      geom_point(size = 2.1, show.legend = FALSE) +
      scale_colour_manual(values = PAL) +
      labs(x = expression(paste("Spearman ", rho, ", IS load vs recombination")),
           y = "fastbaps group", title = "Within-group association")
  }

  if (!length(parts)) { message("  skip fig4: no recombination outputs"); return() }
  p <- Reduce(`+`, parts) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  save_fig(p, "fig4_recombination", if (length(parts) > 1) 7.6 else 4.2, 3.1)
}

## ---------- fig5: accessory gene distance ----------------------------
fig5 <- function() {
  p0 <- file.path(ISOUT, "gene_distance_to_IS.tsv")
  if (!file.exists(p0)) { message("  skip fig5: missing ", p0); return() }
  raw <- read_tsv(p0, show_col_types = FALSE, progress = FALSE) |>
    filter(is.finite(dist_IS))
  ## One row per gene PER GENOME. Plotting those directly counts the same
  ## cluster once per genome carrying it, which inflates n roughly 300-fold
  ## and makes the curves look far more precise than the data are.
  ## Collapse to one value per cluster first.
  d <- if ("cluster" %in% names(raw)) {
    raw |> group_by(cluster) |>
      summarise(dist_IS = median(dist_IS),
                selected = any(as.logical(selected)),
                n_genomes = n(), .groups = "drop")
  } else raw
  d <- d |> mutate(set = if_else(as.logical(selected), "under selection",
                                 "not selected"))
  lab <- d |> count(set) |>
    mutate(lab = paste0(set, " (n=", n, " genes)"))
  meds <- d |> group_by(set) |> summarise(m = median(dist_IS), .groups = "drop")

  p <- ggplot(d, aes(dist_IS / 1000, colour = set)) +
    stat_ecdf(geom = "step", linewidth = 0.9) +
    geom_vline(data = meds, aes(xintercept = m / 1000, colour = set),
               linetype = "22", linewidth = 0.35, show.legend = FALSE) +
    ## pseudo_log keeps 0 visible while compressing the long tail; clamp
    ## the lower limit so the transform does not draw its negative half
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = 0.4),
                       breaks = c(0, 1, 10, 100, 1000),
                       labels = c("0", "1", "10", "100", "1000"),
                       limits = c(0, NA)) +
    scale_colour_manual(values = setNames(c(C_FOCAL, C_NEUT),
                                          c("under selection",
                                            "not selected")),
                        labels = lab$lab, breaks = lab$set) +
    labs(x = "median distance to nearest IS element (kb)",
         y = "cumulative fraction of accessory gene clusters",
         title = "Proximity of accessory genes to IS elements",
         ## fig5 and fig6 run in OPPOSITE directions with respect to IS
         ## association; without saying so they look contradictory
         caption = "\u2190 closer to IS                         further from IS \u2192") +
    theme(plot.caption = element_text(hjust = 0.5, size = 6.5,
                                      colour = MUTED))
  save_fig(p, "fig5_selection", 4.8, 3.2)
}

## ---------- fig6: shell genome IS frequency ---------------------------
fig6 <- function() {
  p0 <- file.path(ISOUT, "shell_cluster_IS_association.tsv")
  if (!file.exists(p0)) { message("  skip fig6: missing ", p0); return() }
  d <- read_tsv(p0, show_col_types = FALSE, progress = FALSE) |>
    mutate(set = if_else(as.logical(selected), "under selection",
                         "other shell genes"))
  if (min(table(d$set)) < 3) { message("  skip fig6: too few clusters"); return() }
  cols <- setNames(c(C_FOCAL, C_NEUT),
                   c("under selection", "other shell genes"))
  lab <- d |> count(set) |> mutate(lab = paste0(set, " (n=", n, ")"))

  pA <- ggplot(d, aes(prox_frac, colour = set)) +
    stat_ecdf(geom = "step", linewidth = 0.9) +
    scale_colour_manual(values = cols, labels = lab$lab, breaks = lab$set) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "carrying genomes with an IS nearby (fraction)",
         y = "cumulative fraction of shell genes",
         title = "IS proximity of shell genes",
         caption = "\u2190 less IS-associated          more IS-associated \u2192") +
    theme(plot.caption = element_text(hjust = 0.5, size = 6.5,
                                      colour = MUTED))

  ## separate panel rather than an inset: an inset over the ECDF hides
  ## exactly the region where the curves separate
  pB <- ggplot(d, aes(set, prox_frac, colour = set, fill = set))
  pB <- if (has_ggdist) {
    pB + ggdist::stat_halfeye(width = 0.5, .width = 0, justification = -0.2,
                              alpha = 0.35, show.legend = FALSE) +
      do.call(ggdist::stat_pointinterval,
              c(list(.width = c(0.5, 0.95), point_size = 1.8,
                     show.legend = FALSE), gd_args)) +
      geom_jitter(width = 0.06, size = 0.35, alpha = 0.18, show.legend = FALSE)
  } else {
    pB + geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.3,
                      show.legend = FALSE)
  }
  pB <- pB + scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_x_discrete(labels = c("other", "selected")) +
    labs(x = NULL, y = "proximity fraction", title = "Distribution")

  pC <- NULL
  if ("freq" %in% names(d) && dplyr::n_distinct(d$freq) > 4) {
    dd <- d |> mutate(band = cut(freq, breaks = 4, include.lowest = TRUE)) |>
      group_by(band, set) |>
      summarise(m = mean(prox_frac), se = sd(prox_frac) / sqrt(n()),
                n = n(), .groups = "drop")
    pC <- ggplot(dd, aes(band, m, colour = set, group = set)) +
      geom_line(position = position_dodge(0.25), linewidth = 0.5) +
      geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se),
                    width = 0.12, position = position_dodge(0.25),
                    linewidth = 0.4) +
      geom_point(position = position_dodge(0.25), size = 2) +
      ## counts on two rows, one per series, so the dodged labels cannot
      ## run into each other
      geom_text(aes(y = Inf, label = n,
                    vjust = if_else(set == "under selection", 1.6, 3.2)),
                size = 1.9, show.legend = FALSE) +
      scale_colour_manual(values = cols) +
      expand_limits(y = 0) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
      labs(x = "gene prevalence band", y = "mean proximity fraction",
           title = "Stratified by prevalence") +
      theme(axis.text.x = element_text(size = 6))
  }

  p <- if (is.null(pC)) {
    pA + pB + plot_layout(widths = c(1.4, 0.8))
  } else {
    pA + pB + pC + plot_layout(widths = c(1.35, 0.7, 1.35))
  }
  p <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  save_fig(p, "fig6_shell_IS_frequency", if (is.null(pC)) 7.0 else 9.6, 3.2)
}


## ---------- fig7: insertion spectrum vs neutral reference -------------
fig7 <- function() {
  p0 <- file.path(ISOUT, "sfs_projected.tsv")
  if (!file.exists(p0)) {
    message("  skip fig7: missing ", p0, " (run 07_neutral_sfs_comparison.py)")
    return()
  }
  d <- read_tsv(p0, show_col_types = FALSE, progress = FALSE)
  need_cols <- c("minor_allele_count", "synonymous_prop", "insertions_prop")
  if (!all(need_cols %in% names(d))) {
    message("  skip fig7: unexpected columns"); return()
  }

  long <- d |>
    select(k = minor_allele_count, synonymous = synonymous_prop,
           insertions = insertions_prop) |>
    pivot_longer(-k, names_to = "set", values_to = "prop") |>
    mutate(set = factor(set, levels = c("synonymous", "insertions"),
                        labels = c("synonymous SNPs", "IS insertions")))
  cols <- setNames(c(C_NEUT, C_FOCAL), c("synonymous SNPs", "IS insertions"))

  pA <- ggplot(long, aes(k, prop, fill = set)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = cols) +
    scale_y_log10() +
    labs(x = "minor allele count (projected)", y = "proportion of sites",
         title = "Folded site frequency spectra")

  rt <- d |>
    mutate(ratio = if_else(synonymous_prop > 0,
                           insertions_prop / synonymous_prop, NA_real_)) |>
    filter(is.finite(ratio), ratio > 0)
  pB <- ggplot(rt, aes(minor_allele_count, ratio)) +
    geom_ribbon(aes(ymin = 1, ymax = pmax(ratio, 1)), fill = C_FOCAL_L,
                alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "22", colour = MUTED,
               linewidth = 0.3) +
    geom_line(colour = C_FOCAL, linewidth = 0.8) +
    geom_point(colour = C_FOCAL, size = 1.5) +
    scale_y_log10() +
    labs(x = "minor allele count (projected)",
         y = "insertions / synonymous",
         title = "Departure from neutrality",
         caption = "above 1: excess of insertions at that frequency") +
    theme(plot.caption = element_text(hjust = 0, size = 6.5, colour = MUTED))

  cp <- file.path(ISOUT, "sfs_masking_comparison.tsv")
  if (file.exists(cp)) {
    cc <- read_tsv(cp, show_col_types = FALSE, progress = FALSE)
    pC <- ggplot(cc, aes(enrichment, reference)) +
      geom_vline(xintercept = 1, linetype = "22", colour = MUTED,
                 linewidth = 0.3) +
      geom_segment(aes(x = 1, xend = enrichment, yend = reference),
                   colour = C_NEUT, linewidth = 0.5) +
      geom_point(colour = C_FOCAL, size = 2.4) +
      geom_text(aes(label = sprintf("%.1fx (%d sites)", enrichment,
                                    n_syn_sites)),
                hjust = -0.15, size = 2.1, colour = MUTED) +
      expand_limits(x = max(cc$enrichment) * 1.7) +
      labs(x = "high-frequency enrichment", y = NULL,
           title = "Neutral reference")
    p <- pA + pB + pC + plot_layout(widths = c(1.4, 1.2, 0.8))
    wd <- 10.2
  } else {
    p <- pA + pB + plot_layout(widths = c(1.3, 1.1))
    wd <- 7.4
  }
  p <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  save_fig(p, "fig7_neutral_sfs", wd, 3.2)
}

## ---------- run ------------------------------------------------------
message("writing figures to ", OUT, "/")
for (nm in c("fig1", "fig2", "fig3", "fig3c", "fig4", "fig5",
              "fig6", "fig7")) {
  if (!want(nm)) next
  res <- try(get(nm)(), silent = TRUE)
  if (inherits(res, "try-error"))
    message("  ", nm, " failed: ", conditionMessage(attr(res, "condition")))
}
message("done")
