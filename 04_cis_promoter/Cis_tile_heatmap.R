###############################################################################
#  Cis_tile_heatmap.R
#  --------------------------------------------------------------------------
#  Single R script that produces both the QrLBD and QsLBD cis-regulatory
#  element tile heatmaps (Figure 6a and 6b) in a journal-friendly layout
#  matching the user's curated chapter figure: numbers-in-cells + heat
#  colour fill, with a separate "5-category totals" panel on the right.
#
#  Layout per species:
#     LEFT panel   : 22 curated CREs grouped under 5 functional category
#                    banners (Hormone / Stress / Light/circadian /
#                    Promoter-Core / T.F. Binding) -- NO MERISTEM.
#     RIGHT panel  : per-gene totals collapsed to the 5 categories.
#     Cell content : raw count (no z-score, no prefilter), printed in
#                    each cell, plus a sequential white -> orange -> red
#                    fill that scales with the count.
#
#  Input files (sit next to this script):
#     Cis_curated_QrLBD42.xlsx   sheets: 'detail' (42 x 22) + 'category_totals' (42 x 5)
#     Cis_curated_QsLBD37.xlsx   sheets: 'detail' (37 x 22) + 'category_totals' (37 x 5)
#
#  Outputs written to this folder:
#     Cis_QrLBD_tile.png  /  .pdf
#     Cis_QsLBD_tile.png  /  .pdf
#
#  Run from the cispromo/ folder:
#       Rscript Cis_tile_heatmap.R
#  or in RStudio:  setwd() to this folder, then source().
#
#  Required packages: ggplot2, readxl, dplyr, tidyr, scales, patchwork
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

WORK_DIR <- "."

# ---------------------------------------------------------------------------
#  Read the 'detail' sheet (gene x 22 motifs with category banner in row 1)
# ---------------------------------------------------------------------------
read_detail <- function(path) {
  body  <- read_excel(path, sheet = "detail", skip = 2, col_names = FALSE)
  ncols <- ncol(body)
  cats  <- as.character(unlist(read_excel(path, sheet = "detail", n_max = 1,
                                          col_names = FALSE)))
  mots  <- as.character(unlist(read_excel(path, sheet = "detail",
                                          range = cell_rows(2),
                                          col_names = FALSE)))
  if (length(cats) < ncols) cats <- c(cats, rep(NA, ncols - length(cats)))
  if (length(mots) < ncols) mots <- c(mots, rep(NA, ncols - length(mots)))
  for (i in seq_along(cats)) {
    if (is.na(cats[i]) || cats[i] == "") {
      cats[i] <- if (i > 1) cats[i - 1] else NA
    }
  }
  colnames(body) <- mots
  names(body)[1] <- "Gene"
  body <- body[!is.na(body$Gene) & body$Gene != "", , drop = FALSE]

  long <- body %>%
    tidyr::pivot_longer(-Gene, names_to = "Motif", values_to = "Count") %>%
    dplyr::mutate(Count = as.numeric(Count))
  cat_lookup <- data.frame(Motif = mots[-1], Category = cats[-1],
                           stringsAsFactors = FALSE)
  long <- dplyr::left_join(long, cat_lookup, by = "Motif")
  long$Gene  <- factor(long$Gene,  levels = body$Gene)             # phylogeny order
  long$Motif <- factor(long$Motif, levels = mots[-1])
  long$Category <- factor(long$Category, levels = unique(cats[-1]))
  long
}

# ---------------------------------------------------------------------------
#  Read the 'category_totals' sheet (gene x 5 categories)
# ---------------------------------------------------------------------------
read_totals <- function(path) {
  d <- read_excel(path, sheet = "category_totals")
  long <- d %>%
    tidyr::pivot_longer(-Gene, names_to = "Motif", values_to = "Count") %>%
    dplyr::mutate(Count = as.numeric(Count))
  long$Gene  <- factor(long$Gene,  levels = d$Gene)
  long$Motif <- factor(long$Motif, levels = setdiff(colnames(d), "Gene"))
  long$Category <- long$Motif        # placeholder so make_tile() works
  long
}

# ---------------------------------------------------------------------------
#  Tile heatmap with cell numbers + heat-colour fill
# ---------------------------------------------------------------------------
make_tile <- function(df, value_col = "Count",
                      facet = TRUE, max_global = NULL,
                      label_fontsize = 2.6,
                      cell_ratio = 1) {
  if (is.null(max_global)) max_global <- max(df[[value_col]], na.rm = TRUE)
  txt_thresh <- max_global * 0.55     # white text above this fill intensity
  df$txt_col <- ifelse(df[[value_col]] >= txt_thresh, "white", "grey20")
  df$lab     <- ifelse(df[[value_col]] == 0, "0",
                       formatC(df[[value_col]], format = "d"))

  p <- ggplot(df, aes(x = Motif, y = Gene)) +
    geom_tile(aes(fill = Count), colour = "grey45", linewidth = 0.30) +
    geom_text(aes(label = lab, colour = txt_col),
              size = label_fontsize, family = "sans", fontface = "bold") +
    scale_colour_identity() +
    scale_fill_gradientn(
      colours = c("white", "#FFF3CC", "#FED976",
                  "#FB8A4F", "#E6371F", "#A50026"),
      name   = "Count",
      limits = c(0, max_global),
      oob    = scales::squish,
      breaks = pretty(c(0, max_global), n = 5)
    ) +
    scale_y_discrete(limits = rev(levels(df$Gene))) +    # phylogeny top-down
    labs(x = NULL, y = NULL) +
    # Cell aspect controlled via figure width in configs (coord_fixed is
    # incompatible with facet_grid(space="free_x")).
    theme_minimal(base_size = 9) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                      size = 7, colour = "grey15"),
      axis.text.y      = element_text(size = 7, colour = "grey15"),
      strip.text       = element_text(face = "bold", size = 7.5,
                                      colour = "grey10",
                                      lineheight = 0.9),
      strip.background = element_rect(fill = "grey88", colour = "grey60",
                                      linewidth = 0.3),
      panel.spacing    = unit(0.10, "lines"),
      legend.position  = "bottom",
      legend.key.width = unit(20, "pt"),
      legend.key.height= unit(7,  "pt"),
      legend.title     = element_text(face = "bold")
    )
  if (facet) {
    # label_wrap_gen wraps long category banners across two lines so the
    # short panels (e.g. "Promoter-Core related") don't get truncated.
    p <- p + facet_grid(. ~ Category,
                        scales = "free_x", space = "free_x",
                        labeller = labeller(Category = label_wrap_gen(width = 14)))
  }
  p
}

# ---------------------------------------------------------------------------
#  Compose detail (left) + 5-category totals (right) for one species
# ---------------------------------------------------------------------------
make_panel <- function(xlsx, species_label) {
  cat("\n==========================================================\n")
  cat("  Tile heatmap:", species_label, "\n")
  cat("==========================================================\n")

  detail <- read_detail(xlsx)
  totals <- read_totals(xlsx)

  cat("  detail   :", nlevels(detail$Gene), "x", nlevels(detail$Motif), "\n")
  cat("  totals   :", nlevels(totals$Gene), "x", nlevels(totals$Category), "\n")

  # Shared colour scale max = max of either matrix (so the two panels are
  # on the same fill scale, like the user's figure)
  global_max <- max(max(detail$Count, na.rm = TRUE),
                    max(totals$Count, na.rm = TRUE), na.rm = TRUE)

  p_left  <- make_tile(detail, max_global = global_max, label_fontsize = 2.3) +
             theme(axis.text.y = element_blank())     # gene names only on right
  p_right <- make_tile(totals, facet = FALSE, max_global = global_max,
                       label_fontsize = 2.6) +
             theme(legend.position = "none",
                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  # Left wider (~22 cols) than right (~5 cols)
  combined <- p_left + p_right + plot_layout(widths = c(22, 5)) +
              plot_annotation(title = species_label,
                              theme = theme(plot.title =
                                element_text(face = "bold", hjust = 0.5,
                                             size = 11)))
  combined
}

# ---------------------------------------------------------------------------
#  Run for both species
# ---------------------------------------------------------------------------
configs <- list(
  Qr = list(xlsx = "Cis_curated_QrLBD42.xlsx",
            label = "(A)  Q. rubra LBD promoter cis-elements (n = 42)",
            short = "(A)  Q. rubra (n = 42)",
            out_png = "Cis_QrLBD_tile.png", out_pdf = "Cis_QrLBD_tile.pdf",
            width = 9, height = 8),
  Qs = list(xlsx = "Cis_curated_QsLBD37.xlsx",
            label = "(B)  Q. suber LBD promoter cis-elements (n = 37)",
            short = "(B)  Q. suber (n = 37)",
            out_png = "Cis_QsLBD_tile.png", out_pdf = "Cis_QsLBD_tile.pdf",
            width = 9, height = 7)
)

# ---------------------------------------------------------------------------
#  PASS 1 - Supplementary: full per-motif detail per species
# ---------------------------------------------------------------------------
cat("\n*****  PASS 1 - Detailed (supplementary) figures  *****\n")
for (sp in names(configs)) {
  cfg <- configs[[sp]]
  if (!file.exists(cfg$xlsx)) {
    warning("Missing input: ", cfg$xlsx, " (skipping)"); next
  }
  fig <- make_panel(cfg$xlsx, cfg$label)
  ggsave(cfg$out_png, fig, width = cfg$width, height = cfg$height,
         dpi = 300, bg = "white")
  ggsave(cfg$out_pdf, fig, width = cfg$width, height = cfg$height,
         bg = "white")
  cat("  wrote:", cfg$out_png, "\n")
}

# ---------------------------------------------------------------------------
#  PASS 2 - Main figure: 5-category summary, both species side-by-side,
#  shared colour scale so Qr and Qs are directly comparable.
# ---------------------------------------------------------------------------
cat("\n*****  PASS 2 - Compact main-text summary figure  *****\n")

qr_tot <- read_totals(configs$Qr$xlsx)
qs_tot <- read_totals(configs$Qs$xlsx)
shared_max <- max(qr_tot$Count, qs_tot$Count, na.rm = TRUE)
cat("  shared colour-scale max =", shared_max, "\n")

main_left  <- make_tile(qr_tot, facet = FALSE, max_global = shared_max,
                        label_fontsize = 3.4) +
              labs(subtitle = configs$Qr$short) +
              theme(legend.position = "none",
                    plot.subtitle  = element_text(face = "bold", size = 11,
                                                  hjust = 0.5))
main_right <- make_tile(qs_tot, facet = FALSE, max_global = shared_max,
                        label_fontsize = 3.4) +
              labs(subtitle = configs$Qs$short) +
              theme(plot.subtitle  = element_text(face = "bold", size = 11,
                                                  hjust = 0.5))

main_fig <- main_left + main_right + plot_layout(widths = c(1, 1))
ggsave("Cis_main_summary.png", main_fig, width = 11, height = 9,
       dpi = 300, bg = "white")
ggsave("Cis_main_summary.pdf", main_fig, width = 11, height = 9,
       bg = "white")
cat("  wrote:  Cis_main_summary.png\n")
cat("  wrote:  Cis_main_summary.pdf\n")

cat("\n==========================================================\n")
cat("  Done.\n")
cat("    Main figure (use in paper)        :  Cis_main_summary.png/.pdf\n")
cat("    Supplementary detail (Qr / Qs)    :  Cis_QrLBD_tile.png /\n")
cat("                                          Cis_QsLBD_tile.png\n")
cat("==========================================================\n")
