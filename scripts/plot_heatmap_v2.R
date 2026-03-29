# plot_heatmap.R
# Per-nucleotide mutation enrichment heatmap using ggplot2
#
# Reads wtnorm_fc.py TSV output.
# Produces the same layout as plot_heatmap.py:
#   heatmap (4 bases x positions, color = mean_log2fc)
#   WT base row below (colored tiles with base letter + position number)
#
# Usage:
#   Rscript plot_heatmap.R \
#     --tsv    scnv_wtnorm_fc.tsv \
#     --title  "scnv XR2" \
#     --out    scnv_heatmap.pdf \
#     --vmax   4
#
# Requirements:
#   install.packages(c("ggplot2", "patchwork", "optparse", "dplyr"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(optparse)
})

# ── CLI -----------------------------------------------------------------------
option_list <- list(
  make_option("--tsv",       type="character", help="wtnorm_fc TSV"),
  make_option("--title",     type="character", default="Mutation enrichment heatmap"),
  make_option("--out",       type="character", default="heatmap.pdf"),
  make_option("--vmax",      type="double",    default=4.0,
              help="Color scale ±max [default 4.0]"),
  make_option("--region_start", type="integer", default=NULL,
              help="Start position (1-based, optional)"),
  make_option("--region_end",   type="integer", default=NULL,
              help="End position (1-based, optional)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$tsv)) stop("--tsv required")

BASES      <- c("A","C","G","T")
BASE_COLORS <- c(A="#E64B35", C="#4DBBD5", G="#3C9A3E", T="#3C5488")

# ── load data -----------------------------------------------------------------
df <- read.table(opt$tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# exclude WT rows
df <- df[df$mut_base != df$wt_base, ]

# optional region filter
if (!is.null(opt$region_start)) df <- df[df$pos >= opt$region_start, ]
if (!is.null(opt$region_end))   df <- df[df$pos <= opt$region_end,   ]

# clamp color scale
df$l2fc_clamped <- pmax(pmin(df$mean_log2fc, opt$vmax), -opt$vmax)

# factor ordering
df$mut_base <- factor(df$mut_base, levels=rev(BASES))  # A top → T bottom

positions <- sort(unique(df$pos))

# WT base per position
wt_df <- df %>%
  group_by(pos) %>%
  summarise(wt_base = first(wt_base), .groups="drop")

# ── heatmap -------------------------------------------------------------------
p_heat <- ggplot(df, aes(x=factor(pos), y=mut_base, fill=l2fc_clamped)) +
  geom_tile(color="white", linewidth=0.3) +
  # WT cell outline
  geom_tile(
    data = df %>%
      inner_join(wt_df, by="pos") %>%
      filter(mut_base == wt_base.y),
    aes(x=factor(pos), y=mut_base),
    fill=NA, color="black", linewidth=0.9
  ) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-opt$vmax, opt$vmax),
    name     = paste0("mean log₂(WT-norm FC)\n±", opt$vmax, " scale")
  ) +
  scale_x_discrete(breaks=factor(positions[positions %% 5 == 0]),
                   labels=positions[positions %% 5 == 0]) +
  labs(x=NULL, y="mutation\nto base") +
  theme_minimal(base_size=10) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(family="mono", size=9, face="bold"),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.title     = element_text(size=8),
    legend.text      = element_text(size=7),
    plot.title       = element_text(face="bold", size=12, hjust=0)
  ) +
  ggtitle(opt$title)

# ── WT base row ---------------------------------------------------------------
p_wt <- ggplot(wt_df, aes(x=factor(pos), y=1, fill=wt_base)) +
  geom_tile(color="white", linewidth=0.3) +
  geom_text(aes(label=wt_base), size=2.2, fontface="bold",
            color="white", vjust=0.8) +
  # position number every 5th
  geom_text(
    data = wt_df %>% filter(pos %% 5 == 0 | pos == min(pos)),
    aes(label=pos), size=1.8, color="white", vjust=2.4
  ) +
  scale_fill_manual(values=BASE_COLORS, guide="none") +
  scale_x_discrete(breaks=factor(positions[positions %% 5 == 0]),
                   labels=positions[positions %% 5 == 0]) +
  labs(x="position", y="WT") +
  theme_minimal(base_size=9) +
  theme(
    axis.text.x      = element_text(size=7, angle=0, hjust=0.5),
    axis.text.y      = element_text(size=8, face="bold"),
    axis.ticks.x     = element_line(linewidth=0.3),
    panel.grid       = element_blank(),
    plot.margin      = margin(0, 5, 2, 5)
  )

# ── combine with patchwork ----------------------------------------------------
n_pos  <- length(positions)
width  <- max(8, n_pos * 0.22)

combined <- p_heat / p_wt +
  plot_layout(heights=c(4, 0.5))

ggsave(opt$out, combined,
       width=width, height=5.5, units="in",
       device=if (grepl("\\.pdf$", opt$out)) "pdf" else "png",
       dpi=200)

cat("Saved:", opt$out, "\n")
