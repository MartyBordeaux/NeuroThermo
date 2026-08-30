suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})

COL_WT <- "#2B6CB0"
COL_SCA3 <- "#C44E52"
COL_EXP <- "#222222"
COL_MODEL <- "#777777"
COL_EARLY <- "#3C8D73"
COL_COUPLED <- "#4C5F9E"
COL_LATE <- "#B36A3C"
COL_KAPPA <- "#7A5195"
COL_J <- "#EF5675"
COL_COMBINED <- "#003F5C"
COL_INTERACTION <- "#A0A0A0"

stage_cols <- c("WT exit"="#4D4D4D", "Balance"="#8C6D31", "SCA3 entry"="#8B1A1A")
stage_lty  <- c("WT exit"="dotted", "Balance"="dashed", "SCA3 entry"="dotdash")
path_cols <- c("Drive early"=COL_EARLY, "Coupled"=COL_COUPLED, "Drive late"=COL_LATE)
group_cols <- c("WT"=COL_WT, "SCA3"=COL_SCA3)

qfun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(data.frame(y=NA_real_, ymin=NA_real_, ymax=NA_real_))
  data.frame(y=median(x), ymin=unname(quantile(x, .25)), ymax=unname(quantile(x, .75)))
}

theme_nt <- function(base_size=9) {
  theme_classic(base_size=base_size) +
    theme(
      axis.title = element_text(face="plain"),
      axis.text = element_text(color="black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face="bold"),
      legend.title = element_blank(),
      legend.key.height = grid::unit(0.35, "cm"),
      plot.margin = margin(5.5, 7, 5.5, 5.5)
    )
}

legend_inset_theme <- function(x=.03, y=.03, size=9.5, justification=c(0, 0)) {
  theme(
    legend.position = "inside",
    legend.position.inside = c(x, y),
    legend.justification = justification,
    legend.direction = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.88), color = "#D0D0D0", linewidth = 0.25),
    legend.key.height = grid::unit(0.48, "cm"),
    legend.key.width = grid::unit(0.62, "cm"),
    legend.text = element_text(size = size),
    legend.margin = margin(3, 4, 3, 4)
  )
}

panel_tag <- function(tag) {
  annotate("text", x=-Inf, y=Inf, label=tag, hjust=-0.25, vjust=1.25,
           fontface="bold", size=3.5)
}

panel_tag_once <- function(tag) {
  labs(tag=tag) + theme(
    plot.tag=element_text(face="bold", size=10),
    plot.tag.position=c(0, 1)
  )
}

save_pub <- function(p, stem, width, height) {
  dir.create(file.path("results", "figures"), showWarnings=FALSE, recursive=TRUE)
  ggsave(file.path("results", "figures", paste0(stem, ".pdf")), p,
         width=width, height=height, units="in", device="pdf")
  ggsave(file.path("results", "figures", paste0(stem, ".png")), p,
         width=width, height=height, units="in", dpi=600, bg="white")
}

read_stage_positions <- function() {
  s <- read_csv("data/figure_source/fig2_primary_isi_staging.csv", show_col_types=FALSE) %>%
    filter(path_family == "coupled", subset == "core_secure_pairs") %>%
    mutate(stage = recode(metric,
      wt_exit_p_isi="WT exit",
      balance_p_isi="Balance",
      sca3_entry_p_isi="SCA3 entry"
    )) %>%
    filter(!is.na(stage))
  s
}

add_stage_vlines <- function(p, alpha=.55) {
  s <- read_stage_positions()
  p + geom_vline(data=s, aes(xintercept=median, linetype=stage),
                 color="#666666", linewidth=.45, alpha=alpha,
                 show.legend=FALSE) +
    scale_linetype_manual(values=stage_lty, breaks=names(stage_lty))
}
