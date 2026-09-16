source("code/figures/R/00_theme.R")

make_figS1 <- function() {
  g <- read_csv("data/figure_source/supp_s1_dynamic_group.csv", show_col_types=FALSE) %>%
    mutate(
      group=factor(group, levels=c("WT","SCA3")),
      source=recode(source, experiment="Experiment", model_best="HR model"),
      source=factor(source, levels=c("Experiment","HR model"))
    )
  cells <- read_csv("data/figure_source/supp_s1_dynamic_cells.csv", show_col_types=FALSE) %>%
    mutate(group=factor(group, levels=c("WT","SCA3")))
  sup <- read_csv("data/figure_source/supp_s1_support_summary.csv", show_col_types=FALSE) %>%
    mutate(group=factor(group, levels=c("WT","SCA3")), qf=factor(q, levels=c(.25,.5,.75)))

  common <- list(
    scale_color_manual(values=group_cols),
    scale_linetype_manual(values=c("Experiment"="solid", "HR model"="22")),
    scale_shape_manual(values=c("Experiment"=16, "HR model"=1)),
    theme_nt()
  )

  pA <- ggplot(g, aes(q, firing_rate_hz_median, color=group, linetype=source, shape=source)) +
    geom_line(linewidth=.8) + geom_point(size=2.0) + common +
    labs(x="Support-restricted current coordinate, q", y="Firing rate (Hz)") + panel_tag("A")
  pB <- ggplot(g, aes(q, mean_isi_ms_median, color=group, linetype=source, shape=source)) +
    geom_line(linewidth=.8) + geom_point(size=2.0) + common +
    labs(x="Support-restricted current coordinate, q", y="Mean ISI (ms)") + panel_tag("B")

  pC <- ggplot(cells, aes(q, firing_rate_hz, group=cell_id, color=group)) +
    geom_line(alpha=.32, linewidth=.45) + geom_point(alpha=.55, size=1.1) +
    facet_wrap(~group) + scale_color_manual(values=group_cols) +
    labs(x="q", y="Experimental firing rate (Hz)") + theme_nt() + theme(legend.position="none") + panel_tag_once("C")

  pD <- ggplot(sup, aes(qf, fraction_supported, fill=group)) +
    geom_col(position=position_dodge(width=.75), width=.65) +
    geom_text(aes(label=paste0(n_supported,"/",n_cells)), position=position_dodge(width=.75),
              vjust=-.25, size=2.7) +
    scale_fill_manual(values=group_cols) + scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1.12)) +
    labs(x="q", y="Cells within observed-current support") + theme_nt() + panel_tag("D")

  p <- (pA | pB) / (pC | pD)
  save_pub(p, "FigS1_support_restricted_dynamics", 7.2, 6.3)
  p
}
