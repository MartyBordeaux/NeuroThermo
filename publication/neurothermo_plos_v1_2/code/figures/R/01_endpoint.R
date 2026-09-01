source("code/figures/R/00_theme.R")

make_fig1 <- function() {
  d <- read_csv("data/figure_source/fig1_endpoint_cells.csv", show_col_types=FALSE) %>%
    mutate(group=factor(group, levels=c("WT","SCA3")))

  dot_panel <- function(var, ylab, logy=FALSE, tag="") {
    p <- ggplot(d, aes(x=group, y=.data[[var]], color=group, group=group)) +
      geom_jitter(width=.08, height=0, size=2.1, alpha=.85) +
      stat_summary(fun.data=qfun, geom="errorbar", width=.16, linewidth=.7, color="black") +
      stat_summary(fun=median, geom="point", shape=95, size=8,
                   color="black") +
      scale_color_manual(values=group_cols) +
      labs(x=NULL, y=ylab) + theme_nt() + theme(legend.position="none") + panel_tag(tag)
    if (logy) p <- p + scale_y_log10()
    p
  }

  pA <- dot_panel("rheobase_J_best", expression(J[rheo]~"(pA/pF)"), TRUE, "A")
  pB <- dot_panel("exp_q75_firing_rate_hz", expression("Experimental firing rate at "*q==0.75*" (Hz)"), FALSE, "B")
  pC <- dot_panel("exp_q75_mean_isi_ms", expression("Experimental mean ISI at "*q==0.75*" (ms)"), FALSE, "C")

  val <- bind_rows(
    d %>% transmute(group, cell_id, metric="Firing rate", experiment=exp_q75_firing_rate_hz, model=model_q75_firing_rate_hz),
    d %>% transmute(group, cell_id, metric="Mean ISI", experiment=exp_q75_mean_isi_ms, model=model_q75_mean_isi_ms)
  ) %>% filter(is.finite(experiment), is.finite(model))

  pD <- ggplot(val, aes(experiment, model, color=group)) +
    geom_abline(slope=1, intercept=0, linewidth=.5, linetype="dashed", color="#555555") +
    geom_point(size=2.0, alpha=.85) +
    facet_wrap(~metric, scales="free") +
    scale_color_manual(values=group_cols) +
    guides(color=guide_legend(override.aes=list(size=3.2, alpha=1))) +
    labs(x="Experiment", y="Best-fit HR model") +
    theme_nt() + panel_tag_once("D") +
    legend_inset_theme(x=.32, y=.08, size=10.4)

  p <- (pA | pB | pC) / pD
  save_pub(p, "Fig1_endpoint_phenotype", 7.2, 6.3)
  p
}
