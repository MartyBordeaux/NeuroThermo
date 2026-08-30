source("code/figures/R/00_theme.R")

make_fig2 <- function() {
  curves <- read_csv("data/figure_source/fig2_core_secure_curves.csv", show_col_types=FALSE) %>%
    mutate(path = recode(path_family,
      drive_early="Drive early", coupled="Coupled", drive_late="Drive late"),
      path=factor(path, levels=c("Drive early","Coupled","Drive late")))
  stages <- read_csv("data/figure_source/fig2_primary_isi_staging.csv", show_col_types=FALSE) %>%
    mutate(
      path=recode(path_family, drive_early="Drive early", coupled="Coupled", drive_late="Drive late"),
      path=factor(path, levels=c("Drive early","Coupled","Drive late")),
      stage=recode(metric, wt_exit_p_isi="WT exit", balance_p_isi="Balance", sca3_entry_p_isi="SCA3 entry"),
      stage=factor(stage, levels=c("WT exit","Balance","SCA3 entry"))
    ) %>% filter(subset=="core_secure_pairs", !is.na(stage))
  ref <- read_csv("data/figure_source/fig2_projection_reference.csv", show_col_types=FALSE) %>%
    filter(projection=="isi_primary_v1_0_frozen")
  thresholds <- tibble(
    stage=factor(c("WT exit","Balance","SCA3 entry"), levels=c("WT exit","Balance","SCA3 entry")),
    A=c(ref$wt_exit_A_threshold[1], .5, ref$sca3_entry_A_threshold[1])
  )

  prot <- tibble(p=seq(0,1,length.out=201)) %>%
    mutate(`Drive early`=2*p-p^2, Coupled=p, `Drive late`=p^2) %>%
    pivot_longer(-p, names_to="path", values_to="p_drive") %>%
    mutate(path=factor(path, levels=c("Drive early","Coupled","Drive late")))

  pA <- ggplot(prot, aes(p, p_drive, color=path)) +
    geom_abline(slope=1, intercept=0, linewidth=.35, color="#BBBBBB") +
    geom_line(linewidth=1.0) +
    annotate("point", x=0, y=0, color=COL_WT, fill=COL_WT, size=2.7) +
    annotate("point", x=1, y=1, color=COL_SCA3, fill=COL_SCA3, size=2.7) +
    annotate("text", x=.18, y=.07, label="WT", color=COL_WT, fontface="bold", hjust=.5, size=3) +
    annotate("text", x=.82, y=.93, label="SCA3", color=COL_SCA3, fontface="bold", hjust=.5, size=3) +
    scale_color_manual(values=path_cols) + coord_equal() +
    labs(x=expression(p[intrinsic]), y=expression(p[drive])) + theme_nt() + panel_tag("A")

  pB <- ggplot(curves, aes(path_progress, A_isi_median, color=path, fill=path)) +
    geom_ribbon(aes(ymin=A_isi_q25, ymax=A_isi_q75), alpha=.12, color=NA) +
    geom_line(linewidth=.9) +
    geom_hline(data=thresholds, aes(yintercept=A, linetype=stage), color="#4D4D4D", linewidth=.45) +
    annotate("point", x=0, y=median(curves$A_isi_median[curves$path_progress == 0]), color=COL_WT, size=2.7) +
    annotate("point", x=1, y=median(curves$A_isi_median[curves$path_progress == 1]), color=COL_SCA3, size=2.7) +
    annotate("text", x=.04, y=median(curves$A_isi_median[curves$path_progress == 0])+.08,
             label="WT", color=COL_WT, fontface="bold", hjust=0, size=3) +
    annotate("text", x=.80, y=1.00,
             label="SCA3", color=COL_SCA3, fontface="bold", hjust=.5, vjust=1, size=3) +
    scale_color_manual(values=path_cols) + scale_fill_manual(values=path_cols) +
    scale_linetype_manual(values=c("WT exit"="dotted","Balance"="dashed","SCA3 entry"="dotdash")) +
    labs(x="Path progress, p", y=expression(A[ISI])) + theme_nt() + panel_tag("B")

  pC <- ggplot(stages, aes(path, median, color=stage)) +
    geom_linerange(aes(ymin=q25, ymax=q75), position=position_dodge(width=.42), linewidth=.7) +
    geom_point(position=position_dodge(width=.42), size=2.3) +
    scale_color_manual(values=stage_cols) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    labs(x=NULL, y="Stage location, p") + theme_nt() + panel_tag("C")

  p <- (pA | pB) / pC
  save_pub(p, "Fig2_transition_staging", 7.2, 6.2)
  p
}
