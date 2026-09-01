source("code/figures/R/00_theme.R")

make_fig3 <- function() {
  surf <- read_csv("data/figure_source/fig3_drive_surface_core_secure.csv", show_col_types=FALSE)
  comp <- read_csv("data/figure_source/fig3_coupled_component_sensitivity_core_secure.csv", show_col_types=FALSE)
  inter <- read_csv("data/figure_source/fig3_interaction_at_boundaries_core_secure.csv", show_col_types=FALSE) %>%
    mutate(stage=factor(recode(stage, WT_exit="WT exit", balance="Balance", SCA3_entry="SCA3 entry"),
                        levels=c("WT exit","Balance","SCA3 entry")))
  ref <- read_csv("data/figure_source/fig2_projection_reference.csv", show_col_types=FALSE) %>%
    filter(projection=="isi_primary_v1_0_frozen")
  thresholds <- c(ref$wt_exit_A_threshold[1], .5, ref$sca3_entry_A_threshold[1])
  stage_pos <- read_stage_positions()

  pA <- ggplot(surf, aes(p_intrinsic, p_drive, fill=A_isi_median)) +
    geom_tile() +
    geom_contour(aes(z=A_isi_median), breaks=thresholds, color="white", linewidth=.45) +
    annotate("point", x=0, y=0, color=COL_WT, size=2.7) +
    annotate("point", x=1, y=1, color=COL_SCA3, size=2.7) +
    annotate("text", x=.04, y=.07, label="WT", color="white", fontface="bold", hjust=0, size=4.0) +
    annotate("text", x=.96, y=.93, label="SCA3", color="white", fontface="bold", hjust=1, size=4.0) +
    annotate("text", x=.04, y=.07, label="WT", color=COL_WT, fontface="bold", hjust=0, size=3) +
    annotate("text", x=.96, y=.93, label="SCA3", color=COL_SCA3, fontface="bold", hjust=1, size=3) +
    scale_fill_gradient2(low="#2B6CB0", mid="white", high="#C44E52", midpoint=.5) +
    coord_equal() + labs(x=expression(p[intrinsic]), y=expression(p[drive]), fill=expression(A[ISI])) +
    theme_nt() + panel_tag("A")

  pB <- ggplot(surf, aes(p_intrinsic, p_drive, fill=drive_dominance_fraction_median)) +
    geom_tile() +
    annotate("point", x=0, y=0, color=COL_WT, size=2.7) +
    annotate("point", x=1, y=1, color=COL_SCA3, size=2.7) +
    annotate("text", x=.04, y=.07, label="WT", color=COL_WT, fontface="bold", hjust=0, size=3) +
    annotate("text", x=.96, y=.93, label="SCA3", color=COL_SCA3, fontface="bold", hjust=1, size=3) +
    scale_fill_gradient2(low="#2B6CB0", mid="white", high="#C44E52", midpoint=.5,
                         limits=c(0,1), oob=squish) +
    coord_equal() + labs(x=expression(p[intrinsic]), y=expression(p[drive]), fill="Drive dominance") +
    theme_nt() + panel_tag("B")

  combined_diag <- surf %>%
    filter(abs(p_intrinsic - p_drive) < 1e-10) %>%
    transmute(p=p_intrinsic, component="Combined drive",
              median=dA_ddrive_median, q25=dA_ddrive_q25, q75=dA_ddrive_q75)

  comp_long <- bind_rows(
    combined_diag,
    comp %>% transmute(p=p_intrinsic, component="kappa_I", median=dA_dkappa_median, q25=dA_dkappa_q25, q75=dA_dkappa_q75),
    comp %>% transmute(p=p_intrinsic, component="J", median=dA_dJ_median, q25=dA_dJ_q25, q75=dA_dJ_q75)
  )
  comp_long$component <- factor(comp_long$component, levels=c("Combined drive","kappa_I","J"))
  cc <- c("Combined drive"=COL_COMBINED, "kappa_I"=COL_KAPPA, "J"=COL_J)
  comp_ribbon <- comp_long %>%
    filter(is.finite(q25), is.finite(q75))

  pC <- ggplot(comp_long, aes(p, median, color=component, group=component)) +
    annotate("rect", xmin=0, xmax=.30, ymin=-Inf, ymax=Inf,
             fill=COL_J, alpha=.035) +
    annotate("rect", xmin=.70, xmax=1, ymin=-Inf, ymax=Inf,
             fill=COL_KAPPA, alpha=.035) +
    geom_hline(yintercept=0, color="#666666", linewidth=.35) +
    geom_ribbon(
      data=comp_ribbon,
      aes(ymin=q25, ymax=q75, fill=component, group=component),
      alpha=.10, color=NA, show.legend=FALSE
    ) +
    geom_line(linewidth=.85, na.rm=TRUE) +
    geom_vline(data=stage_pos, aes(xintercept=median), linetype="dashed", color="#777777", linewidth=.4) +
    scale_color_manual(values=cc, labels=c("Combined drive", "kappa_I", "Applied J")) +
    scale_fill_manual(values=cc, guide="none") +
    annotate("text", x=.15, y=Inf, label="J-shaped", color=COL_J,
             size=2.5, hjust=.5, vjust=1.25) +
    annotate("text", x=.85, y=Inf, label="kappa_I-shaped", color=COL_KAPPA,
             size=2.5, hjust=.5, vjust=1.25) +
    labs(x="Coupled path progress, p", y="Local sensitivity dA_ISI/dp") +
    theme_nt() + panel_tag("C")

  eff <- bind_rows(
    inter %>% transmute(stage, component="Combined", median=median_abs_delta_combined_median, q25=median_abs_delta_combined_q25, q75=median_abs_delta_combined_q75),
    inter %>% transmute(stage, component="kappa_I", median=median_abs_delta_kappa_median, q25=median_abs_delta_kappa_q25, q75=median_abs_delta_kappa_q75),
    inter %>% transmute(stage, component="J", median=median_abs_delta_J_median, q25=median_abs_delta_J_q25, q75=median_abs_delta_J_q75),
    inter %>% transmute(stage, component="Interaction", median=median_abs_interaction_median, q25=median_abs_interaction_q25, q75=median_abs_interaction_q75)
  ) %>% mutate(component=factor(component, levels=c("Combined","kappa_I","J","Interaction")))
  ce <- c("Combined"=COL_COMBINED, "kappa_I"=COL_KAPPA, "J"=COL_J, "Interaction"=COL_INTERACTION)
  pD <- ggplot(eff, aes(stage, median, color=component, group=component)) +
    geom_linerange(aes(ymin=q25, ymax=q75), position=position_dodge(width=.55), linewidth=.65) +
    geom_point(position=position_dodge(width=.55), size=2.1) +
    scale_color_manual(values=ce, labels=c("Combined", "kappa_I", "Applied J", "Interaction")) +
    labs(x=NULL, y=expression("Median |"*Delta*A[ISI]*"| near boundary")) + theme_nt() + panel_tag("D")

  p <- (pA | pB) / (pC | pD)
  save_pub(p, "Fig3_intrinsic_drive_decomposition", 7.2, 6.6)
  p
}
