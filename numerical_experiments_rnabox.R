# R code with numerical experiments from paper (in preparation):
# Wesolowski, J., Wieczorkowski, R., Wojciak, W. (2023?),
# Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints
# on Sample Sizes in Strata.

# Load libraries ----
library(dplyr)
library(microbenchmark)
library(bench)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)

# R codes of helper functions and the RNABOX ----
source("helper_functions.R")

# Generate populations.
# pop1 <- gen_population_boxcnstr(1)
# pop2 <- gen_population_boxcnstr(2)
# pop3 <- gen_population_boxcnstr(3)
# pop11 <- gen_population_boxcnstr(11)
# pop22 <- gen_population_boxcnstr(22)
# save(pop1, pop2, pop3, pop11, pop22, file = "./data/pop.rds")
load(file = "./data/pop.rds")

# Compute optimal allocations, get execution times, etc  ----
# extimes1 <- get_execution_times(pop1)
# extimes2 <- get_execution_times(pop2)
# extimes3 <- get_execution_times(pop3)
# extimes11 <- get_execution_times(pop11)
# extimes22 <- get_execution_times(pop22)
# save(
#   extimes1, extimes2, extimes3, extimes11, extimes22,
#   file = "./data/exec_times.rds"
# )
load(file = "./data/exec_times.rds")

# Create time plots ----

p1_times <- plot_times(extimes1$et, title = NULL, legend.position = "none")
p1_take <- plot_take(extimes1$et, legend.position = "none")
p22_times <- plot_times(extimes22$et, title = NULL, y_lab = NULL)
p22_take <- plot_take(extimes22$et, y_lab = NULL)
fig_1_22 <- p1_times + p22_times + p1_take + p22_take +
  plot_layout(ncol = 2, heights = c(2, 1)) +
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ) &
  theme(legend.justification = "left", legend.text = element_text(face = "italic"))
fig_1_22
ggsave(
  "./figures/fig_rnabox_1_22.pdf",
  fig_1_22,
  device = "pdf",
  dpi = 600,
  width = 12,
  height = 10 / 1.618 # height = 8/1.618
)

p3_times <- plot_times(extimes3$et, title = NULL)
p3_take <- plot_take(extimes3$et)
fig_3 <- p3_times + p3_take + plot_layout(nrow = 2, heights = c(2, 1)) &
  theme(legend.justification = "left", legend.text = element_text(face = "italic"))
fig_3
ggsave(
  "./figures/fig_rnabox_3.pdf",
  fig_3,
  device = "pdf",
  dpi = 600,
  width = 12,
  height = 10 / 1.618 # height = 8/1.618
)

# Compute variances V and V0 ----

options(digits = 6)
(tab1_var <- get_variances_rounding(pop1))
(tab2_var <- get_variances_rounding(pop2))
(tab3_var <- get_variances_rounding(pop3))
(tab11_var <- get_variances_rounding(pop11))
(tab22_var <- get_variances_rounding(pop22))

cap <- "\\label{tab:}
\\footnotesize Variances $V$ and $V_0$ are based on optimal (not necessarily integer)
and optimal integer allocations respectively.
For variances $\\tilde{V}$, based on rounded optimal non-integer allocation,
we systematically get $\\tilde{V}/V_0=1$ (up to five decimal digits)."

tab_rvar <- bind_cols(
  dplyr::select(tab1_var, f, n1 = n, rv1 = rv_rnabox),
  dplyr::select(tab2_var, n2 = n, rv2 = rv_rnabox),
  dplyr::select(tab11_var, n11 = n, rv11 = rv_rnabox),
  dplyr::select(tab22_var, n22 = n, rv22 = rv_rnabox)
)
colnames(tab_rvar) <- c(
  "fraction",
  "$n$", "$V/V_0$",
  "$n$", "$V/V_0$",
  "$n$", "$V/V_0$",
  "$n$", "$V/V_0$"
)
tab_rvar_tex <- knitr::kable(tab_rvar,
  format = "latex", digits = 6, align = "r",
  caption = cap
)
cat(tab_rvar_tex, file = "./data/tab_ratio_variances.tex")

# additional lines to insert after line with text: \begin{tabular}{|r|r|r|r|r|}
# \hline
# sample  & \multicolumn{2}{|c|}{$H=373$} & \multicolumn{2}{|c|}{$H=691$} \\
# \cline{2-5}
