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
# tab_et1 <- get_execution_times(pop1)
# tab_et2 <- get_execution_times(pop2)
# tab_et3 <- get_execution_times(pop3)
# tab_et11 <- get_execution_times(pop11)
# tab_et22 <- get_execution_times(pop22)
# save(tab_et1, tab_et2, tab_et3, tab_et11, tab_et22, file = "./data/tab_exec_times.rds")
load(file = "./data/tab_exec_times.rds")

# Create time plots ----

p1_times <- plot_times(tab_et1, title = NULL, legend.position = "none")
p2_times <- plot_times(tab_et2, title = NULL, legend.position = "none")
p3_times <- plot_times(tab_et3, title = NULL, legend.position = "none")
p11_times <- plot_times(tab_et11, title = NULL, legend.position = "none")
p22_times <- plot_times(tab_et22, title = NULL, legend.position = "none")

p1_take <- plot_take(tab_et1, legend.position = "none")
p2_take <- plot_take(tab_et2, legend.position = "none")
p3_take <- plot_take(tab_et3, legend.position = "none")
p11_take <- plot_take(tab_et11, legend.position = "none")
p22_take <- plot_take(tab_et22, legend.position = "none")

# fig 1
fig_3_2 <- p3_times + p2_times + p3_take + p2_take +
  plot_layout(ncol = 2, heights = c(2, 1)) +
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ) &
  theme(legend.justification = "left", legend.text = element_text(face = "italic"))
fig_3_2
ggsave(
  "./figures/fig_rnabox_3_2.pdf",
  fig_3_2,
  device = "pdf",
  dpi = 600,
  width = 12,
  height = 10 / 1.618 # height = 8/1.618
)

# fig 2
fig_11_22 <- p11_times + p22_times + p11_take + p22_take +
  plot_layout(ncol = 2, heights = c(2, 1)) +
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ) &
  theme(legend.justification = "left", legend.text = element_text(face = "italic"))
fig_11_22
ggsave(
  "./figures/fig_rnabox_11_22.pdf",
  fig_11_22,
  device = "pdf",
  dpi = 600,
  width = 12,
  height = 10 / 1.618 # height = 8/1.618
)

# Compute variances V and V0 ----

options(digits = 6)
# (tab_var1 <- get_variances_rounding(pop1))
# (tab_var2 <- get_variances_rounding(pop2))
# (tab_var3 <- get_variances_rounding(pop3))
# (tab_var11 <- get_variances_rounding(pop11))
# (tab_var22 <- get_variances_rounding(pop22))
# save(tab_var1, tab_var2, tab_var3, tab_var11, tab_var22, file = "./data/tab_var.rds")
load(file = "./data/tab_var.rds")

cap <- "\\label{tab:}
\\footnotesize Variances $V$ and $V_0$ are based on optimal (not necessarily integer)
and optimal integer allocations respectively.
For variances $\\tilde{V}$, based on rounded optimal non-integer allocation,
we systematically get $\\tilde{V}/V_0=1$ (up to five decimal digits)."

tab_rvar <- bind_cols(
  dplyr::select(tab_var1, f, n1 = n, rv1 = rv_rnabox),
  dplyr::select(tab_var2, n2 = n, rv2 = rv_rnabox),
  dplyr::select(tab_var11, n11 = n, rv11 = rv_rnabox),
  dplyr::select(tab_var22, n22 = n, rv22 = rv_rnabox)
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
