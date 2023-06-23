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
rnabox <- stratallo::rnabox

# Generate populations, compute optimal alloc., compare times ----
tab1 <- get_execution_times(pop_n = 1) # population generated with Nrep=100
tab2 <- get_execution_times(pop_n = 2) # population generated with Nrep=200
tab3 <- get_execution_times(pop_n = 3) # a small population

tab11 <- get_execution_times(pop_n = 11) # population generated with Nrep=100
tab22 <- get_execution_times(pop_n = 22) # population generated with Nrep=200

# saveRDS(tab1, "./data/tab1.rds")
# saveRDS(tab2, "./data/tab2.rds")
# saveRDS(tab3, "./data/tab3.rds")
# saveRDS(tab11, "./data/tab11.rds")
# saveRDS(tab22, "./data/tab22.rds")
# tab1 <- readRDS("./data/tab1.rds")
# tab2 <- readRDS("./data/tab2.rds")
# tab3 <- readRDS("./data/tab3.rds")
# tab11 <- readRDS("./data/tab11.rds")
# tab22 <- readRDS("./data/tab22.rds")

# Create time plots ----

p1_times <- plot_times(tab1, title = NULL, legend.position = "none")
p1_take <- plot_take(tab1, legend.position = "none")
p22_times <- plot_times(tab22, title = NULL, y_lab = NULL)
p22_take <- plot_take(tab22, y_lab = NULL)
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

p3_times <- plot_times(tab3, title = NULL)
p3_take <- plot_take(tab3)
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
(tab1_var <- get_variances_rounding(pop_n = 1))
(tab2_var <- get_variances_rounding(pop_n = 2))
(tab3_var <- get_variances_rounding(pop_n = 3))
(tab11_var <- get_variances_rounding(pop_n = 11))
(tab22_var <- get_variances_rounding(pop_n = 22))

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
