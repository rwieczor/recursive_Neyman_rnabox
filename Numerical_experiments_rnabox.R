
# R code with numerical experiments from paper (in preparation):
# Wesolowski J., Wieczorkowski R., Wojciak W. (2023?),
# Recursive Neyman Algorithm for Optimum Sample Allocation under Box Constraints on Sample Sizes in Strata.

# Load libraries ----
library(dplyr)
library(microbenchmark)
library(bench)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)

# R codes of the algorithms ----
source("CapacityScaling.R")
source("SimpleGreedy.R")
source("fpia.R")
source("gen_population.R")
rnabox <- stratallo::rnabox

# Custom functions ----

# random rounding
ran_round <- function(x) {
  floor(x) + (runif(length(x)) < (x - floor(x)))
}

# rounding based on article:
# Rama Count, Massoud Heidari, Optimal rounding under integer constraints
# December 2014, arxiv
round_oric <- function(x) {
  n <- round(sum(x))
  m <- floor(x)
  y <- x - m
  Ix <- sum(y)

  if (Ix == 0) {
    return(x)
  } else {
    iy <- order(-y)
    u <- unique(y[iy])
    z <- integer(length(x))
    for (i in 1:length(u)) z[iy] <- z[iy] + (y[iy] == u[i]) * i
    iy2 <- order(-y, z, -m)
    # m[iy][iy2][1:Ix] <- ceiling(x[iy][iy2][1:Ix])
    m[iy2][1:Ix] <- (m[iy2][1:Ix]) + 1
    return(m)
  }
}

# generates artificial populations and box constraints.
# pop_n - population number: 1, 2 or 3.
generate_population <- function(pop_n = 1) {
  if (pop_n %in% c(1, 2)) {
    if (pop_n == 1) {
      Nrep <- 100
      set.seed(2234)
    }
    if (pop_n == 2) {
      Nrep <- 200
      set.seed(876)
    }

    pop <- gen_population(Nrep = Nrep)
    Nh <- pop$Nh
    Sh <- pop$Sh
    NROW(Nh)
    # plot(Nh*Sh)
    # hist(Sh*Nh)

    # Generation of lower and upper bounds - based on method used in article:
    # Jacek Wesołowski, Robert Wieczorkowski, Wojciech Wójciak,
    # Optimality of the recursive Neyman allocation, Journal of Survey Statistics and
    #   Methodology (2022) 10, 1263–1275.
    # (added part with lower bounds)

    # improving population to have allocation with values greater than 0
    # using integer allocation algorithm
    mh <- rep(0, length(Nh)) # lower bounds
    Mh <- Nh # upper bounds
    (n <- round(0.1 * sum(Nh)))
    (alc <- CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh))
    sum(alc)

    length(alc[alc > 2]) / length(alc)
    ix <- which(alc > 2)
    Nh <- Nh[ix]
    Sh <- Sh[ix]
    dh <- Nh * Sh
    mh <- mh[ix]
    Mh <- Mh[ix]
    mh <- pmin(100, Mh) # additional lower constraints

    # removing cases with mh=Mh
    ix <- which(!(mh == Mh))
    Nh <- Nh[ix]
    Sh <- Sh[ix]
    dh <- Nh * Sh
    mh <- mh[ix]
    Mh <- Mh[ix]

    (s1 <- sum(mh))
    (s2 <- sum(Mh))
  } else if (pop_n == 3) {
    # Generation of a simple population used in article published in JSSaM (2021, fig.2)
    # added lower constraints
    set.seed(4321)
    Nh <- rep(1000, 20)
    Sh <- 10^(1:20)
    Sh <- sample(Sh, 20)
    mh <- rep(100, 20)
    Mh <- Nh
    dh <- Sh * Nh
  } else {
    print("pop_n must take on a single value from {1, 2, 3}")
  }

  list(Nh = Nh, Sh = Sh, mh = mh, Mh = Mh, dh = dh)
}


# Creates data with variances for selected algorithms and different fractions
# (examination of rounding effect)
get_variances_rounding <- function(pop_n = 1) {
  # pop_n - number of population (1,2 or 3)
  pop <- generate_population(pop_n)
  Nh <- pop$Nh
  Sh <- pop$Sh
  mh <- pop$mh
  Mh <- pop$Mh
  dh <- pop$dh
  
  
  # variance estimation for given allocation
  varal <- function(Nh,Sh,nh)
  {
    return( sum(Nh * (Nh - nh) * Sh*Sh / nh) )
  }
  
  
  tab <- NULL
  
  sum_m <- sum(mh)
  sum_M <- sum(Mh)
  N <- sum(Nh)
  dh <- Nh * Sh
  
  for (f in seq(0.1,0.5,0.1)) {
  #for (f in seq(sum_m / N, sum_M / N, 0.1)) {
    print(paste("Fraction:", f))
    N <- sum(Nh)
    n <- round(f * N)
    
    if (n > sum_m && n < sum_M) {
      alc <- CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
      V0<-varal(Nh,Sh,alc)
      
      al_rnabox <- (stratallo::rnabox(n, dh, mh, Mh))
      al_rnabox_round <- round_oric(al_rnabox)
      #al_fpi <- fpia(n, Nh, Sh, mh, Mh)$nh
      v_rnabox <- varal(Nh,Sh,al_rnabox)
      v_rnabox_round <- varal(Nh,Sh,al_rnabox_round)
  
      tabi <- data.frame(N=N,H=length(Nh),f=f,n=n,n_round=sum(al_rnabox_round),
                         rv_rnabox=v_rnabox/V0,
                         rv_rnabox_round=v_rnabox_round/V0
                         )
      
      tab <- bind_rows(tab, tabi)
    }
  }
  
  
  return(tab)
}


options(digits = 6)
(tab1 <- get_variances_rounding(pop_n=1))
(tab2 <- get_variances_rounding(pop_n=2))
(tab3 <- get_variances_rounding(pop_n=3))


tab <- bind_cols(dplyr::select(tab1,f,n1=n,rv1=rv_rnabox),
       dplyr::select(tab2,n2=n,rv2=rv_rnabox))
colnames(tab)<-c("fraction","$n$","$V/V_0$","$n$","$V/V_0$")
tab_tex <- kable(tab,format="latex",digits=6,align="r",
                 caption = "\\label{tab:} 	
\\footnotesize Variances $V$ and $V_0$ are based on optimal non-integer
and optimal integer allocations respectively. 
For variances $\\tilde{V}$, based on rounded optimal non-integer allocation,  
we systematically get $\\tilde{V}/V_0=1$ (up to five decimal digits).")

cat(tab_tex,file="tab_ratio_variances.tex")

# additional lines to insert after line with text: \begin{tabular}{|r|r|r|r|r|}
# \hline
# sample  & \multicolumn{2}{|c|}{$H=373$} & \multicolumn{2}{|c|}{$H=691$} \\
# \cline{2-5}




# Creates data with times for selected algorithms and different fractions
get_execution_times <- function(pop_n = 1) {
  # pop_n - number of population (1,2 or 3)
  pop <- generate_population(pop_n)
  Nh <- pop$Nh
  Sh <- pop$Sh
  mh <- pop$mh
  Mh <- pop$Mh
  dh <- pop$dh

  tab <- NULL

  sum_m <- sum(mh)
  sum_M <- sum(Mh)
  N <- sum(Nh)
  dh <- Nh * Sh

  # for (f in seq(0.01,0.8,0.05)) {
  for (f in seq(sum_m / N, sum_M / N, 0.1)) {
    print(paste("Fraction:", f))
    N <- sum(Nh)
    n <- round(f * N)

    if (n > sum_m && n < sum_M) {
      alc <- CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
      n_take_min <- sum(alc <= mh) # number of take-all strata
      n_take_max <- sum(alc >= Mh) # number of take-all strata
      n_take_min_max <- sum(mh == Mh)
      n_take_min <- n_take_min - n_take_min_max
      n_take_max <- n_take_max - n_take_min_max
      n_take_Neyman <- sum(alc > mh & alc < Mh)

      al_rnabox <- round(rnabox(n, dh, mh, Mh))
      if (max(abs((al_rnabox - alc))) > 1) {
        stop("Bad allocation rnabox!")
      }

      al_fpi <- round(fpia(n, Nh, Sh, mh, Mh)$nh)
      if (max(abs((al_fpi - alc))) > 1) {
        stop("Bad allocation fpia!")
      }

      # options(digits = 10)
      ex <- microbenchmark(
        fpia = fpia(n, Nh, Sh, mh, Mh)$nh,
        RNABOX = stratallo::rnabox(n, dh, mh, Mh),
        times = 100,
        unit = "ns"
      )
      # summary(ex)
      # autoplot(ex)
      exi <- rename(ex, Algorithm = expr)
      exi <- exi %>%
        group_by(Algorithm) %>%
        summarise(
          Median_time = median(time) / 1e6, # from nanoseconds to miliseconds
          Mean_time = mean(time) / 1e6
        ) %>%
        mutate(
          sum_m = sum_m,
          sum_M = sum_M,
          N = N,
          f = f,
          H = length(Nh),
          n_take_min = n_take_min,
          n_take_max = n_take_max,
          n_take_min_max = n_take_min_max,
          n_take_Neyman = n_take_Neyman
        )

      tab <- bind_rows(tab, exi)
    }
  }
  tab
}

plot_times <- function(data,
                       legend.position = "right",
                       title = "Time comparison of selected algorithms",
                       # title = " ",
                       y_lab = "Time [miliseconds]") {
  data <- data %>% mutate(
    population = paste(H, "strata, N =", N),
    f = round(f, 2),
    flab = paste0("[", as.character(n_take_min), ",", as.character(n_take_max), "]")
  )

  data %>%
    group_by(Algorithm) %>%
    summarise(mean(Median_time), mean(Mean_time))

  p <-
    ggplot(data = data, aes(x = f, y = Median_time)) +
    geom_point(size = 2, aes(shape = Algorithm)) +
    geom_line(aes(linetype = Algorithm)) +
    facet_wrap(~population, scales = "free_y") +
    labs(
      x = "Sample fraction",
      y = y_lab,
      color = "Algorithms: ",
      title = title
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = legend.position
    ) +
    coord_cartesian(xlim = c(0.1, 1)) +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6))
  p
}

plot_take <- function(data, legend.position = "right", y_lab = "Number of strata") {
  data <- data %>% mutate(
    population = paste(H, "strata, N =", N),
    f = round(f, 2),
    flab = paste0("[", as.character(n_take_min), ",", as.character(n_take_max), "]")
  )
  H <- unique(data$H)

  tab_take <- tidyr::gather(
    data[data$Algorithm == "RNABOX", ],
    "series",
    "value",
    c("n_take_max", "n_take_Neyman", "n_take_min", "n_take_min_max"),
    factor_key = TRUE
  )
  tab_take <- tab_take[tab_take$value != 0, ]
  tab_take$pct <- round((tab_take$value / tab_take$H) * 100, 1)

  # tab_take$series must be a factor with levels: "n_take_max", "n_take_Neyman", "n_take_max".
  p <- ggplot(data = tab_take, mapping = aes(x = f, y = value)) +
    geom_bar(mapping = aes(fill = series), position = "stack", stat = "identity", width = 0.03) +
    facet_wrap(~population) +
    scale_fill_manual(
      drop = TRUE,
      name = "Series",
      values = c("grey20", "grey55", "grey80", "green"),
      labels = c("take-max", "take-Neyman", "take-min", "take-min/take-max \n(m = M)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = legend.position,
      panel.spacing.x = unit(2.5, "lines"),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
    ) +
    labs(x = "Sample fraction", y = y_lab) +
    coord_cartesian(xlim = c(0.1, 1)) +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(10),
      sec.axis = sec_axis(
        trans = ~ (. / H) * 100,
        breaks = seq(0, 100, 20),
        labels = function(x) scales::percent(x, scale = 1)
      )
    )
  p
}

# Generate Populations, Compute Optimal Alloc., Compare Times ----
tab1 <- get_execution_times(pop_n = 1) # population generated with Nrep=100
tab2 <- get_execution_times(pop_n = 2) # population generated with Nrep=200
tab3 <- get_execution_times(pop_n = 3) # a small population

# saveRDS(tab1, "tab1.rds")
# saveRDS(tab2, "tab2.rds")
# saveRDS(tab3, "tab3.rds")
# tab1 <- readRDS("tab1.rds")
# tab2 <- readRDS("tab2.rds")
# tab3 <- readRDS("tab3.rds")

# Create plots ----

p1_times <- plot_times(tab1, title = NULL, legend.position = "none")
p1_take <- plot_take(tab1, legend.position = "none")
p2_times <- plot_times(tab2, title = NULL, y_lab = NULL)
p2_take <- plot_take(tab2, y_lab = NULL)
fig_12 <- p1_times + p2_times + p1_take + p2_take +
  plot_layout(ncol = 2, heights = c(2, 1)) +
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ) &
  theme(legend.justification = "left", legend.text = element_text(face = "italic"))
fig_12
ggsave(
  "rnabox_fig_12.pdf",
  fig_12,
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
  "rnabox_fig_3.pdf",
  fig_3,
  device = "pdf",
  dpi = 600,
  width = 12,
  height = 10 / 1.618 # height = 8/1.618
)
