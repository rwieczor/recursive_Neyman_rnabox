#' Generation of Data for Experiments with Allocation Algorithms
#'
#' @description Artificial population is generated with many strata for
#'   experiments with optimal allocation algorithms in stratified sampling.
#'   Simulated population is used for computing parameters needed for allocation
#'   algorithms.
#'
#' @param Nrep - number of iterations for sequential generation of population,
#'
#' @param seed - seed for reproducible random numbers.
#'
#' @return list with parameters of generated population i.e.
#'  Nh - vector of population sizes in strata,
#'  Sh - vector of population standard deviations in strata.
#'
#' @export
#' @examples
#' pop <- gen_population(Nrep = 20)
#' Nh <- pop$Nh
#' Sh <- pop$Sh
#  NROW(Nh) # number of generated strata
#'
gen_population <- function(Nrep = 10, seed = NULL) {
  require(dplyr)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  d <- NULL
  Ls <- 10 # number of strata in single generation step
  for (i in 1:Nrep) {
    sigma <- log(i + 1)
    cat("iteration ", i, " sigma ", sigma, "\n")
    di <- data.frame(x = rlnorm(10000, sdlog = sigma))
    # di <- data.frame(x= sqrt(abs(rcauchy(10000, scale=sigma))))
    ## di <- data.frame(x=abs(rnorm(10000, sd=1)))

    strata_boundaries <- stratification::strata.geo(x = di$x, CV = 0.05, Ls = Ls)$bh
    ## strata_boundaries<-strata.cumrootf(x=di$x, CV=0.05, Ls=Ls)$bh
    di$h <- (i - 1) * Ls + as.numeric(cut(di$x, breaks = c(0, strata_boundaries, Inf)))
    d <- bind_rows(d, di)
  }
  d$h <- as.factor(d$h)

  pop <- group_by(d, h) %>%
    summarise(Sh = sqrt(var(x)), Nh = n()) %>%
    na.omit()

  # additional random shuffling of data
  s <- sample(1:nrow(pop), nrow(pop))
  pop <- pop[s, ]

  list(Nh = pop$Nh, Sh = pop$Sh)
}

# generates artificial populations and box constraints.
# pop_n - population number: 1, 11, 2, 22, or 3.
gen_population_boxcnstr <- function(pop_n = 1) {
  stopifnot(pop_n %in% c(1, 2, 3, 11, 22))

  if (pop_n == 3) {
    # Simple population as used in article in JSSaM (2021, fig.2).
    set.seed(4321)
    Nh <- rep(1000, 20)
    Sh <- 10^(1:20)
    Sh <- sample(Sh, 20)
    mh <- rep(100, 20)
    Mh <- Nh
  } else {
    pop <- if (pop_n == 1 || pop_n == 11) {
      gen_population(Nrep = 100, seed = 2234)
    } else { # i.e. pop_n = 2 || pop_n = 22
      gen_population(Nrep = 200, seed = 876)
    }
    Nh <- pop$Nh
    Sh <- pop$Sh
    Mh <- Nh # upper bounds

    if (pop_n == 1 || pop_n == 2) {
      mh <- integer(length(Nh)) # initial lower bounds = 0
      n <- round(0.1 * sum(Nh))
      alc <- stratallo:::CapacityScaling(n, Nh * Sh, mh = mh, Mh = Mh)
      mh <- pmin(100, Mh) # additional lower constraints
      ix <- which(alc > 2 & mh != Mh)
    } else {
      mh <- pmin(100, Mh) # lower bounds
      ix <- (mh >= 2) & (mh != Mh)
    }
    Nh <- Nh[ix]
    Sh <- Sh[ix]
    mh <- mh[ix]
    Mh <- Mh[ix]
  }
  Ah <- Sh * Nh
  data.frame(Nh = Nh, Sh = Sh, Ah = Ah, mh = mh, Mh = Mh)
}

# Creates data with variances for selected algorithms and different fractions
# (examination of rounding effect)
get_variances_rounding <- function(pop, fractions = seq(0.1, 0.9, 0.1)) {
  Nh <- pop$Nh
  Sh <- pop$Sh
  Ah <- pop$Ah
  mh <- pop$mh
  Mh <- pop$Mh
  N <- sum(Nh)

  tab <- NULL
  for (f in fractions) {
    print(paste("Fraction:", f))
    n <- round(f * N)

    if (n >= sum(mh) && n <= sum(Mh)) {
      # CapacityScaling
      alloc_cs <- stratallo:::CapacityScaling(n, Ah, mh = mh, Mh = Mh)
      var_cs <- stratallo::var_st_tsi(alloc_cs, Nh, Sh)

      # RNABOX
      alloc_rnabox <- stratallo::rnabox(n, Ah, Mh, mh)
      alloc_rnabox_round <- stratallo::round_oric(alloc_rnabox)
      var_rnabox <- stratallo::var_st_tsi(alloc_rnabox, Nh, Sh)
      var_rnabox_round <- stratallo::var_st_tsi(alloc_rnabox_round, Nh, Sh)

      tabi <- data.frame(
        N = N, H = length(Ah), f = f, n = n, n_round = sum(alloc_rnabox_round),
        rv_rnabox = var_rnabox / var_cs,
        rv_rnabox_round = var_rnabox_round / var_cs
      )
      tab <- bind_rows(tab, tabi)
    }
  }
  tab
}

# Creates data with times for selected algorithms and different fractions
get_execution_times <- function(pop, time_unit = "milliseconds") {
  Ah <- pop$Ah
  Nh <- pop$Nh
  mh <- pop$mh
  Mh <- pop$Mh

  N <- sum(Nh)
  sum_m <- sum(mh)
  sum_M <- sum(Mh)

  tab_exec <- NULL
  tab_niter <- NULL
  samples_fractions <- seq(sum_m / N, sum_M / N, 0.1)
  for (f in samples_fractions) {
    print(paste("Fraction:", f))
    n <- round(f * N)

    if (n > sum_m && n < sum_M) {
      alc <- stratallo:::CapacityScaling(n, Ah, mh = mh, Mh = Mh)
      n_take_min <- sum(alc <= mh) # number of take-min strata
      n_take_max <- sum(alc >= Mh) # number of take-max strata
      n_take_minmax <- sum(mh == Mh)
      n_take_min <- n_take_min - n_take_minmax
      n_take_max <- n_take_max - n_take_minmax
      n_take_Neyman <- sum(alc > mh & alc < Mh)

      # RNABOX
      alc_rnabox <- round(stratallo::rnabox(n, Ah, Mh, mh))
      if (max(abs((alc_rnabox - alc))) > 1) {
        stop("Bad allocation rnabox!")
      }
      alc_rnabox_debug <- stratallo:::rnabox_debug(n, Ah, Mh, mh)
      tab_rnabox_niter <- stratallo:::rnabox_debug_summary(alc_rnabox_debug$details)
      tab_rnabox_niter <- data.frame(
        Algorithm = "RNABOX",
        f = f,
        n_iter = paste0("(", paste0(tab_rnabox_niter$tB1_iter, collapse = ","), ")")
      )

      # FPIA
      alc_fpi <- stratallo:::fpia(n, Ah, mh, Mh)
      alc_fpi_nh <- round(alc_fpi$nh)
      if (max(abs((alc_fpi_nh - alc))) > 1) {
        stop("Bad allocation stratallo:::fpia!")
      }
      tab_fpia_niter <- data.frame(
        Algorithm = "FPIA",
        f = f,
        n_iter = paste0("(", alc_fpi$iter, ")")
      )

      ex <- microbenchmark(
        FPIA = stratallo:::fpia(n, Ah, mh, Mh)$nh,
        RNABOX = stratallo::rnabox(n, Ah, Mh, mh),
        times = 100,
        unit = time_unit
      )

      # summary(ex)
      # autoplot(ex)
      exi <- rename(
        summary(ex)[, c("expr", "mean", "median")],
        Algorithm = expr,
        Mean_time = mean,
        Median_time = median
      ) %>%
        group_by(Algorithm) %>%
        mutate(
          sum_m = sum_m,
          sum_M = sum_M,
          N = N,
          f = f,
          H = length(Ah),
          n_take_min = n_take_min,
          n_take_max = n_take_max,
          n_take_minmax = n_take_minmax,
          n_take_Neyman = n_take_Neyman
        )
      tab_exec <- bind_rows(tab_exec, exi)
      tab_niter <- bind_rows(tab_niter, tab_rnabox_niter, tab_fpia_niter)
    } else {
      sprintf("Skipped infeasbile f = %f", f)
    }
  }
  left_join(tab_exec, tab_niter, by = c("Algorithm", "f"))
}

plot_times <- function(data,
                       legend.position = "right",
                       title = "Time comparison of selected algorithms",
                       x_lab = "Total sample size (as fraction of N)",
                       y_lab = "Median Time [miliseconds]") {
  data <- data %>% mutate(
    population = paste(H, "strata, N =", N),
    f = round(f, 2),
    flab = paste0("[", as.character(n_take_min), ",", as.character(n_take_max), "]")
  )

  p <- ggplot(data = data, aes(x = f, y = Median_time)) +
    geom_point(size = 2, aes(shape = Algorithm)) +
    geom_line(aes(linetype = Algorithm)) +
    geom_text_repel(data = data, aes(x = f, y = Median_time, label = n_iter)) +
    facet_wrap(~population, scales = "free_y") +
    labs(
      x = x_lab,
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

plot_take <- function(data, legend.position = "right", x_lab = "Total sample size (as fraction of N)",
                      y_lab = "Number of strata") {
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
    c("n_take_max", "n_take_Neyman", "n_take_min", "n_take_minmax"),
    factor_key = TRUE
  )
  tab_take <- tab_take[tab_take$value != 0, ]
  tab_take$pct <- round((tab_take$value / tab_take$H) * 100, 1)

  # tab_take$series must be a factor with levels: "n_take_max", "n_take_Neyman", "n_take_max".
  p <- ggplot(data = tab_take, mapping = aes(x = f, y = value)) +
    geom_bar(
      mapping = aes(fill = series), position = "stack", stat = "identity", width = 0.03
    ) +
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
    labs(x = x_lab, y = y_lab) +
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
