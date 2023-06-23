source("gen_population.R")

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
      alc <- stratallo:::CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
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
  dh <- Sh * Nh
  list(Nh = Nh, Sh = Sh, dh = dh, mh = mh, Mh = Mh)
}

# Creates data with variances for selected algorithms and different fractions
# (examination of rounding effect)
get_variances_rounding <- function(pop_n = 1) {
  # pop_n - number of population (1, 2 or 3)
  pop <- gen_population_boxcnstr(pop_n)
  Nh <- pop$Nh
  Sh <- pop$Sh
  dh <- pop$dh
  mh <- pop$mh
  Mh <- pop$Mh
  N <- sum(Nh)

  tab <- NULL
  for (f in seq(0.1, 0.5, 0.1)) {
    print(paste("Fraction:", f))
    n <- round(f * N)

    if (n > sum(mh) && n < sum(Mh)) {
      # CapacityScaling
      alc <- stratallo:::CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
      var_cs <- stratallo::var_st_tsi(alc, Nh, Sh)

      # RNABOX
      alc_rnabox <- rnabox(n, dh, mh, Mh)
      alc_rnabox_round <- stratallo::round_oric(alc_rnabox)
      var_rnabox <- stratallo::var_st_tsi(alc_rnabox, Nh, Sh)
      var_rnabox_round <- stratallo::var_st_tsi(alc_rnabox_round, Nh, Sh)

      tabi <- data.frame(
        N = N, H = length(Nh), f = f, n = n, n_round = sum(alc_rnabox_round),
        rv_rnabox = var_rnabox / var_cs,
        rv_rnabox_round = var_rnabox_round / var_cs
      )
      tab <- bind_rows(tab, tabi)
    }
  }
  tab
}

# Creates data with times for selected algorithms and different fractions
get_execution_times <- function(pop_n = 1, time_unit = "milliseconds") {
  # pop_n - number of population (1,2 or 3)
  pop <- gen_population_boxcnstr(pop_n)
  Nh <- pop$Nh
  Sh <- pop$Sh
  mh <- pop$mh
  Mh <- pop$Mh
  dh <- pop$dh

  N <- sum(Nh)
  sum_m <- sum(mh)
  sum_M <- sum(Mh)

  tab <- NULL
  samples_fractions <- seq(sum_m / N, sum_M / N, 0.1)
  for (f in samples_fractions) {
    print(paste("Fraction:", f))
    n <- round(f * N)

    if (n > sum_m && n < sum_M) {
      alc <- stratallo:::CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
      n_take_min <- sum(alc <= mh) # number of take-min strata
      n_take_max <- sum(alc >= Mh) # number of take-max strata
      n_take_minmax <- sum(mh == Mh)
      n_take_min <- n_take_min - n_take_minmax
      n_take_max <- n_take_max - n_take_minmax
      n_take_Neyman <- sum(alc > mh & alc < Mh)

      alc_rnabox <- round(rnabox(n, dh, mh, Mh))
      if (max(abs((alc_rnabox - alc))) > 1) {
        stop("Bad allocation rnabox!")
      }

      alc_fpi <- round(stratallo:::fpia(n, Nh, Sh, mh, Mh)$nh)
      if (max(abs((alc_fpi - alc))) > 1) {
        stop("Bad allocation stratallo:::fpia!")
      }

      ex <- microbenchmark(
        FPIA = stratallo:::fpia(n, Nh, Sh, mh, Mh)$nh,
        RNABOX = rnabox(n, dh, mh, Mh),
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
          H = length(Nh),
          n_take_min = n_take_min,
          n_take_max = n_take_max,
          n_take_minmax = n_take_minmax,
          n_take_Neyman = n_take_Neyman
        )
      tab <- bind_rows(tab, exi)
    } else {
      sprintf("Skipped infeasbile f = %f", f)
    }
  }
  tab
}

plot_times <- function(data,
                       legend.position = "right",
                       title = "Time comparison of selected algorithms",
                       y_lab = "Time [miliseconds]") {
  data <- data %>% mutate(
    population = paste(H, "strata, N =", N),
    f = round(f, 2),
    flab = paste0("[", as.character(n_take_min), ",", as.character(n_take_max), "]")
  )

  p <- ggplot(data = data, aes(x = f, y = Median_time)) +
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