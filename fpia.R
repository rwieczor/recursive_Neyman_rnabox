#' @title Optimal univariate allocation under box constraints for stratified sampling
#' @description Algorithm for optimal allocation in stratified sampling with lower and upper #' constraints based on fixed point iteration
#'
#' @param n target sample size for allocation
#' @param Nh population sizes in strata
#' @param Sh standard deviations for given variable in strata
#' @param mh lower constraints for sample sizes in strata
#' @param Mh upper constraints for sample sizes in strata
#' @param lambda0 initial parameter 'lambda' (optional)
#' @param maxiter maximal number of iterations for algorithm
#' @param tol the desired accuracy (convergence tolerance)
#'
#' @return  vector of optimal allocation sizes,
#'  and number of iterations
#'
#' @references Ralf T. Munnich, Ekkehard W. Sachs, Matthias Wagner (2012),
#' Numerical solution of optimal allocation problems in stratified sampling under box
#'  constraints, AStA Adv Stat Anal (2012) 96:435?450, DOI 10.1007/s10182-011-0176-z
#'
#' @export

fpia <- function(n, Nh, Sh, mh = NULL, Mh = NULL, lambda0 = NULL, maxiter = 100, tol = .Machine$double.eps * 1000) {
  H <- seq_along(Nh)
  dh <- Sh * Nh

  lambda <- if (is.null(lambda0)) {
    (sum(dh)^2) / (n^2) # according to article MSW

    # # initial interval for searching 'lambda' - according to J.Wesolowski
    # r <- c(dh / mh, dh / Mh)
    # a <- min(r)^2
    # b <- max(r)^2
    # lambda <- uniroot(
    #   function(x) glambda(x, n, Nh, Sh, mh, Mh), lower = a, upper = b, maxiter = maxiter
    # )$root
    # lambda <- (a + b) / 2 # the simplest starting value for bisection
  } else {
    lambda0
  }

  iter <- 0
  while (1) {
    iter <- iter + 1
    L <- which((dh^2) / (mh^2) <= lambda)
    U <- which((dh^2) / (Mh^2) >= lambda)
    Hc <- H[-c(L, U)]

    lambda_n <- (sum(dh[Hc]) / (n - sum(mh[L]) - sum(Mh[U])))^2
    if (iter > maxiter || abs(lambda_n - lambda) < tol || is.nan(lambda_n)) {
      break
    }
    lambda <- lambda_n
    # cat("iteracja ",iter," lambda ",lambda,"\n")
  }

  nh <- dh / sqrt(lambda)
  nh[L] <- mh[L]
  nh[U] <- Mh[U]

  # v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
  list(nh = nh, iter = iter)
}

# auxiliary functions

glambda <- function(lambda, n, Nh, Sh, mh = NULL, Mh = NULL) {
  H <- seq_along(Nh)
  dh <- Sh * Nh

  L <- which((dh^2) / (mh^2) <= lambda)
  U <- which((dh^2) / (Mh^2) >= lambda)
  Hc <- H[-c(L, U)]

  nh <- dh / sqrt(lambda)
  nh[L] <- mh[L]
  nh[U] <- Mh[U]

  sum(nh) - n
}

philambda <- function(lambda, n, Nh, Sh, mh = NULL, Mh = NULL) {
  H <- seq_along(Nh)
  dh <- Sh * Nh

  L <- which((dh^2) / (mh^2) <= lambda)
  U <- which((dh^2) / (Mh^2) >= lambda)
  Hc <- H[-c(L, U)]

  # cat("L: ",L," U: ",U,"\n")
  (sum(dh[Hc]) / (n - sum(mh[L]) - sum(Mh[U])))^2 # lambda_n
}
