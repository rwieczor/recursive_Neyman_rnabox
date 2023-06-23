# Example that shows that fixed point iteration may not converge.

fpia <- stratallo:::fpia
SimpleGreedy <- stratallo:::SimpleGreedy
CapacityScaling <- stratallo:::CapacityScaling
glambda <- stratallo:::glambda
philambda <- stratallo:::philambda

# Generate some artificial parameters.
seed <- 308
set.seed(seed)
cat("SEED = ", seed, "\n")
Nh <- rep(100, 4)
Sh <- round(abs(rcauchy(4, 1)), 1)
Sh[Sh == 0] <- 1
mh <- rep(10, 4)
Mh <- rep(50, 4)
dh <- Sh * Nh
N <- sum(Nh)
(s1 <- sum(mh))
(s2 <- sum(Mh))
(n <- round(0.2 * N))

# Valid optimal allocation: (13, 10,  10,  47).
(alloc0 <- CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh))

(lam0 <- (sum(Sh * Nh)^2) / (n^2))

# Default starting value lambda0 - wrong allocation (no convergence).
(al_test <- fpia(n, Nh, Sh, mh, Mh, lambda0 = lam0, maxiter = 100)$nh)

# Values causing oscillations of the algorithm.
lam1 <- 1444
lam2 <- 739.84

# Starting value, better approximated.
# lam<-uniroot(function(x) glambda(x,n,Nh,Sh,mh,Mh),lower = lam0/1e8, upper=lam0*1e8)$root
lam <- uniroot(function(x) glambda(x, n, Nh, Sh, mh, Mh),
  lower = lam0 / 10, upper = lam0 * 10,
  extendInt = "yes"
)$root
lam

# Initial interval for bisection, according to JW.
rm <- dh / mh
rM <- dh / Mh
r <- c(rm, rM)
a <- min(r)^2
b <- max(r)^2

lam <- uniroot(function(x) glambda(x, n, Nh, Sh, mh, Mh), lower = a, upper = b)$root
lam # lambda_{*}

# Correct optimal allocation.
(al_test <- fpia(n, Nh, Sh, mh, Mh, lambda0 = lam)$nh)

# Pictures of g(lambda) function, the root of which we search.
fg <- function(x) {
  glambda(x, n, Nh, Sh, mh, Mh)
}

curve(Vectorize(fg)(x), 1e-6, 1700, ylim = c(-15, 15))
abline(h = 0, col = 2)
points(lam, fg(lam), cex = 1.5, pch = 19, col = 3)
points(lam0, fg(lam0), cex = 1.5, pch = 19, col = 2) # wrong starting lambda
points(lam1, fg(lam1), cex = 1.5, pch = 19, col = 4)
points(lam2, fg(lam2), cex = 1.5, pch = 19, col = 4)

library(ggplot2)
library(latex2exp)
library(patchwork)

df0 <- data.frame(x = lam0, y = fg(lam0))
dfstar <- data.frame(x = lam, y = fg(lam))
df1 <- data.frame(x = lam1, y = fg(lam1))
df2 <- data.frame(x = lam2, y = fg(lam2))

p1 <-
  ggplot(data.frame(x = seq(500, 1700, by = 0.1)), aes(x)) +
  geom_function(fun = Vectorize(fg), n = 20000, linewidth = 1) +
  geom_hline(yintercept = 0, color = "red") +
  geom_point(data = df0, aes(x = x, y = y), size = 4) +
  geom_linerange(data = df0, aes(ymin = 0, ymax = y), linetype = "dotted") +
  geom_linerange(data = df1, aes(ymin = 0, ymax = y), linetype = "dotted") +
  geom_linerange(data = df2, aes(ymin = 0, ymax = y), linetype = "dotted") +
  geom_point(data = dfstar, aes(x = x, y = y), size = 4) +
  geom_point(data = df1, aes(x = x, y = y), size = 4, shape = 21) +
  geom_point(data = df2, aes(x = x, y = y), size = 4, shape = 21) +
  geom_linerange(data = dfstar, aes(ymin = 0, ymax = y), linetype = "dotted") +
  ylim(-15, 10) +
  labs(x = TeX("$\\lambda$"), y = TeX("$g(\\lambda)$")) +
  # title = TeX("Example of function $g(\\lambda)$")) +
  theme_bw(base_size = 12)

p1 <- p1 +
  annotate("text", x = lam0 - 10, y = -2, label = TeX("$\\lambda_0$"), size = 5) +
  annotate("text", x = lam, y = -2, label = TeX("$\\lambda^{*}$"), size = 5) +
  annotate("text", x = lam1, y = -2, label = TeX("$\\lambda_{1}$"), size = 5) +
  annotate("text", x = lam2, y = -2, label = TeX("$\\lambda_{2}$"), size = 5)

fphi <- function(x) {
  philambda(x, n, Nh, Sh, mh, Mh)
}

curve(Vectorize(fphi)(x), 500, 1700)
abline(a = 0, b = 1, col = 2)
points(lam, fphi(lam), cex = 1.5, pch = 19, col = 3)
points(lam0, fphi(lam0), cex = 1.5, pch = 19, col = 2) # wrong starting lambda

x <- seq(1 - 6, 1700, by = 10)
y <- Vectorize(fphi)(x)
x1 <- min(x[which(is.infinite(y))])
x2 <- max(x[which(is.infinite(y))])
df <- data.frame(x = x, y = y)
dfa <- df[df$x >= x1 - 0.05 & df$x <= x2 + 0.05, ]
df0 <- data.frame(x = lam0, y = fphi(lam0), miny = min(y))
df1 <- data.frame(x = lam1, y = fphi(lam1), miny = min(y))
df2 <- data.frame(x = lam2, y = fphi(lam2), miny = min(y))
dfstar <- data.frame(x = lam, y = fphi(lam), miny = min(y))

p2 <-
  ggplot(df, aes(x = x, y = y)) +
  geom_function(fun = Vectorize(fphi), n = 20000, na.rm = TRUE, linewidth = 1) +
  # geom_line(aes(x = x,y = y)) +
  geom_line(data = dfa, aes(x = x, y = y), color = "white", linewidth = 1) +
  # xlim(x1-100,x2+100) + ylim(5000,25000) +
  geom_abline(intercept = 1e-6, slope = 1, color = "red") +
  geom_linerange(data = df0, aes(ymin = miny, ymax = y), linetype = "dotted") +
  geom_point(data = dfstar, aes(x = x, y = y), size = 4) +
  geom_point(data = df0, aes(x = x, y = y), size = 4) +
  geom_point(data = df1, aes(x = x, y = y), size = 4, shape = 21) +
  geom_point(data = df2, aes(x = x, y = y), size = 4, shape = 21) +
  geom_linerange(data = dfstar, aes(ymin = miny, ymax = y), linetype = "dotted") +
  geom_linerange(data = df1, aes(ymin = miny, ymax = y), linetype = "dotted") +
  geom_linerange(data = df2, aes(ymin = miny, ymax = y), linetype = "dotted") +
  # coord_cartesian(xlim=c(0,20000),ylim=c(0,20000)) +
  labs(x = TeX("$\\lambda$"), y = TeX("$\\phi(\\lambda)$")) +
  # title = TeX("Example of function $\\phi(\\lambda)$")) +
  theme_bw(base_size = 12)

p2 <- p2 +
  annotate("text", x = lam0, y = 0, label = TeX("$\\lambda_0$"), size = 5) +
  annotate("text", x = lam, y = 0, label = TeX("$\\lambda^{*}$"), size = 5) +
  annotate("text", x = lam1, y = 40, label = TeX("$\\lambda_{1}$"), size = 5) +
  annotate("text", x = lam2, y = 40, label = TeX("$\\lambda_{2}$"), size = 5)

pp <- p1 + p2
# ggsave("fig_fixedpi.jpg",pp,device="jpg", dpi=1200, width = 8, height = 8/1.618)
ggsave("./figures/fig_fpia.pdf", pp, device = "pdf", dpi = 1200, width = 8, height = 8 / 1.618)
