# Innovation laws, the LS, 2S, 3S and J estimators, and the population covariances.

# Innovation laws -------------------------------------------------------------
# Standardized GED(v): density proportional to exp(-|x / lam|^v / 2), variance one.
ged_lam <- function(v) sqrt(2^(-2 / v) * gamma(1 / v) / gamma(3 / v))
ged_abs <- function(v) ged_lam(v) * 2^(1 / v) * gamma(2 / v) / gamma(1 / v)  # E|eta|
ged_kappa <- function(v) gamma(5 / v) * gamma(1 / v) / gamma(3 / v)^2  # E eta^4
dged <- function(x, v) {
  l <- ged_lam(v)
  v * exp(-0.5 * abs(x / l)^v) / (l * 2^(1 + 1 / v) * gamma(1 / v))
}
dlaw <- function(x, v) if (v == 2) dnorm(x) else dged(x, v)
law_c <- function(v) if (v == 2) sqrt(2 / pi) else ged_abs(v)
law_k <- function(v) if (v == 2) 3 else ged_kappa(v)
rlaw <- function(n, v) {
  if (v == 2) return(rnorm(n))
  ged_lam(v) * (2 * rgamma(n, 1 / v))^(1 / v) * sample(c(-1, 1), n, replace = TRUE)
}

# Estimators (Definition 2.1) -------------------------------------------------
EST <- c("nu", "phi", "omega", "alpha", "gamma", "beta")

# Weighted LS of y on (1, y lagged), or on y lagged when nu is known, returning (nu, phi).
ls_ar <- function(y, w = NULL, nu0 = NULL) {
  n <- length(y)
  Y <- y[2:n]
  X <- y[1:(n - 1)]
  if (is.null(w)) w <- rep(1, n - 1)
  if (is.null(nu0)) {
    b <- lm.wfit(cbind(1, X), Y, w)$coefficients
    c(b[1] / (1 - b[2]), b[2])
  } else {
    c(nu0, lm.wfit(as.matrix(X - nu0), Y - nu0, w)$coefficients)
  }
}

# phi = tanh(f), alpha = exp(a), gamma = alpha tanh(g) and beta = plogis(b).
to_u <- function(th) {
  c(th[1], atanh(th[2]), th[3], log(th[4]), atanh(th[5] / th[4]), qlogis(th[6]))
}
from_u <- function(u) {
  al <- exp(u[4])
  c(u[1], tanh(u[2]), u[3], al, al * tanh(u[5]), plogis(u[6]))
}
jac_u <- function(u) {  # d theta / d u, diagonal except the gamma row
  J <- diag(6)
  al <- exp(u[4])
  J[2, 2] <- 1 - tanh(u[2])^2
  J[4, 4] <- al
  J[5, 4] <- al * tanh(u[5])
  J[5, 5] <- al * (1 - tanh(u[5])^2)
  J[6, 6] <- plogis(u[6]) * (1 - plogis(u[6]))
  J
}

# Maximizes the average quasi-likelihood over the free coordinates of theta.
fit_qml <- function(th0, y, c, h0, free, nu_known) {
  u0 <- to_u(th0)
  qll_u <- function(uf) {
    u <- u0
    u[free] <- uf
    list(u = u, r = qll(from_u(u), y, 1, c, h0, nu_known))
  }
  f <- function(uf) {
    r <- qll_u(uf)$r
    if (!is.finite(r$value)) 1e10 else -r$value
  }
  gr <- function(uf) {
    q <- qll_u(uf)
    if (!is.finite(q$r$value)) return(rep(0, length(uf)))
    -(crossprod(jac_u(q$u), q$r$grad))[free]
  }
  ctrl <- list(maxit = MAXIT, reltol = RELTOL)
  o <- optim(u0[free], f, gr, method = "BFGS", control = ctrl)
  for (i in 1:RESTARTS) {  # restart from the last point while the iteration limit binds
    if (o$convergence != 1) break
    o <- optim(o$par, f, gr, method = "BFGS", control = ctrl)
  }
  u <- u0
  u[free] <- o$par
  list(theta = from_u(u), conv = o$convergence)
}

# Filter start: log mean squared LS residual.
start_h0 <- function(y, m1) {
  n <- length(y)
  log(mean((y[2:n] - m1[1] - m1[2] * (y[1:(n - 1)] - m1[1]))^2))
}

# LS, 2S, 3S and J for the AR(1) mean. Variance start (0.1 h0, 0.1, 0, 0.9), J starts at 2S.
estimate_all <- function(y, c, nu0 = NULL) {
  nu_known <- !is.null(nu0)
  m1 <- ls_ar(y, nu0 = nu0)
  h0 <- start_h0(y, m1)
  f2 <- fit_qml(c(m1, 0.1 * h0, 0.1, 0, 0.9), y, c, h0, free = 3:6, nu_known = TRUE)
  th2 <- f2$theta
  m3 <- ls_ar(y, w = exp(-filter_h(th2, y, 1, c, h0)), nu0 = nu0)
  fJ <- fit_qml(th2, y, c, h0, free = if (nu_known) 2:6 else 1:6, nu_known = nu_known)
  list(LS = m1, S2 = th2, S3 = c(m3, th2[3:6]), J = fJ$theta,
       conv = c(S2 = f2$conv, J = fJ$conv), h0 = h0)
}

# Population covariances (Sections 3 and 4) ----------------------------------
# Covariances of J, LS, 2S and 3S from the population moments pm.
pop_cov <- function(pm, kappa4, nu_known = FALSE) {
  k <- 6
  m <- 2
  C <- pm$C
  W <- pm$W
  M <- pm$M
  S <- pm$S
  if (nu_known) {
    C <- C[2:k, 2:k]
    W <- W[-1, -1, drop = FALSE]
    M <- M[-1, -1, drop = FALSE]
    S <- S[-1, -1, drop = FALSE]
    m <- 1
    k <- 5
  }
  DW <- matrix(0, k, k)
  DW[1:m, 1:m] <- W
  I <- (kappa4 - 1) / 4 * C + DW
  J <- C / 2 + DW
  Ji <- solve(J)
  VLS <- solve(M) %*% S %*% solve(M)
  iv <- (m + 1):k
  Jvm <- J[iv, 1:m]
  V2S <- solve(J[iv, iv]) %*% (I[iv, iv] + Jvm %*% VLS %*% t(Jvm)) %*% solve(J[iv, iv])
  list(VJ = Ji %*% I %*% Ji, VLS = VLS, V2S_v = V2S, V3S_m = solve(W))
}

# Covariances with nu estimated (est) and known (kn) at th, from a path of N_POP draws.
population <- function(th, v) {
  th <- unlist(th)
  set.seed(SEED)
  eta <- rlaw(N_POP + BURN, v)
  pm <- pop_moments(eta, BURN, th[["nu"]], th[["phi"]], th[["omega"]], th[["alpha"]],
                    th[["gamma"]], th[["beta"]], law_c(v))
  list(est = pop_cov(pm, law_k(v)), kn = pop_cov(pm, law_k(v), nu_known = TRUE))
}

# Variance ratios over J; level is the mean log variance omega / (1 - beta).
var_ratios <- function(p, th) {
  e <- p$est
  k <- p$kn
  g <- c(1, 0, 0, th[["omega"]] / (1 - th[["beta"]])) / (1 - th[["beta"]])
  c(level = drop(t(g) %*% e$V2S_v %*% g) / drop(t(g) %*% e$VJ[3:6, 3:6] %*% g),
    omega = e$V2S_v[1, 1] / e$VJ[3, 3], alpha = e$V2S_v[2, 2] / e$VJ[4, 4],
    gamma = e$V2S_v[3, 3] / e$VJ[5, 5], beta = e$V2S_v[4, 4] / e$VJ[6, 6],
    omega_k = k$V2S_v[1, 1] / k$VJ[2, 2], gamma_k = k$V2S_v[3, 3] / k$VJ[4, 4],
    S3_nu = e$V3S_m[1, 1] / e$VJ[1, 1], S3_phi = e$V3S_m[2, 2] / e$VJ[2, 2],
    LS_nu = e$VLS[1, 1] / e$VJ[1, 1], LS_phi = e$VLS[2, 2] / e$VJ[2, 2])
}

# Closed forms at gamma0 = 0 (Theorem 3.2) ------------------------------------
# Kurtosis thresholds 5 + tau and 5 + tau (1 - varpi) of Corollary 5.5.
thresholds <- function(d, v) {
  Elaw <- function(f) integrate(function(x) f(x) * dlaw(x, v), -Inf, Inf, rel.tol = 1e-10)$value
  c0 <- law_c(v)
  g <- function(x) d$alpha * (abs(x) - c0) + d$gamma * x
  Esig <- function(s) exp(s * d$omega / (2 * (1 - d$beta))) *
    prod(sapply(0:400, function(j) Elaw(function(x) exp(s * d$beta^j * g(x) / 2))))
  rho <- Elaw(function(x) (d$beta - 0.5 * (d$alpha * abs(x) + d$gamma * x))^2)
  tau <- d$alpha^2 / (1 - rho)
  varpi <- c0^2 * Esig(-1)^2 / Esig(-2)
  c(phi = 5 + tau, nu = 5 + tau * (1 - varpi))
}

# LS/J = E sigma^2 E sigma^-2 / R(kappa4, x) and 3S/J = 1 / R(kappa4, x) for nu, with
# x = tau (1 - varpi), under GED(v) (Corollary 5.2).
nu_ratios <- function(v, om, al, be) {
  Ex <- function(g) {
    2 * integrate(function(x) g(x) * dged(x, v), 0, 60, rel.tol = 1e-12,
                  subdivisions = 1000L)$value
  }
  c <- Ex(abs)
  k4 <- Ex(function(x) x^4)
  Es <- function(s) exp(s * om / (1 - be) + sum(sapply(1:400, function(j)
    log(Ex(function(x) exp(s * be^(j - 1) * al * (abs(x) - c)))))))
  rho <- Ex(function(x) (be - 0.5 * al * abs(x))^2)
  tau <- al^2 / (1 - rho)
  x <- tau * (1 - c^2 * Es(-0.5)^2 / Es(-1))
  R <- (1 + (k4 - 1) * x / 4) / (1 + x / 2)^2
  c(LS = Es(1) * Es(-1) / R, S3 = 1 / R)
}
