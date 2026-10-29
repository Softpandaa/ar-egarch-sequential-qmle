# Application to S&P 500 returns (Section 7).

# Daily log returns in percent over EMP_WINDOW.
sp500 <- function() {
  d <- read.csv(DATA_FILE)
  d[d$date >= EMP_WINDOW[1] & d$date <= EMP_WINDOW[2] & !is.na(d$logret_pct), ]
}

# LS, 2S, 3S and J estimates, standardized J residuals and the GED shape fitted to them.
fit_sp500 <- function(y) {
  n <- length(y)
  c0 <- law_c(2)
  est <- estimate_all(y, c0)
  tab <- rbind(LS = c(est$LS, rep(NA, 4)), `2S` = est$S2, `3S` = est$S3, J = est$J)
  colnames(tab) <- EST
  th <- est$J
  e <- y[2:n] - th[1] - th[2] * (y[1:(n - 1)] - th[1])
  z <- e * exp(-filter_h(th, y, 1, c0, est$h0) / 2)
  nll <- function(v) {
    l <- ged_lam(v)
    -sum(log(v) - 0.5 * abs(z / l)^v - log(l * 2^(1 + 1 / v) * gamma(1 / v)))
  }
  list(tab = tab, vhat = optimize(nll, c(0.5, 4))$minimum, z = z, y = y)
}

# Table 3: mean, standard deviation, skewness and raw kurtosis.
describe <- function(x) {
  u <- (x - mean(x)) / sd(x)
  c(mean = mean(x), sd = sd(x), skewness = mean(u^3), kurtosis = mean(u^4))
}

# Residual eps, regressor xi = d eps / d m, b = xi / sigma, derivative series a of the
# filter and score s at theta, started at h0.
series <- function(th, y, h0) {
  c0 <- law_c(2)
  N <- length(y) - 1
  a <- matrix(0, N, 6)
  b <- xi <- matrix(0, N, 2)
  eta <- eps <- h <- numeric(N)
  at <- rep(0, 6)
  ht <- h0
  for (t in 1:N) {
    e <- y[t + 1] - th[1] - th[2] * (y[t] - th[1])
    x <- c(-(1 - th[2]), -(y[t] - th[1]))
    z <- e * exp(-ht / 2)
    a[t, ] <- at
    b[t, ] <- x * exp(-ht / 2)
    xi[t, ] <- x
    eta[t] <- z
    eps[t] <- e
    h[t] <- ht
    an <- (th[6] - 0.5 * (th[4] * abs(z) + th[5] * z)) * at
    an[1:2] <- an[1:2] + (th[4] * sign(z) + th[5]) * exp(-ht / 2) * x
    an[3:6] <- an[3:6] + c(1, abs(z) - c0, z, ht)
    ht <- th[3] + th[4] * (abs(z) - c0) + th[5] * z + th[6] * ht
    at <- an
  }
  s <- -0.5 * (1 - eta^2) * a
  s[, 1:2] <- s[, 1:2] - eta * b
  list(a = a, b = b, xi = xi, eps = eps, h = h, s = s)
}

# Table 4 standard errors, the 3S/J and OLS/J variance ratios of the mean with moving-block
# bootstrap intervals, and Hausman z-statistics. J is estimated by mean(a a') / 2 + diag(W)
# (Proposition 3.1), OLS is heteroskedasticity robust, 2S includes the first-step correction
# of (4.6), and 3S is weighted LS with the 2S weights.
sandwich <- function(E) {
  tab <- E$tab
  y <- E$y
  N <- length(y) - 1
  m1 <- tab["LS", 1:2]
  h0 <- start_h0(y, m1)
  iv <- 3:6
  Jhat <- function(q) {
    J <- crossprod(q$a) / N / 2
    J[1:2, 1:2] <- J[1:2, 1:2] + crossprod(q$b) / N
    J
  }
  se <- function(V) sqrt(diag(V) / N)
  qJ <- series(tab["J", ], y, h0)
  Ji <- solve(Jhat(qJ))
  VJ <- Ji %*% (crossprod(qJ$s) / N) %*% Ji
  psiJ <- (qJ$s %*% t(Ji))[, 1:2]
  qO <- series(c(m1, tab["2S", iv]), y, h0)
  psiO <- -(qO$xi * qO$eps) %*% solve(crossprod(qO$xi) / N)
  VLS <- crossprod(psiO) / N
  q2 <- series(tab["2S", ], y, h0)
  J2 <- Jhat(q2)
  Jvvi <- solve(J2[iv, iv])
  Jvm <- J2[iv, 1:2]
  Esp <- crossprod(q2$s[, iv], psiO) / N
  V2 <- Jvvi %*% (crossprod(q2$s[, iv]) / N - Jvm %*% t(Esp) - Esp %*% t(Jvm) +
                    Jvm %*% VLS %*% t(Jvm)) %*% Jvvi
  q3 <- series(c(tab["3S", 1:2], tab["2S", iv]), y, h0)
  w <- exp(-q2$h)
  psi3 <- -(q3$xi * q3$eps * w) %*% solve(crossprod(q3$xi * sqrt(w)) / N)
  V3 <- crossprod(psi3) / N
  ses <- rbind(LS = c(se(VLS), rep(NA, 4)), `2S` = c(se(VLS), se(V2)), `3S` = c(se(V3), se(V2)),
               J = se(VJ))
  colnames(ses) <- EST
  ratios <- function(i) {
    vJ <- colMeans(psiJ[i, ]^2)
    c(S3_J = colMeans(psi3[i, ]^2) / vJ, LS_J = colMeans(psiO[i, ]^2) / vJ)
  }
  set.seed(SEED)
  boot <- replicate(BLOCK_REPS, {
    st <- sample(1:(N - BLOCK_LEN + 1), ceiling(N / BLOCK_LEN), replace = TRUE)
    ratios(as.vector(outer(0:(BLOCK_LEN - 1), st, "+"))[1:N])
  })
  z <- function(d, pa, pb) d / sqrt(colMeans((pa - pb)^2) / N)
  list(se = ses,
       ratios = cbind(estimate = ratios(1:N), t(apply(boot, 1, quantile, c(0.025, 0.975)))),
       hausman = rbind(LS_J = z(m1 - tab["J", 1:2], psiO, psiJ),
                       LS_3S = z(m1 - tab["3S", 1:2], psiO, psi3),
                       S3_J = z(tab["3S", 1:2] - tab["J", 1:2], psi3, psiJ)))
}

# J estimate with omega shifted by alpha (c - c0) for the law GED(v), which leaves the filter
# unchanged.
shift_omega <- function(th, v) {
  th["omega"] <- th["omega"] + th["alpha"] * (law_c(v) - law_c(2))
  th
}

# Table 5 and the Section 7 comparisons: population ratios at the J estimate under the
# Gaussian law, under the fitted GED and at gamma0 = 0, and the omega ratio with nu known.
plugin_ratios <- function(E) {
  th <- E$tab["J", ]
  thg <- shift_omega(th, E$vhat)
  th0 <- th
  th0["gamma"] <- 0
  sp <- population(th, 2)
  r <- rbind(Gaussian = var_ratios(sp, th), GED = var_ratios(population(thg, E$vhat), thg),
             `gamma0 = 0` = var_ratios(population(th0, 2), th0))
  r[, colnames(r) != "gamma_k"]
}

# Plug-in ratios with innovations resampled from the standardized J residuals, by the general
# sandwich formulas with the E eta^3 cross terms, and under their symmetrized law.
skew_ratios <- function(E) {
  th <- E$tab["J", ]
  z <- (E$z - mean(E$z)) / sd(E$z)
  ratios <- function(eta) {
    ce <- mean(abs(eta))
    om <- th[3] + th[4] * (ce - law_c(2))
    g <- pop_moments(eta, BURN, th[1], th[2], om, th[4], th[5], th[6], ce)
    iv <- 3:6
    DW <- matrix(0, 6, 6)
    DW[1:2, 1:2] <- g$W
    J <- g$C / 2 + DW
    Ji <- solve(J)
    VJ <- Ji %*% g$I %*% Ji
    Mi <- solve(g$M)
    VLS <- Mi %*% g$S %*% Mi
    Jvm <- J[iv, 1:2]
    Esp <- -g$G[iv, 1:2] %*% Mi  # E[s_v p'] with p = -M^-1 xi eps
    U <- g$I[iv, iv] - Jvm %*% t(Esp) - Esp %*% t(Jvm) + Jvm %*% VLS %*% t(Jvm)
    V2 <- solve(J[iv, iv]) %*% U %*% solve(J[iv, iv])
    c(setNames(diag(V2) / diag(VJ)[iv], EST[iv]), S3_phi = diag(solve(g$W))[2] / VJ[2, 2])
  }
  set.seed(SEED)
  emp <- ratios(sample(z, N_POP + BURN, replace = TRUE))
  zs <- c(z, -z)
  rbind(empirical = emp, symmetrized = ratios(sample(zs / sd(zs), N_POP + BURN, replace = TRUE)))
}

# 3S/J variance ratio for phi implied by Corollary 5.5 at gamma0 = 0 and kurtosis kappa,
# with tau at the J estimate under the Gaussian law.
kurtosis_ratio <- function(th, kappa) {
  be <- th[["beta"]]
  al <- th[["alpha"]]
  rho <- integrate(function(x) (be - 0.5 * al * abs(x))^2 * dnorm(x), -Inf, Inf)$value
  x <- al^2 / (1 - rho)
  (1 + x / 2)^2 / (1 + (kappa - 1) * x / 4)
}
