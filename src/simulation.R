# Monte Carlo study and population values (Section 6, Figure 1, Appendix E).

mc_file <- function(dn, Tn) file.path(CACHE_DIR, sprintf("mc_%s_T%d.rds", dn, Tn))

# Sample r of design dn at size index i, with its own seed.
mc_sample <- function(dn, i, r) {
  d <- DESIGNS[[dn]]
  set.seed(SEED + 1e5 * (10 * match(dn, names(DESIGNS)) + i) + r)
  p <- d$par
  sim_ar_egarch(rlaw(T_GRID[i] + BURN, d$v), BURN, p$nu, p$phi, p$omega, p$alpha, p$gamma,
                p$beta, law_c(d$v))
}

# LS, 2S, 3S and J on common samples, nu estimated (_e) and nu known (_k).
simulate_design <- function(dn) {
  nm <- c(paste0("LS_", EST[1:2]), paste0("S2_", EST), paste0("S3_", EST), paste0("J_", EST),
          "conv_S2", "conv_J")
  nm <- c(paste0(nm, "_e"), paste0(nm, "_k"))
  c0 <- law_c(DESIGNS[[dn]]$v)
  for (i in seq_along(T_GRID)) {
    M <- do.call(rbind, parallel::mclapply(1:N_REP, function(r) {
      y <- mc_sample(dn, i, r)
      tryCatch({
        a <- estimate_all(y, c0)
        b <- estimate_all(y, c0, nu0 = DESIGNS[[dn]]$par$nu)
        c(a$LS, a$S2, a$S3, a$J, a$conv, b$LS, b$S2, b$S3, b$J, b$conv)
      }, error = function(e) rep(NA_real_, length(nm)))
    }, mc.cores = CORES))
    colnames(M) <- nm
    saveRDS(list(est = M), mc_file(dn, T_GRID[i]))
  }
}

mc_all <- function(dn, Tn) readRDS(mc_file(dn, Tn))$est
converged <- function(M) complete.cases(M) & rowSums(M[, grep("^conv", colnames(M))] != 0) == 0

# Discarded samples, and those of them with an estimate on the boundary |gamma| = alpha.
discarded <- function() {
  t(sapply(names(DESIGNS), function(dn) {
    M <- mc_all(dn, T_GRID[1])
    bad <- !converged(M)
    q <- sapply(c("S2_%s_e", "J_%s_e", "S2_%s_k", "J_%s_k"), function(f)
      abs(M[, sprintf(f, "gamma")]) / M[, sprintf(f, "alpha")])
    setNames(c(sum(bad), sum(bad & apply(q, 1, max) > BOUNDARY, na.rm = TRUE),
               sum(!converged(mc_all(dn, T_GRID[2])))),
             c(paste("discarded, T =", T_GRID[1]), "on the boundary",
               paste("discarded, T =", T_GRID[2])))
  }))
}

# Smallest margin of the filter start over (omega - alpha c) / (1 - beta) at the 2S and J
# estimates of every sample.
start_margin <- function() {
  min(sapply(names(DESIGNS), function(dn) sapply(seq_along(T_GRID), function(i) {
    d <- DESIGNS[[dn]]
    c0 <- law_c(d$v)
    M <- mc_all(dn, T_GRID[i])
    h0 <- sapply(1:N_REP, function(r) {
      y <- mc_sample(dn, i, r)
      c(e = start_h0(y, ls_ar(y)), k = start_h0(y, ls_ar(y, nu0 = d$par$nu)))
    })
    lower <- function(s) (M[, s["omega"]] - M[, s["alpha"]] * c0) / (1 - M[, s["beta"]])
    cols <- function(est, nk) setNames(paste0(est, "_", EST, "_", nk), EST)
    min(h0["e", ] - lower(cols("S2", "e")), h0["e", ] - lower(cols("J", "e")),
        h0["k", ] - lower(cols("S2", "k")), h0["k", ] - lower(cols("J", "k")), na.rm = TRUE)
  })))
}

# Table 1: Design 1, nu estimated. Bias, Monte Carlo s.d. and asymptotic s.d. sqrt(V / T).
table_bias <- function(pop) {
  e <- pop$d1$est
  rows <- list(c("nu", "LS", e$VLS[1, 1]), c("nu", "S3", e$V3S_m[1, 1]),
               c("nu", "J", e$VJ[1, 1]), c("phi", "LS", e$VLS[2, 2]),
               c("phi", "S3", e$V3S_m[2, 2]), c("phi", "J", e$VJ[2, 2]))
  for (k in 1:4) rows <- c(rows, list(c(EST[k + 2], "S2", e$V2S_v[k, k]),
                                      c(EST[k + 2], "J", e$VJ[k + 2, k + 2])))
  X <- lapply(T_GRID, function(Tn) {
    M <- mc_all("d1", Tn)
    M[converged(M), ]
  })
  out <- t(sapply(rows, function(r) {
    x <- lapply(X, function(M) M[, sprintf("%s_%s_e", r[2], r[1])])
    unlist(lapply(seq_along(T_GRID), function(i) c(bias = mean(x[[i]]) - BASE[[r[1]]],
      sd = sd(x[[i]]), ad = sqrt(as.numeric(r[3]) / T_GRID[i]))))
  }))
  rownames(out) <- sapply(rows, function(r) paste(r[1], r[2]))
  colnames(out) <- paste(rep(c("bias", "sd", "ad"), 2), rep(T_GRID, each = 3))
  out
}

# Tables 2 and E.1: variance ratios over J, population and Monte Carlo with bootstrap s.e.
ratio_cols <- list(gamma = c("S2_gamma_e", "J_gamma_e"), gamma_k = c("S2_gamma_k", "J_gamma_k"),
                   omega = c("S2_omega_e", "J_omega_e"), omega_k = c("S2_omega_k", "J_omega_k"),
                   S3_nu = c("S3_nu_e", "J_nu_e"), S3_phi = c("S3_phi_e", "J_phi_e"),
                   LS_nu = c("LS_nu_e", "J_nu_e"), LS_phi = c("LS_phi_e", "J_phi_e"),
                   alpha = c("S2_alpha_e", "J_alpha_e"), beta = c("S2_beta_e", "J_beta_e"))
table_ratios <- function(pop) {
  boot_se <- function(a, b) {
    n <- length(a)
    sd(replicate(BOOT_REPS, {
      i <- sample.int(n, n, TRUE)
      var(a[i]) / var(b[i])
    }))
  }
  set.seed(SEED)
  lapply(setNames(names(DESIGNS), names(DESIGNS)), function(dn) {
    d <- DESIGNS[[dn]]
    pv <- var_ratios(pop[[dn]], d$par)[names(ratio_cols)]
    if (d$par$gamma == 0) {  # closed forms for the intercept at gamma0 = 0
      cf <- nu_ratios(d$v, d$par$omega, d$par$alpha, d$par$beta)
      pv[c("LS_nu", "S3_nu")] <- cf[c("LS", "S3")]
    }
    out <- rbind(population = pv)
    for (Tn in T_GRID) {
      M <- mc_all(dn, Tn)
      M <- M[converged(M), ]
      v <- sapply(ratio_cols, function(z) c(var(M[, z[1]]) / var(M[, z[2]]),
                                            boot_se(M[, z[1]], M[, z[2]])))
      out <- rbind(out, v)
      rownames(out)[nrow(out) - 1:0] <- paste(c("T", "se"), Tn)
    }
    out
  })
}

# Figure 1: 2S/J for the mean log variance against gamma0, nu estimated, Gaussian law.
pop_curve <- function() {
  curve <- expand.grid(gamma = CURVE_GAMMA, beta = CURVE_BETA)
  curve$ratio <- mapply(function(g, b) {
    d <- modifyList(BASE, list(gamma = g, beta = b))
    var_ratios(population(d, 2), d)[["level"]]
  }, curve$gamma, curve$beta)
  curve
}
