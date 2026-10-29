# Single entry point. Run from the repository root: Rscript src/run.R
# Monte Carlo samples, population values and S&P 500 estimates are cached in CACHE_DIR.
suppressPackageStartupMessages(library(Rcpp))
source(file.path("src", "config.R"))
sourceCpp(file.path("src", "egarch.cpp"))
for (f in c("estimators.R", "simulation.R", "empirical.R", "contraction.R", "plots.R"))
  source(file.path("src", f))
for (dir in c(CACHE_DIR, FIG_DIR)) dir.create(dir, showWarnings = FALSE)
cached <- function(name, compute) {
  path <- file.path(CACHE_DIR, name)
  if (!file.exists(path)) saveRDS(compute(), path)
  readRDS(path)
}
section <- function(title) cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
put <- function(x, d) {  # d decimals, "-0.0" printed as "0.0"
  print(noquote(sub("^-(0\\.0+)$", "\\1", formatC(x, format = "f", digits = d))))
}

# Monte Carlo study (Section 6) -----------------------------------------------
for (dn in names(DESIGNS))
  if (!all(file.exists(sapply(T_GRID, mc_file, dn = dn)))) simulate_design(dn)
pop <- cached("population.rds", function() lapply(DESIGNS, function(d) population(d$par, d$v)))
curve <- cached("pop_curve.rds", pop_curve)

section("Section 6: discarded samples")
print(discarded())
cat("smallest margin of the filter start over (omega - alpha c) / (1 - beta):",
    round(start_margin(), 2), "\n")
cat(paste0("kappa4 of GED ", DESIGNS$d5$v, ":"), round(law_k(DESIGNS$d5$v), 2), "\n")
cat("Corollary 5.5 thresholds in Design 2:",
    round(thresholds(DESIGNS$d2$par, DESIGNS$d2$v), 2), "\n")

section("Table 1: bias, Monte Carlo SD and asymptotic SD, Design 1")
put(table_bias(pop), 4)

section("Tables 2 and E.1: variance ratios over J")
tabs <- table_ratios(pop)
for (dn in names(tabs)) {  # OLS ratios at two decimals, standard errors at three
  x <- tabs[[dn]]
  out <- formatC(x, format = "f", digits = 3)
  i <- !grepl("^se", rownames(x))
  out[i, c("LS_nu", "LS_phi")] <- formatC(x[i, c("LS_nu", "LS_phi")], format = "f", digits = 2)
  cat("\nDesign", match(dn, names(DESIGNS)), "\n")
  print(noquote(out))
}
tl <- paste("T", T_GRID)
z <- sapply(tabs, function(x)
  max(abs(x[tl[2], ] - x["population", ]) / x[paste("se", T_GRID[2]), ]))
w <- sapply(tabs, function(x) x["population", "omega"] - x[tl[1], "omega"])
cat(sprintf("\nlargest |simulated - population| / s.e. at T = %d: %.1f\n", T_GRID[2], max(z)))
cat(sprintf("largest shortfall of the simulated omega ratio at T = %d: %.2f\n", T_GRID[1], max(w)))

section("Figure 1: population 2S/J ratio for the mean log variance")
m <- matrix(curve$ratio, length(CURVE_GAMMA), dimnames = list(CURVE_GAMMA, CURVE_BETA))
put(m, 2)
plot_curve(curve)

section(paste0("Remark 2.2: E log Lambda^(", INV_R, ") at the designs"))
put(sapply(DESIGNS, function(d) elog_lip(d$par, d$v)), 2)

# Application to S&P 500 returns (Section 7) ----------------------------------
d <- sp500()
E <- cached("empirical.rds", function() fit_sp500(d$logret_pct))

section("Table 3: descriptive statistics")
put(rbind(returns = describe(E$y), `J residuals` = describe(E$z)), 3)

section("Table 4: estimates and standard errors")
S <- sandwich(E)
put(E$tab, 4)
put(S$se, 4)
cat("GED shape", round(E$vhat, 2), " kappa4", round(ged_kappa(E$vhat), 2), "\n")
th <- E$tab["J", ]
cat(paste0("E log Lambda^(", INV_R, ") at the J estimate, Gaussian and GED:"),
    round(c(elog_lip(th, 2), elog_lip(shift_omega(th, E$vhat), E$vhat)), 2), "\n")

section("Table 5 and Section 7: population ratios at the J estimate")
r <- plugin_ratios(E)
put(r[, colnames(r) != "omega_k"], 3)
cat("omega with nu known:", round(r["Gaussian", "omega_k"], 4), "\n")

section("Section 7: sample variance ratios of the mean, 95% block bootstrap interval")
put(S$ratios, 2)
section("Section 7: Hausman z-statistics")
put(S$hausman, 1)
section("Section 7: ratios under the resampled and symmetrized residuals")
put(skew_ratios(E), 2)
cat("3S/J for phi implied by Corollary 5.5 at the residual kurtosis:",
    round(kurtosis_ratio(th, describe(E$z)[["kurtosis"]]), 2), "\n")

# Notes -----------------------------------------------------------------------
section(paste("Notes: largest absolute autocorrelation of returns, lags 1 to", ACF_LAGS))
put(plot_facts(d), 3)
