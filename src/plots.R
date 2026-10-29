# Figure 1 of the report and the figure of the notes, written to FIG_DIR.

# Figure 1: population 2S/J variance ratio for the mean log variance against gamma0.
plot_curve <- function(curve) {
  pdf(file.path(FIG_DIR, "fig_level.pdf"), width = 5.5, height = 3.6, pointsize = 10)
  par(mar = c(4, 4.2, 0.6, 0.6), las = 1)
  plot(NA, xlim = range(CURVE_GAMMA), ylim = c(1, max(curve$ratio) + 0.05),
       xlab = expression(gamma[0]), ylab = "Variance ratio, 2S over J")
  lt <- c(3, 2, 4, 1, 6)
  lw <- c(1, 1, 1, 1, 1.6)
  for (i in seq_along(CURVE_BETA)) {
    s <- curve[curve$beta == CURVE_BETA[i], ]
    lines(s$gamma, s$ratio, lty = lt[i], lwd = lw[i])
  }
  for (i in 1:4) {  # Designs 1 to 4
    d <- DESIGNS[[i]]$par
    y <- curve$ratio[abs(curve$gamma - d$gamma) < 1e-9 & curve$beta == d$beta]
    points(d$gamma, y, pch = 19, cex = 0.7)
    text(d$gamma, y, i, pos = if (i == 4) 1 else 3, cex = 0.85)
  }
  legend("topright", legend = as.expression(lapply(rev(CURVE_BETA), function(b)
    bquote(beta[0] == .(sprintf("%.2f", b))))), lty = rev(lt), lwd = rev(lw), bty = "n")
  invisible(dev.off())
}

# Figure of the notes: returns, quantiles and autocorrelation functions.
plot_facts <- function(d) {
  r <- d$logret_pct
  dates <- as.Date(d$date)
  oi <- c(black = "#000000", blue = "#0072B2", red = "#D55E00")
  yr <- seq(as.Date(paste0(substr(EMP_WINDOW[1], 1, 4), "-01-01")), max(dates), by = "5 years")
  band <- 1.96 / sqrt(length(r))
  pdf(file.path(FIG_DIR, "fig_facts.pdf"), width = 7.5, height = 6.2)
  layout(matrix(1:4, 2, byrow = TRUE))
  par(family = "serif", mar = c(4, 4.2, 1.5, 1), cex.lab = 1.05, cex.axis = 0.95, las = 1)
  plot(dates, r, type = "l", lwd = 0.4, col = oi["black"], xlab = "Year",
       ylab = "Daily return (%)", xaxt = "n", main = "(a) Daily returns")
  axis(1, at = yr[seq(1, length(yr), 2)], labels = format(yr[seq(1, length(yr), 2)], "%Y"))
  q <- qqnorm((r - mean(r)) / sd(r), plot.it = FALSE)
  plot(q$x, q$y, pch = 16, cex = 0.3, xlab = "Normal quantile",
       ylab = "Standardized return quantile", main = "(b) Quantiles against the normal law")
  abline(0, 1, lty = 2, col = oi["red"])
  acf_panel <- function(x, main) {
    a <- acf(x, lag.max = ACF_LAGS, plot = FALSE)$acf[-1]
    plot(1:ACF_LAGS, a, type = "h", lwd = 1, col = oi["blue"], ylim = range(c(a, -band, band, 0)),
         xlab = "Lag (days)", ylab = "Autocorrelation", main = main)
    abline(h = 0)
    abline(h = c(-1, 1) * band, lty = 2, col = oi["red"])
    a
  }
  a <- acf_panel(r, "(c) ACF of returns")
  acf_panel(abs(r), "(d) ACF of absolute returns")
  invisible(dev.off())
  max(abs(a))
}
