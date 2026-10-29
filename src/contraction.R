# E log Lambda^(r) of Remark 2.2, with Lambda^(r) the Lipschitz coefficient of r composed
# filter maps on [(omega - alpha c) / (1 - beta), Inf), approximated on a grid.
elog_lip <- function(th, v) {
  c <- law_c(v)
  set.seed(SEED)  # common draws, independent of the order of evaluation
  th <- unlist(th)
  om <- th[["omega"]]
  al <- th[["alpha"]]
  ga <- th[["gamma"]]
  be <- th[["beta"]]
  eta <- rlaw(INV_N, v)
  h <- numeric(INV_N)
  h[1] <- om / (1 - be)
  for (t in 2:INV_N) h[t] <- om + al * (abs(eta[t - 1]) - c) + ga * eta[t - 1] + be * h[t - 1]
  eps <- exp(h / 2) * eta
  x0 <- (om - al * c) / (1 - be) + INV_GRID
  mean(sapply(sample(INV_BURN:(INV_N - INV_R - 1), INV_DRAWS, TRUE), function(s) {
    x <- x0
    d <- rep(1, length(x0))
    for (e in eps[s:(s + INV_R - 1)]) {
      ex <- exp(-x / 2)
      d <- d * (be - 0.5 * (al * abs(e) + ga * e) * ex)
      x <- om + al * (abs(e) * ex - c) + ga * e * ex + be * x
    }
    log(max(abs(d)))
  }))
}
