# Every parameter of the project, declared once.

# Paths -----------------------------------------------------------------------
DATA_FILE <- file.path("data", "SPX_daily.csv")
CACHE_DIR <- "cache"
FIG_DIR   <- "latex"

# Monte Carlo designs (Section 6) ---------------------------------------------
# Design 1 is the base design, and Designs 2 to 5 change one factor each.
# v is the GED shape of the innovations, with v = 2 the Gaussian law.
BASE <- list(nu = 0, phi = 0.05, omega = 0, alpha = 0.15, gamma = -0.08, beta = 0.98)
DESIGNS <- list(
  d1 = list(par = BASE, v = 2),
  d2 = list(par = modifyList(BASE, list(gamma = 0)), v = 2),
  d3 = list(par = modifyList(BASE, list(gamma = -0.10)), v = 2),
  d4 = list(par = modifyList(BASE, list(beta = 0.95)), v = 2),
  d5 = list(par = BASE, v = 1.5))
T_GRID    <- c(1000, 10000)
N_REP     <- 1000
BURN      <- 1000
SEED      <- 2026
CORES     <- if (.Platform$OS.type == "windows") 1 else 2  # mclapply forks on Unix only
BOUNDARY  <- 0.99  # |gamma| / alpha above this counts as on the boundary of Theta
BOOT_REPS <- 500  # bootstrap of the Monte Carlo variance ratios

# Quasi-likelihood maximization -----------------------------------------------
MAXIT    <- 1000
RELTOL   <- 1e-12
RESTARTS <- 10

# Population values -----------------------------------------------------------
N_POP       <- 2e7  # path length of the derivative recursion
CURVE_GAMMA <- seq(0, -0.14, by = -0.01)  # Figure 1
CURVE_BETA  <- c(0.90, 0.95, 0.97, 0.98, 0.99)

# Contraction E log Lambda^(r) (Remark 2.2) -----------------------------------
INV_R     <- 20
INV_N     <- 2e5   # simulated path length
INV_BURN  <- 2000  # first admissible block start
INV_DRAWS <- 20000 # blocks averaged
INV_GRID  <- seq(0, 25, length.out = 300)  # offsets above (omega - alpha c) / (1 - beta)

# Application to S&P 500 returns (Section 7) ----------------------------------
EMP_WINDOW <- c("1970-01-02", "2026-08-31")
BLOCK_LEN  <- 50    # moving-block bootstrap of the mean variance ratios
BLOCK_REPS <- 1999
ACF_LAGS   <- 100   # figure of the notes
