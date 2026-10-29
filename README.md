# Efficiency Loss of Sequential Quasi-Maximum Likelihood Estimation in the AR(p)-EGARCH(1,1) Model

This study compares joint and sequential Gaussian quasi-maximum likelihood estimation of the AR(p)-EGARCH(1,1) model. The two-step (2S) estimator fits the mean by ordinary least squares and then the variance equation by quasi-maximum likelihood, and the three-step (3S) estimator refits the mean by feasible generalized least squares. The asymptotic covariances of the three estimators are derived from the information matrices, evaluated in a Monte Carlo study and applied to daily S&P 500 returns from January 1970 to August 2026. An accompanying note introduces AR-EGARCH models and quasi-maximum likelihood estimation.

## Findings

The information matrix is not block diagonal between the mean and the variance parameters, even without a leverage effect, through a single entry that links the intercept of the mean to the coefficient on the lagged innovation in the log variance. Without leverage, the 2S estimator therefore loses efficiency only for this coefficient, and only when the intercept is estimated. With leverage, the loss extends to the intercept and the persistence of the log variance at second order in the leverage coefficient. It is largest for the mean log variance, whose asymptotic variance exceeds that of the joint estimator by 15% in the base design of the Monte Carlo study and by 17% to 21% at the S&P 500 estimate. The 3S estimator attains the efficiency of generalized least squares with known variances, and its variance for the mean lies within 9% of that of the joint estimator in every design.

## Layout

```
data/     committed S&P 500 prices and returns
src/      analysis modules, every parameter declared once in config.R
report.pdf
notes.pdf
```

The report and the note are distributed as compiled PDFs. Their typesetting source is not included.

## Data

`data/SPX_daily.csv` holds the daily close of the S&P 500 index (^GSPC) from Yahoo Finance and its daily log return in percent, from December 1927 to September 2026.

## Reproducing

R 4.3.3 with `Rcpp` 1.0.12 and a C++ compiler (Rtools on Windows).

```
install.packages("Rcpp")
Rscript src/run.R
```

This prints every number quoted in the report and the note and writes the two figure files they use to `latex/`. The Monte Carlo samples, population values and S&P 500 estimates are cached in `cache/`. Both folders are created on the first run and are not tracked.
