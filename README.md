# Efficiency Loss of Sequential Quasi-Maximum Likelihood Estimation in the AR(p)-EGARCH(1,1) Model

This study examines which parameters of the AR(p)-EGARCH(1,1) model are estimated less efficiently by sequential than by joint Gaussian quasi-maximum likelihood estimation. The two-step (2S) estimator fits the mean by ordinary least squares and then the variance equation by quasi-maximum likelihood, and the three-step (3S) estimator refits the mean by feasible generalized least squares.

## Findings

The information matrix is not block diagonal between the mean and the variance parameters, even without a leverage effect, through a single entry that links the mean return to the coefficient on the lagged innovation in the log variance. Without leverage, the 2S estimator therefore loses efficiency only for this coefficient, and only when the mean return is estimated.

With leverage, the loss extends to the intercept and the persistence of the log variance at second order in the leverage coefficient. It is largest for the mean log variance, about 15% in the base Monte Carlo design and 17% to 21% at the S&P 500 estimate. The 3S estimator recovers most of the efficiency for the mean parameters.

## Layout

```
data/     committed S&P 500 prices and returns
src/      analysis modules, every parameter declared once in config.R
report.pdf
notes.pdf an introduction to AR-EGARCH models and quasi-maximum likelihood
```

The report and the note are distributed as compiled PDFs. Their typesetting source is not included.

## Data

`data/SPX_daily.csv` holds the daily close of the S&P 500 index (^GSPC) from Yahoo Finance and its daily log return in percent, from December 1927 to September 2026. The report uses January 1970 to August 2026.

## Reproducing

R 4.3.3 with `Rcpp` 1.0.12 and a C++ compiler (Rtools on Windows).

```
install.packages("Rcpp")
Rscript src/run.R
```

This prints every number quoted in the report and the note and draws their two figures, caching the Monte Carlo samples, population values and S&P 500 estimates on the first run.
