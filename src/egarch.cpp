// AR(p)-EGARCH(1,1): simulation, Gaussian quasi-likelihood with analytic gradient, filter,
// and population moments from the derivative recursion at the true parameter.
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static inline double sgn(double x) { return (x > 0) - (x < 0); }

// Path y of model (2.1) driven by eta, without the first burn values.
// [[Rcpp::export]]
NumericVector sim_ar_egarch(NumericVector eta, int burn, double nu, NumericVector phi,
                            double omega, double alpha, double gamma, double beta, double c) {
  int N = eta.size(), p = phi.size();
  NumericVector y(N), h(N);
  h[0] = omega / (1.0 - beta);
  for (int t = 0; t < N; ++t) {
    if (t > 0)
      h[t] = omega + alpha * (std::fabs(eta[t-1]) - c) + gamma * eta[t-1] + beta * h[t-1];
    double m = nu;
    for (int i = 1; i <= p; ++i)
      if (t - i >= 0) m += phi[i-1] * (y[t-i] - nu);
    y[t] = m + std::exp(h[t] / 2.0) * eta[t];
  }
  return y[Range(burn, N - 1)];
}

// Average quasi-log-likelihood (2.3) and its almost-everywhere gradient, filter started at
// h0 for t = p, theta = (nu, phi, omega, alpha, gamma, beta). nu_known fixes nu = theta[0].
// [[Rcpp::export]]
List qll(NumericVector theta, NumericVector y, int p, double c, double h0, bool nu_known) {
  int n = y.size(), k = p + 5;
  double nu = theta[0];
  double om = theta[p+1], al = theta[p+2], ga = theta[p+3], be = theta[p+4];
  double sphi = 0;
  for (int i = 0; i < p; ++i) sphi += theta[1+i];
  std::vector<double> dh(k, 0.0), deps(k, 0.0), g(k, 0.0), dhn(k);
  double h = h0, L = 0;
  for (int t = p; t < n; ++t) {
    double e = y[t] - nu;
    for (int i = 1; i <= p; ++i) e -= theta[i] * (y[t-i] - nu);
    std::fill(deps.begin(), deps.end(), 0.0);
    deps[0] = -(1.0 - sphi);
    for (int i = 1; i <= p; ++i) deps[i] = -(y[t-i] - nu);
    double eh = std::exp(-h), z = e * std::exp(-h / 2.0);
    L += -0.5 * (h + e * e * eh);
    for (int j = 0; j < k; ++j)
      g[j] += -0.5 * (1.0 - z * z) * dh[j] - z * std::exp(-h / 2.0) * deps[j];
    double A = be - 0.5 * (al * std::fabs(z) + ga * z);
    double s = (al * sgn(z) + ga) * std::exp(-h / 2.0);
    for (int j = 0; j < k; ++j) dhn[j] = A * dh[j] + s * deps[j];
    dhn[p+1] += 1.0;
    dhn[p+2] += std::fabs(z) - c;
    dhn[p+3] += z;
    dhn[p+4] += h;
    h = om + al * (std::fabs(z) - c) + ga * z + be * h;
    dh = dhn;
    if (!std::isfinite(h) || std::fabs(h) > 700)
      return List::create(_["value"] = R_NegInf, _["grad"] = NumericVector(k));
  }
  double T = n - p;
  NumericVector G(k);
  for (int j = 0; j < k; ++j) G[j] = g[j] / T;
  if (nu_known) G[0] = 0.0;
  return List::create(_["value"] = L / T, _["grad"] = G);
}

// Feasible filter (2.2) at theta, started at h0, for t = p, ..., n - 1.
// [[Rcpp::export]]
NumericVector filter_h(NumericVector theta, NumericVector y, int p, double c, double h0) {
  int n = y.size();
  double nu = theta[0], om = theta[p+1], al = theta[p+2], ga = theta[p+3], be = theta[p+4];
  NumericVector H(n - p);
  double h = h0;
  for (int t = p; t < n; ++t) {
    double e = y[t] - nu;
    for (int i = 1; i <= p; ++i) e -= theta[i] * (y[t-i] - nu);
    H[t-p] = h;
    double z = e * std::exp(-h / 2.0);
    h = om + al * (std::fabs(z) - c) + ga * z + be * h;
  }
  return H;
}

// Moments at theta0 over a long path driven by eta: C = E[a a'], I = E[s s'],
// G = E[s xi' eps], W = E[sigma^-2 xi xi'], M = E[xi xi'], S = E[sigma^2 xi xi'],
// with xi = d eps / d m, a the derivative series and s the score (3.1).
// [[Rcpp::export]]
List pop_moments(NumericVector eta, int burn, double nu, NumericVector phi,
                 double omega, double alpha, double gamma, double beta, double c) {
  int N = eta.size(), p = phi.size(), k = p + 5, m = p + 1;
  double sphi = 0;
  for (int i = 0; i < p; ++i) sphi += phi[i];
  std::vector<double> y(N), h(N), a(k, 0.0), an(k), xi(m), s(k);
  h[0] = omega / (1.0 - beta);
  NumericMatrix C(k, k), I(k, k), G(k, m), W(m, m), M(m, m), S(m, m);
  double cnt = 0;
  for (int t = 0; t < N; ++t) {
    if (t > 0)
      h[t] = omega + alpha * (std::fabs(eta[t-1]) - c) + gamma * eta[t-1] + beta * h[t-1];
    double mu = nu;
    for (int i = 1; i <= p; ++i)
      if (t - i >= 0) mu += phi[i-1] * (y[t-i] - nu);
    double sig = std::exp(h[t] / 2.0), et = eta[t];
    y[t] = mu + sig * et;
    if (t < p) continue;
    xi[0] = -(1.0 - sphi);
    for (int i = 1; i <= p; ++i) xi[i] = -(y[t-i] - nu);
    if (t >= burn) {
      for (int j = 0; j < k; ++j)
        s[j] = -0.5 * (1 - et * et) * a[j] - (j < m ? et * xi[j] / sig : 0.0);
      for (int r = 0; r < k; ++r) {
        for (int q = 0; q < k; ++q) {
          C(r, q) += a[r] * a[q];
          I(r, q) += s[r] * s[q];
        }
        for (int q = 0; q < m; ++q) G(r, q) += s[r] * xi[q] * sig * et;
      }
      for (int r = 0; r < m; ++r)
        for (int q = 0; q < m; ++q) {
          double x = xi[r] * xi[q];
          W(r, q) += x / (sig * sig);
          M(r, q) += x;
          S(r, q) += x * sig * sig;
        }
      cnt += 1;
    }
    // a_{t+1} = w_t + A_t a_t
    double A = beta - 0.5 * (alpha * std::fabs(et) + gamma * et);
    double sc = (alpha * sgn(et) + gamma) / sig;
    for (int j = 0; j < k; ++j) an[j] = A * a[j];
    for (int j = 0; j < m; ++j) an[j] += sc * xi[j];
    an[p+1] += 1.0;
    an[p+2] += std::fabs(et) - c;
    an[p+3] += et;
    an[p+4] += h[t];
    a = an;
  }
  for (int r = 0; r < k; ++r) {
    for (int q = 0; q < k; ++q) {
      C(r, q) /= cnt;
      I(r, q) /= cnt;
    }
    for (int q = 0; q < m; ++q) G(r, q) /= cnt;
  }
  for (int r = 0; r < m; ++r)
    for (int q = 0; q < m; ++q) {
      W(r, q) /= cnt;
      M(r, q) /= cnt;
      S(r, q) /= cnt;
    }
  return List::create(_["C"] = C, _["I"] = I, _["G"] = G, _["W"] = W, _["M"] = M, _["S"] = S);
}
