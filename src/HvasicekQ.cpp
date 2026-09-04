#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace {

const double PI_ = 3.141592653589793238462643383279502884;
const double LOG_PI_ = 1.144729885849400174143427351353058711;
const double LOG_TWO_ = 0.693147180559945309417232121458176568;

inline double qhs_(const double u) {
    return std::log(std::tan(0.5 * PI_ * u));
}

inline double H_(const double z) {
    if (z > 0.0) {
        return 1.0 - (2.0 / PI_) * std::atan(std::exp(-z));
    }
    return (2.0 / PI_) * std::atan(std::exp(z));
}

inline double logcosh_(const double z) {
    const double az = std::fabs(z);
    return az + std::log1p(std::exp(-2.0 * az)) - LOG_TWO_;
}

inline double logh_(const double z) {
    return -LOG_PI_ - logcosh_(z);
}

inline double logH_(const double z, const bool lower_tail) {
    double tail;
    if (z >= 0.0) {
        tail = (2.0 / PI_) * std::atan(std::exp(-z));
        return lower_tail ? std::log1p(-tail) : std::log(tail);
    }
    tail = (2.0 / PI_) * std::atan(std::exp(z));
    return lower_tail ? std::log(tail) : std::log1p(-tail);
}

inline void check_len_(const int len, const int n, const char* name) {
    if (len != 1 && len != n) {
        Rcpp::stop("Length of '%s' must be 1 or equal to the input length.", name);
    }
}

inline double logpdf_hsvasique(const double x, const double mu,
                               const double sigma, const double quantile) {
    if (!(x > 0.0 && x < 1.0)) return R_NegInf;
    if (!(mu > 0.0 && mu < 1.0)) return NA_REAL;
    if (!(sigma > 0.0 && sigma < 1.0)) return NA_REAL;

    const double a = std::sqrt((1.0 - sigma) / sigma);
    const double qx = qhs_(x);
    const double z = a * (qx - qhs_(mu)) + qhs_(quantile);
    return std::log(a) + logh_(z) - logh_(qx);
}

inline double logcdf_hsvasique(const double x, const double mu,
                               const double sigma, const double quantile,
                               const bool lower_tail) {
    if (x <= 0.0) return lower_tail ? R_NegInf : 0.0;
    if (x >= 1.0) return lower_tail ? 0.0 : R_NegInf;

    const double a = std::sqrt((1.0 - sigma) / sigma);
    const double z = a * (qhs_(x) - qhs_(mu)) + qhs_(quantile);
    return logH_(z, lower_tail);
}

inline double quantile_hsvasique(const double p, const double mu,
                                 const double sigma, const double quantile,
                                 const bool lower_tail,
                                 const bool log_p) {
    double probability;
    if (log_p) {
        if (lower_tail) {
            if (p == R_NegInf) return 0.0;
            if (p == 0.0) return 1.0;
            probability = std::exp(p);
        } else {
            if (p == R_NegInf) return 1.0;
            if (p == 0.0) return 0.0;
            probability = -std::expm1(p);
        }
    } else {
        probability = lower_tail ? p : 1.0 - p;
    }

    if (probability <= 0.0) return 0.0;
    if (probability >= 1.0) return 1.0;

    const double scale = std::sqrt(sigma / (1.0 - sigma));
    const double value = qhs_(mu) + scale * (qhs_(probability) - qhs_(quantile));
    return H_(value);
}

} // namespace

// [[Rcpp::export]]
NumericVector cpp_dHVASIQ(const NumericVector x,
                           const NumericVector mu,
                           const NumericVector sigma,
                           const double quantile = 0.5,
                           const bool logprob = false) {
    const int n = x.length();
    const int nm = mu.length();
    const int ns = sigma.length();
    check_len_(nm, n, "mu");
    check_len_(ns, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        const double value = logpdf_hsvasique(
            x[i], mu[nm == 1 ? 0 : i], sigma[ns == 1 ? 0 : i], quantile
        );
        out[i] = logprob ? value : std::exp(value);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_pHVASIQ(const NumericVector q,
                           const NumericVector mu,
                           const NumericVector sigma,
                           const double quantile = 0.5,
                           const bool lower_tail = true,
                           const bool log_p = false) {
    const int n = q.length();
    const int nm = mu.length();
    const int ns = sigma.length();
    check_len_(nm, n, "mu");
    check_len_(ns, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        const double value = logcdf_hsvasique(
            q[i], mu[nm == 1 ? 0 : i], sigma[ns == 1 ? 0 : i], quantile,
            lower_tail
        );
        out[i] = log_p ? value : std::exp(value);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_qHVASIQ(const NumericVector p,
                           const NumericVector mu,
                           const NumericVector sigma,
                           const double quantile = 0.5,
                           const bool lower_tail = true,
                           const bool log_p = false) {
    const int n = p.length();
    const int nm = mu.length();
    const int ns = sigma.length();
    check_len_(nm, n, "mu");
    check_len_(ns, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        out[i] = quantile_hsvasique(
            p[i], mu[nm == 1 ? 0 : i], sigma[ns == 1 ? 0 : i], quantile,
            lower_tail, log_p
        );
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_rHVASIQ(const int n,
                           const NumericVector mu,
                           const NumericVector sigma,
                           const double quantile = 0.5) {
    const int nm = mu.length();
    const int ns = sigma.length();
    check_len_(nm, n, "mu");
    check_len_(ns, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        out[i] = quantile_hsvasique(
            R::runif(0.0, 1.0), mu[nm == 1 ? 0 : i],
            sigma[ns == 1 ? 0 : i], quantile, true, false
        );
    }
    return out;
}
