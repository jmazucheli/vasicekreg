#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline void check_len_quant(int len, int n, const char* name)
{
    if (len != 1 && len != n)
        Rcpp::stop("Length of '%s' must be 1 or equal to the data length.",
                   name);
}

// log-pdf Vasicek - quant parameterization

inline double logpdf_NvasicekQ(double x, double mu, double theta, double quantile)
{
    double qnormmu = R::qnorm(mu, 0.0, 1.0, TRUE, FALSE);
    double qnorm_quantile = R::qnorm(quantile, 0.0, 1.0, TRUE, FALSE);
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double qnormalpha = -sqrt(theta) * qnorm_quantile +
        qnormmu * sqrt(0.1e1 - theta);
    double t2 = 0.1e1 - theta;
    double t3 = log(t2);
    double t5 = log(theta);
    double t7 = qnormx * qnormx;
    double t9 = sqrt(t2);
    double t13 = pow(t9 * qnormx - qnormalpha, 0.2e1);
    return(0.5e0 * (t3 - t5) + 0.5e0 * t7 - 0.5e0 * t13 / theta);
    }

// [[Rcpp::export]]
NumericVector cpp_dNVASIQ(const NumericVector x,
                               const NumericVector mu,
                               const NumericVector theta,
                               const double quantile,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nmu = mu.length(); 
    const int ntheta = theta.length();
    check_len_quant(nmu, n, "mu");
    check_len_quant(ntheta, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0))
        Rcpp::stop("'quantile' must be in (0, 1).");
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = logpdf_NvasicekQ(
            x[i],
            mu[nmu == 1 ? 0 : i],
            theta[ntheta == 1 ? 0 : i],
            quantile
        );
    
    if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf Vasicek - quant parameterization

inline double cdf_NvasicekQ(double x, double mu, double theta, double quantile,
                               bool lowertail, bool logprob)
{
    double qnormmu = R::qnorm(mu, 0.0, 1.0, TRUE, FALSE);
    double qnorm_quantile = R::qnorm(quantile, 0.0, 1.0, TRUE, FALSE);
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double qnormalpha = -sqrt(theta) * qnorm_quantile +
        qnormmu * sqrt(0.1e1 - theta);
    return(R::pnorm(
        (sqrt(0.1e1 - theta) * qnormx - qnormalpha) / sqrt(theta),
        0.0, 1.0, lowertail, logprob
    ));
}

// [[Rcpp::export]]
NumericVector cpp_pNVASIQ(const NumericVector x,
                               const NumericVector mu,
                               const NumericVector theta,
                               const double quantile,
                               const bool lowertail = true,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nmu = mu.length(); 
    const int ntheta = theta.length();
    check_len_quant(nmu, n, "mu");
    check_len_quant(ntheta, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0))
        Rcpp::stop("'quantile' must be in (0, 1).");
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = cdf_NvasicekQ(
            x[i],
            mu[nmu == 1 ? 0 : i],
            theta[ntheta == 1 ? 0 : i],
            quantile,
            lowertail,
            logprob
        );
    
    return(out);
}

// inv-cdf NvasicekQ

inline double invcdf_NvasicekQ(double p, double mu, double theta,
                                  double quantile, bool lowertail, bool logprob)
{
    double qnormx = R::qnorm(p, 0.0, 1.0, lowertail, logprob);
    double qnormmu = R::qnorm(mu, 0.0, 1.0, TRUE, FALSE);
    double qnorm_quantile = R::qnorm(quantile, 0.0, 1.0, TRUE, FALSE);
    double qnormalpha = -sqrt(theta) * qnorm_quantile +
        qnormmu * sqrt(0.1e1 - theta);
    return(R::pnorm((qnormalpha + sqrt(theta) * qnormx) / sqrt(0.1e1 - theta),0.0, 1.0, TRUE, FALSE)); 
}

// [[Rcpp::export]]
NumericVector cpp_qNVASIQ(const NumericVector x,
                                const NumericVector mu,
                                const NumericVector theta,
                                const double quantile,
                                const bool lowertail = true,
                                const bool logprob = false)
{
    const int n = x.length(); 
    const int nmu = mu.length(); 
    const int ntheta = theta.length();
    check_len_quant(nmu, n, "mu");
    check_len_quant(ntheta, n, "sigma");
    if (!(quantile > 0.0 && quantile < 1.0))
        Rcpp::stop("'quantile' must be in (0, 1).");
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = invcdf_NvasicekQ(
            x[i],
            mu[nmu == 1 ? 0 : i],
            theta[ntheta == 1 ? 0 : i],
            quantile,
            lowertail,
            logprob
        );
    
    return(out);
}
