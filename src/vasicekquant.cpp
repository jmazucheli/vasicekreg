#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// log-pdf Vasicek - quant parameterization

inline double logpdf_vasicekquant(double x, double mu, double theta, double tau)
{
    double qnormmu = R::qnorm(mu, 0.0, 1.0, TRUE, FALSE);
    double qnormtau = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE);
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double alpha = R::pnorm(-sqrt(theta) * qnormtau + qnormmu * sqrt(0.1e1 - theta), 0.0, 1.0, TRUE, FALSE);
    double qnormalpha = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    double t2 = 0.1e1 - theta;
    double t3 = log(t2);
    double t5 = log(theta);
    double t7 = qnormx * qnormx;
    double t9 = sqrt(t2);
    double t13 = pow(t9 * qnormx - qnormalpha, 0.2e1);
    return(0.5e0 * (t3 - t5) + 0.5e0 * t7 - 0.5e0 * t13 / theta);
    }

// [[Rcpp::export]]
NumericVector cpp_dvasicekquant(const NumericVector x,
                               const NumericVector mu,
                               const NumericVector theta,
                               const NumericVector tau,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nmu = mu.length(); 
    const int ntheta = theta.length();
    const int ntau = tau.length();
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = logpdf_vasicekquant(x[i], mu[i % nmu], theta[i % ntheta], tau[i % ntau]);
    
    if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf Vasicek - quant parameterization

inline double cdf_vasicekquant(double x, double mu, double theta, double tau)
{
    double qnormmu = R::qnorm(mu, 0.0, 1.0, TRUE, FALSE);
    double qnormtau = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE);
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double alpha = R::pnorm(-sqrt(theta) * qnormtau + qnormmu * sqrt(0.1e1 - theta), 0.0, 1.0, TRUE, FALSE);
    double qnormalpha = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    return(R::pnorm((sqrt(0.1e1 - theta) * qnormx - qnormalpha) / sqrt(theta), 0.0, 1.0, TRUE, FALSE));
}

// [[Rcpp::export]]
NumericVector cpp_pvasicekquant(const NumericVector x,
                               const NumericVector mu,
                               const NumericVector theta,
                               const NumericVector tau,
                               const bool lowertail = true,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nmu = mu.length(); 
    const int ntheta = theta.length();
    const int ntau = tau.length();
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = cdf_vasicekquant(x[i], mu[i % nmu], theta[i % ntheta], tau[i % ntau]);
    
    if (!lowertail) out = 0.1e1 - out;
    if (logprob) out = Rcpp::log(out);
    
    return(out);
}

// inv-cdf vasicekquant

inline double invcdf_vasicekquant(double x, double mu, double theta, double tau)
{
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double qnormmu = R::qnorm(mu, 0.0, 1.0, TRUE, FALSE);
    double qnormtau = R::qnorm(tau, 0.0, 1.0, TRUE, FALSE);
    double alpha = R::pnorm(-sqrt(theta) * qnormtau + qnormmu * sqrt(0.1e1 - theta), 0.0, 1.0, TRUE, FALSE);
    double qnormalpha = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE); ;
    return(R::pnorm((qnormalpha + sqrt(theta) * qnormx) / sqrt(0.1e1 - theta),0.0, 1.0, TRUE, FALSE)); 
}

// [[Rcpp::export]]
NumericVector cpp_qvasicekquant(const NumericVector x,
                                const NumericVector mu,
                                const NumericVector theta,
                                const NumericVector tau,
                                const bool lowertail = true,
                                const bool logprob = false)
{
    const int n = x.length(); 
    const int nmu = mu.length(); 
    const int ntheta = theta.length();
    const int ntau = tau.length();
    NumericVector out(n);
    
    if(lowertail)
    {
        for(int i = 0; i < n; i++)
            out[i] = invcdf_vasicekquant(x[i], mu[i % nmu], theta[i % ntheta], tau[i % ntau] );
    }
    else
    {
        for(int i = 0; i < n; i++)
            out[i] = invcdf_vasicekquant(0.1e1 - x[i], mu[i % nmu], theta[i % ntheta], tau[i % ntau]);
    }
    
    if(logprob) return(Rcpp::log(out)); else return(out);
}
