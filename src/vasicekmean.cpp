#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// log-pdf Vasicek - mean parameterization

inline double logpdf_vasicekmean(double x, double alpha, double theta)
{
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double qnorma = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    double t1 = 0.1e1 - theta;
    double t2 = log(t1);
    double t4 = log(theta);
    double t6 = qnormx * qnormx;
    double t8 = sqrt(t1);
    double t11 = pow(qnormx * t8 - qnorma, 0.2e1);
    return(0.05e1 * t2 - 0.05e1 * t4 + 0.05e1 * t6 - 0.05e1 * t11 / theta);
}

// [[Rcpp::export]]
NumericVector cpp_dvasicekmean(const NumericVector x,
                               const NumericVector alpha,
                               const NumericVector theta,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nalpha = alpha.length(); 
    const int ntheta = theta.length();
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = logpdf_vasicekmean(x[i], alpha[i % nalpha], theta[i % ntheta]);
    
    if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf Vasicek - mean parameterization

inline double cdf_vasicekmean(double x, double alpha, double theta)
{
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double qnorma = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    double t1 = (qnormx * pow(0.1e1 - theta, 0.05e1) - qnorma) / pow(theta, 0.05e1);
    return(R::pnorm(t1, 0.0, 1.0, TRUE, FALSE));
}

// [[Rcpp::export]]
NumericVector cpp_pvasicekmean(const NumericVector x,
                               const NumericVector alpha,
                               const NumericVector theta,
                               const bool lowertail = true,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nalpha = alpha.length(); 
    const int ntheta = theta.length();
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = cdf_vasicekmean(x[i], alpha[i % nalpha], theta[i % ntheta]);
    
    if (!lowertail) out = 0.1e1 - out;
    if (logprob) out = Rcpp::log(out);
    
    return(out);
}

 // inv-cdf vasicekmean

 inline double invcdf_vasicekmean(double x, double alpha, double theta)
 {
     double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
     double qnorma = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
     double t1 = (qnorma + qnormx * pow(theta, 0.05e1)) / pow(0.1e1 - theta, 0.05e1);
     return(R::pnorm(t1, 0.0, 1.0, TRUE, FALSE));
 }
 
 // [[Rcpp::export]]
 NumericVector cpp_qvasicekmean(const NumericVector x,
                                const NumericVector alpha,
                                const NumericVector theta,
                                const bool lowertail = true,
                                const bool logprob = false)
 {
     const int n = x.length(); 
     const int nalpha = alpha.length(); 
     const int ntheta = theta.length();
     NumericVector out(n);
     
     if(lowertail)
     {
     for(int i = 0; i < n; i++)
        out[i] = invcdf_vasicekmean(x[i], alpha[i % nalpha], theta[i % ntheta]);
     }
     else
     {
     for(int i = 0; i < n; i++)
        out[i] = invcdf_vasicekmean(0.1e1 - x[i], alpha[i % nalpha], theta[i % ntheta]);
     }
     
     if(logprob) return(Rcpp::log(out)); else return(out);
 }
