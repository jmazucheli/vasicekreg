#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline void check_len_mean(int len, int n, const char* name)
{
    if (len != 1 && len != n)
        Rcpp::stop("Length of '%s' must be 1 or equal to the data length.",
                   name);
}

// log-pdf Vasicek - mean parameterization

inline double logpdf_NvasicekM(double x, double alpha, double theta)
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
NumericVector cpp_dNVASIM(const NumericVector x,
                               const NumericVector alpha,
                               const NumericVector theta,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nalpha = alpha.length(); 
    const int ntheta = theta.length();
    check_len_mean(nalpha, n, "mu");
    check_len_mean(ntheta, n, "sigma");
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = logpdf_NvasicekM(
            x[i],
            alpha[nalpha == 1 ? 0 : i],
            theta[ntheta == 1 ? 0 : i]
        );
    
    if(logprob) return(out); else return(Rcpp::exp(out));
}

// cdf Vasicek - mean parameterization

inline double cdf_NvasicekM(double x, double alpha, double theta,
                              bool lowertail, bool logprob)
{
    double qnormx = R::qnorm(x, 0.0, 1.0, TRUE, FALSE);
    double qnorma = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
    double t1 = (qnormx * pow(0.1e1 - theta, 0.05e1) - qnorma) / pow(theta, 0.05e1);
    return(R::pnorm(t1, 0.0, 1.0, lowertail, logprob));
}

// [[Rcpp::export]]
NumericVector cpp_pNVASIM(const NumericVector x,
                               const NumericVector alpha,
                               const NumericVector theta,
                               const bool lowertail = true,
                               const bool logprob = false)
{
    const int n = x.length(); 
    const int nalpha = alpha.length(); 
    const int ntheta = theta.length();
    check_len_mean(nalpha, n, "mu");
    check_len_mean(ntheta, n, "sigma");
    NumericVector out(n);
    
    for(int i = 0; i < n; i++)
        out[i] = cdf_NvasicekM(
            x[i],
            alpha[nalpha == 1 ? 0 : i],
            theta[ntheta == 1 ? 0 : i],
            lowertail,
            logprob
        );
    
    return(out);
}

 // inv-cdf NvasicekM

 inline double invcdf_NvasicekM(double p, double alpha, double theta,
                                  bool lowertail, bool logprob)
 {
     double qnormx = R::qnorm(p, 0.0, 1.0, lowertail, logprob);
     double qnorma = R::qnorm(alpha, 0.0, 1.0, TRUE, FALSE);
     double t1 = (qnorma + qnormx * pow(theta, 0.05e1)) / pow(0.1e1 - theta, 0.05e1);
     return(R::pnorm(t1, 0.0, 1.0, TRUE, FALSE));
 }
 
 // [[Rcpp::export]]
 NumericVector cpp_qNVASIM(const NumericVector x,
                                const NumericVector alpha,
                                const NumericVector theta,
                                const bool lowertail = true,
                                const bool logprob = false)
 {
     const int n = x.length(); 
     const int nalpha = alpha.length(); 
     const int ntheta = theta.length();
     check_len_mean(nalpha, n, "mu");
     check_len_mean(ntheta, n, "sigma");
     NumericVector out(n);
     
     for(int i = 0; i < n; i++)
         out[i] = invcdf_NvasicekM(
             x[i],
             alpha[nalpha == 1 ? 0 : i],
             theta[ntheta == 1 ? 0 : i],
             lowertail,
             logprob
         );
     
     return(out);
 }
