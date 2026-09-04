// =============================================================================
//  L-Vasicek: distribuicao de Vasicek com NUCLEO LOGISTICO no intervalo (0,1)
//  -----------------------------------------------------------------------------
//  Duas parametrizacoes:
//    (A) LOCACAO-DISPERSAO  (alpha, theta): parametrizacao "original".
//        ATENCAO: alpha e' um parametro de LOCACAO, NAO a media. Sob o nucleo
//        logistico E[Y] != alpha em geral (a identidade E[Y]=alpha so vale no
//        nucleo NORMAL). Use a camada por quantil abaixo para regressao.
//    (B) QUANTIL           (mu, theta; quantile fixo): mu e' o quantile-esimo quantil
//        condicional exato. E' a parametrizacao recomendada para regressao
//        (analoga a NVASIQ do vasicekreg, porem com nucleo logistico).
//
//  Nucleo (parametrizacao de locacao):
//    logit(u)   = log(u/(1-u)),   Lambda(z) = 1/(1+exp(-z)),  lambda = Lambda(1-Lambda)
//    z2(x)      = ( sqrt(1-theta)*logit(x) - logit(alpha) ) / sqrt(theta)
//    CDF:  F(x) = Lambda( z2(x) )
//    PDF:  f(x) = sqrt((1-theta)/theta) * 1/(x(1-x)) * lambda( z2(x) )
//    QF:   Q(p) = Lambda( ( logit(alpha) + sqrt(theta)*logit(p) ) / sqrt(1-theta) )
//
//  Mapa de reparametrizacao (quantil -> locacao), fechado, sem raiz numerica:
//    logit(alpha) = sqrt(1-theta)*logit(mu) - sqrt(theta)*logit(quantile)
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Auxiliares numericamente estaveis
// -----------------------------------------------------------------------------

// logit(u) = log(u/(1-u))
inline double logit_(double u) {
    return std::log(u) - std::log1p(-u);
}

// softplus(z) = log(1 + exp(z)), estavel para |z| grande
inline double softplus_(double z) {
    if (z > 0.0)  return z + std::log1p(std::exp(-z));
    else          return std::log1p(std::exp(z));
}

// log Lambda(z)      = -softplus(-z) = -log(1+exp(-z))
inline double log_Lambda_(double z)   { return -softplus_(-z); }
// log(1 - Lambda(z)) = -softplus(z)  = -log(1+exp(z))
inline double log_1m_Lambda_(double z){ return -softplus_(z);  }
// log lambda(z)      = log Lambda + log(1-Lambda) = -(softplus(z)+softplus(-z))
inline double log_lambda_(double z)   { return -(softplus_(z) + softplus_(-z)); }
// Lambda(z)
inline double Lambda_(double z)       { return 1.0 / (1.0 + std::exp(-z)); }

// -----------------------------------------------------------------------------
// Nucleo — parametrizacao de LOCACAO (alpha) e DISPERSAO (theta)
// -----------------------------------------------------------------------------

// LOG-PDF
inline double logpdf_lvasicek(double x, double alpha, double theta) {
    if (!(x > 0.0 && x < 1.0))          return R_NegInf;
    if (!(alpha > 0.0 && alpha < 1.0))  return R_NegInf;
    if (!(theta > 0.0 && theta < 1.0))  return R_NegInf;

    const double lx = logit_(x);
    const double la = logit_(alpha);
    const double st  = std::sqrt(theta);
    const double s1t = std::sqrt(1.0 - theta);

    const double z = (s1t * lx - la) / st;                 // z2(x)

    const double log_const    = 0.5 * (std::log1p(-theta) - std::log(theta)); // 0.5*log((1-theta)/theta)
    const double log_jacobian = -std::log(x) - std::log1p(-x);                // -log[x(1-x)]
    const double log_kernel   = log_lambda_(z);                              // log lambda(z2)

    return log_const + log_jacobian + log_kernel;
}

// CDF (retorna tambem em log e cauda superior, em log-espaco quando pedido)
inline double logcdf_lvasicek(double x, double alpha, double theta, bool lower) {
    // limites
    if (x <= 0.0) return lower ? R_NegInf : 0.0;   // log(0) ou log(1)
    if (x >= 1.0) return lower ? 0.0      : R_NegInf;

    const double lx = logit_(x);
    const double la = logit_(alpha);
    const double st  = std::sqrt(theta);
    const double s1t = std::sqrt(1.0 - theta);
    const double z = (s1t * lx - la) / st;

    return lower ? log_Lambda_(z) : log_1m_Lambda_(z);
}

// QUANTIL
inline double invcdf_lvasicek(double p, double alpha, double theta,
                              bool lowertail = true,
                              bool logprob = false) {
    const double la = logit_(alpha);
    const double st  = std::sqrt(theta);
    const double s1t = std::sqrt(1.0 - theta);
    const double quantile_logit =
        R::qlogis(p, 0.0, 1.0, lowertail, logprob);

    const double z = (la + st * quantile_logit) / s1t;
    return Lambda_(z);
}

// -----------------------------------------------------------------------------
// Reparametrizacao quantil -> locacao:  alpha(mu, theta, quantile)
// -----------------------------------------------------------------------------
inline double alpha_from_mu(double mu, double theta, double quantile) {
    const double s1t = std::sqrt(1.0 - theta);
    const double st  = std::sqrt(theta);
    const double log_alpha = s1t * logit_(mu) - st * logit_(quantile); // = logit(alpha)
    return Lambda_(log_alpha);
}

// -----------------------------------------------------------------------------
// Guarda de reciclagem: comprimento deve ser 1 ou n
// -----------------------------------------------------------------------------
inline void check_len(int len, int n, const char* nm) {
    if (len != 1 && len != n)
        Rcpp::stop("Comprimento de '%s' deve ser 1 ou igual a length(x).", nm);
}

// =============================================================================
//  EXPORTS — Parametrizacao de LOCACAO (alpha, theta)
// =============================================================================

// [[Rcpp::export]]
NumericVector cpp_dlvasicek(const NumericVector x,
                            const NumericVector alpha,
                            const NumericVector theta,
                            const bool logprob = false) {
    const int n = x.length();
    const int na = alpha.length(), nt = theta.length();
    check_len(na, n, "alpha"); check_len(nt, n, "theta");

    NumericVector out(n);
    for (int i = 0; i < n; ++i)
        out[i] = logpdf_lvasicek(x[i], alpha[na == 1 ? 0 : i], theta[nt == 1 ? 0 : i]);

    if (logprob) return out;
    return NumericVector(Rcpp::exp(out));
}

// [[Rcpp::export]]
NumericVector cpp_plvasicek(const NumericVector x,
                            const NumericVector alpha,
                            const NumericVector theta,
                            const bool lowertail = true,
                            const bool logprob = false) {
    const int n = x.length();
    const int na = alpha.length(), nt = theta.length();
    check_len(na, n, "alpha"); check_len(nt, n, "theta");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double lp = logcdf_lvasicek(x[i], alpha[na == 1 ? 0 : i],
                                    theta[nt == 1 ? 0 : i], lowertail);
        out[i] = logprob ? lp : std::exp(lp);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_qlvasicek(const NumericVector p,
                            const NumericVector alpha,
                            const NumericVector theta,
                            const bool lowertail = true,
                            const bool logprob = false) {
    const int n = p.length();
    const int na = alpha.length(), nt = theta.length();
    check_len(na, n, "alpha"); check_len(nt, n, "theta");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        out[i] = invcdf_lvasicek(
            p[i],
            alpha[na == 1 ? 0 : i],
            theta[nt == 1 ? 0 : i],
            lowertail,
            logprob
        );
    }
    return out;
}

// Gerador (vetorizado em alpha, theta) — necessario para Monte Carlo com covariaveis
// [[Rcpp::export]]
NumericVector cpp_rlvasicek(const int n,
                            const NumericVector alpha,
                            const NumericVector theta) {
    const int na = alpha.length(), nt = theta.length();
    check_len(na, n, "alpha"); check_len(nt, n, "theta");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double u = R::runif(0.0, 1.0);
        out[i] = invcdf_lvasicek(
            u,
            alpha[na == 1 ? 0 : i],
            theta[nt == 1 ? 0 : i]
        );
    }
    return out;
}

// =============================================================================
//  EXPORTS — Parametrizacao por QUANTIL (mu, theta; quantile fixo)
//  mu = quantil condicional exato de ordem quantile. Recomendada para regressao.
// =============================================================================

// [[Rcpp::export]]
NumericVector cpp_dLVASIQ(const NumericVector x,
                          const NumericVector mu,
                          const NumericVector theta,
                          const double quantile = 0.5,
                          const bool logprob = false) {
    const int n = x.length();
    const int nm = mu.length(), nt = theta.length();
    check_len(nm, n, "mu"); check_len(nt, n, "theta");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double th = theta[nt == 1 ? 0 : i];
        double a  = alpha_from_mu(mu[nm == 1 ? 0 : i], th, quantile);
        out[i] = logpdf_lvasicek(x[i], a, th);
    }
    if (logprob) return out;
    return NumericVector(Rcpp::exp(out));
}

// [[Rcpp::export]]
NumericVector cpp_pLVASIQ(const NumericVector x,
                          const NumericVector mu,
                          const NumericVector theta,
                          const double quantile = 0.5,
                          const bool lowertail = true,
                          const bool logprob = false) {
    const int n = x.length();
    const int nm = mu.length(), nt = theta.length();
    check_len(nm, n, "mu"); check_len(nt, n, "theta");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double th = theta[nt == 1 ? 0 : i];
        double a  = alpha_from_mu(mu[nm == 1 ? 0 : i], th, quantile);
        double lp = logcdf_lvasicek(x[i], a, th, lowertail);
        out[i] = logprob ? lp : std::exp(lp);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_qLVASIQ(const NumericVector p,
                          const NumericVector mu,
                          const NumericVector theta,
                          const double quantile = 0.5,
                          const bool lowertail = true,
                          const bool logprob = false) {
    const int n = p.length();
    const int nm = mu.length(), nt = theta.length();
    check_len(nm, n, "mu"); check_len(nt, n, "theta");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double th = theta[nt == 1 ? 0 : i];
        double a  = alpha_from_mu(mu[nm == 1 ? 0 : i], th, quantile);
        out[i] = invcdf_lvasicek(
            p[i], a, th, lowertail, logprob
        );
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_rLVASIQ(const int n,
                          const NumericVector mu,
                          const NumericVector theta,
                          const double quantile = 0.5) {
    const int nm = mu.length(), nt = theta.length();
    check_len(nm, n, "mu"); check_len(nt, n, "theta");
    if (!(quantile > 0.0 && quantile < 1.0)) Rcpp::stop("'quantile' must be in (0, 1).");

    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double th = theta[nt == 1 ? 0 : i];
        double a  = alpha_from_mu(mu[nm == 1 ? 0 : i], th, quantile);
        double u  = R::runif(0.0, 1.0);
        out[i] = invcdf_lvasicek(u, a, th);
    }
    return out;
}

// =============================================================================
//  Utilitario opcional: media numerica sob parametrizacao de LOCACAO.
//  E[Y] nao tem forma fechada elementar no nucleo logistico; aqui via
//  quadratura de Gauss-Legendre no espaco do quantil (rapida e estavel).
//  Serve para recuperar E[Y] quando se usa a parametrizacao (alpha, theta).
// =============================================================================
// [[Rcpp::export]]
NumericVector cpp_mean_lvasicek(const NumericVector alpha,
                                const NumericVector theta,
                                const int npts = 128) {
    const int na = alpha.length(), nt = theta.length();
    const int n = std::max(na, nt);
    // Regra composta do ponto medio em (0,1), suficientemente precisa para
    // E[Y] = int_0^1 Q(p) dp quando npts e grande.
    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double a  = alpha[na == 1 ? 0 : i % na];
        double th = theta[nt == 1 ? 0 : i % nt];
        double acc = 0.0;
        for (int k = 0; k < npts; ++k) {
            double p = (k + 0.5) / npts;          // ponto medio
            acc += invcdf_lvasicek(p, a, th);
        }
        out[i] = acc / npts;                       // E[Y] = int_0^1 Q(p) dp
    }
    return out;
}
