# test_that("Analitycal derivatives for alpha theta works", {
# library(numDeriv)
#     
# lnfalpha <- function(P){
#     alpha <- P[1]; 
#     theta <- P[2];
#     qnormx <- qnorm(x);
#     t1 <- 0.1e1 - theta;
#     t2 <- log(t1);
#     t4 <- log(theta);
#     t6 <- qnormx * qnormx;
#     t8 <- sqrt(t1);
#     t10 <- qnorm(alpha);
#     t12 <- pow(qnormx * t8 - t10, 0.2e1);
#     return(0.5e0 * t2 - 0.5e0 * t4 + 0.5e0 * t6 - 0.5e0 * t12 / theta);
# }
#     
# x       <- 0.37;
# theta   <- 0.37; 
# alpha   <- 0.69;
# grad(lnfalpha, c(alpha, theta))
# d1alpha(x, theta, alpha)
# d1theta(x, theta, alpha)
# 
# hessian(lnfalpha, c(alpha, theta))
# d2alpha(x, theta, alpha)
# d2theta(x, theta, alpha)
# dalphadtheta(x, theta, alpha)})

#########################################################
# test_that("Analitycal derivatives for mu theta works", {
#     library(numDeriv)
#     
# lnfmu <- function(P){
#     theta <- P[1]; 
#     mu <- P[2];
#     qnormx <- qnorm(x);
#     t1 <- 0.1e1 - theta;
#     t2 <- log(t1);
#     t4 <- log(theta);
#     t6 <- qnormx * qnormx;
#     t8 <- sqrt(t1);
#     t10 <- qnorm(mu);
#     t12 <- qnorm(tau);
#     t13 <- sqrt(theta);
#     t16 <- pow(qnormx * t8 - t10 * t8 + t12 * t13, 0.2e1);
#     return(t2 / 0.2e1 - t4 / 0.2e1 + 0.5e0 * t6 - 0.5e0 * t16 / theta);
# }
# 
# x       <- 0.37;
# theta   <- 0.37; 
# mu      <- 0.69;
# tau     <- 0.96;
# grad(lnfmu, c(theta, mu))
# d1mu(x, theta, mu, tau)
# d1theta_mu(x, theta, mu, tau)
# 
# hessian(lnfmu, c(theta, mu))
# d2mu(x, theta, mu, tau)
# d2theta_mu(x, theta, mu, tau)
# dmudtheta(x, theta, mu, tau)})
# 
