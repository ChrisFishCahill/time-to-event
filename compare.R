# -----------------------------------------------------
# Exploring the relationship between Weibull AFT and log-hazard models:
# -----------------------------------------------------
#
# This script compares four equivalent parameterizations of Weibull models:
#
# 1. AFT model using survreg:
#    - Models log(T) = X * beta + error, where error ~ Extreme Value
#    - survreg internally parameterizes Weibull with:
#        scale = 1 / k, intercept = -log(lambda)
#    - We recover:
#        k = 1 / scale
#        lambda = exp(mu)
#
# 2. AFT model using optim:
#    - Maximizes the log-likelihood directly under the Weibull AFT formulation.
#    - Estimates mu = log(scale) and logsigma = log(1 / shape)
#    - Recover:
#        k = 1 / exp(logsigma)
#        lambda = exp(mu)
#
# 3. Log-hazard model using (a, b):
#    - Starts from the Weibull hazard function:
#        h(t) = (k / lambda) * (t / lambda)^(k - 1)
#      which expands to:
#        log(h(t)) = log(k) - log(lambda) + (k - 1) * log(t)
#    - We model:
#        log(h(t)) = a + b * log(t)
#      and recover:
#        k = b + 1
#        lambda = exp((log(k) - a) / k)
#
# 4. Log-hazard model using (log(shape), log(scale)):
#    - Parameterize directly in terms of log(k) and log(lambda)
#    - From these, derive:
#        a = log(k) - k * log(lambda)
#        b = k - 1
#
# Censoring and the Log-Likelihood:
# ---------------------------------
# Let d_i = 1 if event observed, 0 if right-censored.
# Define:
#   h(t) = hazard function
#   H(t) = cumulative hazard function = ∫₀^t h(s) ds
#
# The full log-likelihood is:
#   ℓ = Σ_i [ d_i * log(h(t_i)) - H(t_i) ]
#
# For the Weibull:
#   h(t) = (k / lambda) * (t / lambda)^(k - 1)
#   H(t) = (t / lambda)^k
#
# So:
#   ℓ = Σ_i [ d_i * (log(k) - log(lambda) + (k - 1) * log(t_i)) - (t_i / lambda)^k ]
#
# All four approaches optimize this same expression, just reparameterized.
# With or without censoring, they yield the same estimates when done correctly.

set.seed(1)
library(survival)

add_censoring <- FALSE
n_simulations <- 1000
n <- 20
true_k <- 1.5
true_lambda <- 2

lambda_values_surv <- numeric(n_simulations)
lambda_values_opt <- numeric(n_simulations)
lambda_values_loghaz <- numeric(n_simulations)
lambda_values_loghaz2 <- numeric(n_simulations)
k_values_surv <- numeric(n_simulations)
k_values_opt <- numeric(n_simulations)
k_values_loghaz <- numeric(n_simulations)
k_values_loghaz2 <- numeric(n_simulations)

fit_surv <- function(tte, status) {
  fit <- survreg(Surv(tte, status) ~ 1, dist = "weibull")
  mu <- coef(fit)
  sigma <- fit$scale
  k <- 1 / sigma
  lambda <- exp(mu)
  return(c(k, lambda))
}

nll_weibull_aft <- function(par, t, d) {
  mu <- par[1]
  logsigma <- par[2]
  sigma <- exp(logsigma)
  k <- 1 / sigma
  lambda <- exp(mu)
  logh <- log(k) + k * log(lambda) + (k - 1) * log(t)
  H <- (lambda * t)^k
  ll <- d * logh - H
  return(-sum(ll))
}

nll_loghaz <- function(par, t, d) {
  a <- par[1]
  b <- par[2]
  logt <- log(t)
  eta <- a + b * logt
  h <- exp(eta)
  H <- exp(a) * t^(b + 1) / (b + 1)
  ll <- d * log(h) - H
  return(-sum(ll))
}

nll_loghaz_shape_scale <- function(par, t, d) {
  k <- exp(par[1])
  lambda <- exp(par[2])
  a <- log(k) - k * log(lambda)
  b <- k - 1
  logt <- log(t)
  eta <- a + b * logt
  h <- exp(eta)
  H <- exp(a) * t^(b + 1) / (b + 1)
  ll <- d * log(h) - H
  return(-sum(ll))
}

for (i in 1:n_simulations) {
  tte <- rweibull(n, shape = true_k, scale = true_lambda)
  status <- rep(1, n)
  
  if (add_censoring) {
    censor_time <- quantile(tte, 0.7)
    status <- ifelse(tte <= censor_time, 1, 0)
    tte <- pmin(tte, censor_time)
  }
  
  fit_surv_results <- fit_surv(tte, status)
  k_surv <- fit_surv_results[1]
  lambda_surv <- fit_surv_results[2]
  
  init <- c(log(lambda_surv), log(1 / k_surv))
  fit_opt <- optim(init, nll_weibull_aft, t = tte, d = status)
  mu_opt <- fit_opt$par[1]
  sigma_opt <- exp(fit_opt$par[2])
  k_opt <- 1 / sigma_opt
  lambda_opt <- exp(-mu_opt)
  
  fit_loghaz <- optim(
    par = c(log(k_surv) - log(lambda_surv), k_surv - 1),
    fn = nll_loghaz,
    t = tte, d = status
  )
  a_loghaz <- fit_loghaz$par[1]
  b_loghaz <- fit_loghaz$par[2]
  k_loghaz <- b_loghaz + 1
  lambda_loghaz <- exp((log(k_loghaz) - a_loghaz) / k_loghaz)
  
  fit_loghaz2 <- optim(c(log(true_k), log(true_lambda)), 
                       nll_loghaz_shape_scale, t = tte, d = status)
  k_loghaz2 <- exp(fit_loghaz2$par[1])
  lambda_loghaz2 <- exp(fit_loghaz2$par[2])
  
  lambda_values_surv[i] <- lambda_surv
  lambda_values_opt[i] <- lambda_opt
  lambda_values_loghaz[i] <- lambda_loghaz
  lambda_values_loghaz2[i] <- lambda_loghaz2
  k_values_surv[i] <- k_surv
  k_values_opt[i] <- k_opt
  k_values_loghaz[i] <- k_loghaz
  k_values_loghaz2[i] <- k_loghaz2
}

par(mfrow = c(4, 2), mar = c(5, 4, 4, 2))

plot(density(lambda_values_surv), col = rgb(0, 0, 1, 0.5),
     main = "lambda (survreg AFT)", xlab = "lambda")
abline(v = true_lambda, col = "black", lty = 2, lwd = 2)

plot(density(k_values_surv), col = rgb(0, 0, 1, 0.5),
     main = "shape k (survreg AFT)", xlab = "k")
abline(v = true_k, col = "black", lty = 2, lwd = 2)

plot(density(lambda_values_opt), col = rgb(0, 1, 0, 0.5),
     main = "lambda (optim AFT)", xlab = "lambda")
abline(v = true_lambda, col = "black", lty = 2, lwd = 2)

plot(density(k_values_opt), col = rgb(0, 1, 0, 0.5),
     main = "shape k (optim AFT)", xlab = "k")
abline(v = true_k, col = "black", lty = 2, lwd = 2)

plot(density(lambda_values_loghaz), col = rgb(1, 0, 0, 0.5),
     main = "lambda (log-hazard ab)", xlab = "lambda")
abline(v = true_lambda, col = "black", lty = 2, lwd = 2)

plot(density(k_values_loghaz), col = rgb(1, 0, 0, 0.5),
     main = "shape k (log-hazard ab)", xlab = "k")
abline(v = true_k, col = "black", lty = 2, lwd = 2)

plot(density(lambda_values_loghaz2), col = rgb(1, 0.4, 0.4, 0.5),
     main = "lambda (log-hazard shape/scale)", xlab = "lambda")
abline(v = true_lambda, col = "black", lty = 2, lwd = 2)

plot(density(k_values_loghaz2), col = rgb(1, 0.4, 0.4, 0.5),
     main = "shape k (log-hazard shape/scale)", xlab = "k")
abline(v = true_k, col = "black", lty = 2, lwd = 2)