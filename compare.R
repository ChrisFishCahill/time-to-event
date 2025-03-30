# -----------------------------------------------------
# Exploring the relationship between Weibull AFT and log-hazard models:
# -----------------------------------------------------
#
# **Weibull AFT Model:**
# The AFT (Accelerated Failure Time) model is about directly modeling the 
# survival time. Here, we assume that the survival times follow a Weibull 
# distribution and we model the **time-to-event** directly.
#
# The AFT model works by transforming the survival times, typically by taking 
# the logarithm, and fitting the following relationship:
# 
#   log(T) = X * beta + epsilon
#
# Here, **T** is the time-to-event, **X** is the covariate matrix (which we’re 
# ignoring in this case), **beta** is the regression coefficient (effect of 
# covariates), and **epsilon** is the error term. 
#
# For the Weibull distribution, the survival function is:
# 
#   S(t) = exp(-(t / lambda)^k)
#
# The **scale parameter (lambda)** controls the time scale, while the **shape 
# parameter (k)** determines the rate of failure—whether events are more likely 
# to happen earlier or later in time.
#
# So, the AFT model estimates these two parameters **k** and **lambda** directly 
# by fitting the data to the Weibull distribution.

# **Log-Hazard Model:**
# The log-hazard model is a bit more indirect. Instead of modeling the survival 
# time directly, this model focuses on the **hazard function**—the rate of failure 
# at a given point in time. We can derive the log-hazard function for a Weibull 
# distribution from the following hazard and cumulative hazard functions:
#
#   h(t) = (k / lambda) * (t / lambda)^(k - 1)
#   H(t) = (t / lambda)^k
#
# The log of the hazard function gives us:
#
#   log(h(t)) = log(k) - log(lambda) + (k - 1) * log(t)
#
# In this model, we’re directly working with the **log(hazard)** as a linear function 
# of **log(t)**, where **log(k)** represents the intercept and **(k-1)** is the 
# coefficient (slope). So, here **k** is derived from the slope, and **lambda** 
# is tied to the intercept.
#
# **Key Difference:**
# The AFT and log-hazard models are estimating the same **underlying Weibull distribution** 
# but approach the data differently. The AFT model focuses on modeling the **time-to-event** 
# directly, whereas the log-hazard model focuses on the **failure rate** (hazard) over time.
#
# In short: AFT = "How long until the event?" while log-hazard = "How likely is the event 
# at this point in time?"

# To prove this, this script simulates weibull data, estimates it using 
# survreg, a custom optim function specifying an AFT weibull, and a custom
# optim function specifying a log hazard weibull. I then transform (i.e., use math)
# to show that these all recover identical estimates of the underlying parameters of 
# interest used to simulate data

#----------------------
# simulate Weibull data and compare accelerated failure time models 
# (i.e., survreg, optim) vs. a log-hazard model
set.seed(1)
library(survival)

# --- settings ---
add_censoring <- FALSE # set to TRUE to deal with nightmare of censoring
n_simulations <- 1000 # number of simulations
n <- 20 # large sample size per simulation
true_k <- 1.5 # true shape parameter
true_lambda <- 2 # true scale parameter

# storage for results
lambda_values_surv <- numeric(n_simulations)
lambda_values_opt <- numeric(n_simulations)
lambda_values_loghaz <- numeric(n_simulations)
k_values_surv <- numeric(n_simulations)
k_values_opt <- numeric(n_simulations)
k_values_loghaz <- numeric(n_simulations)

# --- define likelihood functions outside the simulation loop ---

# survreg AFT Weibull fit
fit_surv <- function(tte, status) {
  fit <- survreg(Surv(tte, status) ~ 1, dist = "weibull")
  mu <- coef(fit) # survreg returns -log(lambda)
  sigma <- fit$scale # sigma = 1/k
  k <- 1 / sigma
  lambda <- exp(mu) # correct transformation for scale (lambda)
  return(c(k, lambda))
}

# custom AFT Weibull likelihood with censoring
nll_weibull_aft <- function(par, t, d) {
  mu <- par[1]
  logsigma <- par[2]
  sigma <- exp(logsigma)
  k <- 1 / sigma
  lambda <- exp(mu)
  logh <- log(k) + k * log(lambda) + (k - 1) * log(t)
  H <- (lambda * t)^k
  ll <- d * logh - H
  return(-sum(ll)) # negative log-likelihood
}

# log-hazard likelihood directly in log(k), log(lambda)
nll_loghaz <- function(par, t, d) {
  a <- par[1] # intercept
  b <- par[2] # slope
  logt <- log(t)
  eta <- a + b * logt
  h <- exp(eta)
  H <- exp(a) * t^(b + 1) / (b + 1)
  ll <- d * log(h) - H
  return(-sum(ll)) # negative log-likelihood
}

# --- simulation loop ---
for (i in 1:n_simulations) {
  # simulate weibull data
  tte <- rweibull(n, shape = true_k, scale = true_lambda)
  status <- rep(1, n) # no censoring in this case

  # --- optional censoring ---
  if (add_censoring) {
    censor_time <- quantile(tte, 0.7)
    status <- ifelse(tte <= censor_time, 1, 0)
    tte <- pmin(tte, censor_time)
  }

  # --- survreg AFT Weibull fit ---
  fit_surv_results <- fit_surv(tte, status)
  k_surv <- fit_surv_results[1]
  lambda_surv <- fit_surv_results[2]

  # --- custom AFT Weibull likelihood with censoring ---
  init <- c(log(lambda_surv), log(1 / k_surv)) # initial guesses
  fit_opt <- optim(init, nll_weibull_aft, t = tte, d = status, hessian = TRUE)
  mu_opt <- fit_opt$par[1]
  sigma_opt <- exp(fit_opt$par[2])
  k_opt <- 1 / sigma_opt
  lambda_opt <- exp(-mu_opt) # correct transformation for lambda

  # --- log-hazard likelihood directly in log(k), log(lambda) ---
  fit_loghaz <- optim(
    par = c(log(k_surv) - log(lambda_surv), k_surv - 1),
    fn = nll_loghaz,
    t = tte, d = status,
    hessian = TRUE
  )
  a_loghaz <- fit_loghaz$par[1]
  b_loghaz <- fit_loghaz$par[2]
  k_loghaz <- b_loghaz + 1
  lambda_loghaz <- exp((log(k_loghaz) - a_loghaz) / k_loghaz)

  # store results
  lambda_values_surv[i] <- lambda_surv
  lambda_values_opt[i] <- lambda_opt
  lambda_values_loghaz[i] <- lambda_loghaz
  k_values_surv[i] <- k_surv
  k_values_opt[i] <- k_opt
  k_values_loghaz[i] <- k_loghaz
}

# --- compare distributions of lambda and k ---
par(mfrow = c(3, 2), mar = c(5, 4, 4, 2))

# lambda comparison (survreg)
plot(density(lambda_values_surv),
  col = rgb(0, 0, 1, 0.5),
  main = "lambda estimates (survreg AFT)",
  xlab = "lambda", ylab = "density"
)
abline(v = true_lambda, col = "black", lty = 2, lwd = 3)

# shape (k) comparison (survreg)
plot(density(k_values_surv),
  col = rgb(0, 0, 1, 0.5),
  main = "shape (k) estimates (survreg AFT)",
  xlab = "shape (k)", ylab = "density"
)
abline(v = true_k, col = "black", lty = 2, lwd = 3)

# lambda comparison (optim)
plot(density(lambda_values_opt),
  col = rgb(0, 1, 0, 0.5),
  main = "lambda estimates (optim AFT)",
  xlab = "lambda", ylab = "density"
)
abline(v = true_lambda, col = "black", lty = 2, lwd = 3)

# shape (k) comparison (optim)
plot(density(k_values_opt),
  col = rgb(0, 1, 0, 0.5),
  main = "shape (k) estimates (optim AFT)",
  xlab = "shape (k)", ylab = "density"
)
abline(v = true_k, col = "black", lty = 2, lwd = 3)

# lambda comparison (loghaz)
plot(density(lambda_values_loghaz),
  col = rgb(1, 0, 0, 0.5),
  main = "lambda estimates (optim loghaz)",
  xlab = "lambda", ylab = "density"
)
abline(v = true_lambda, col = "black", lty = 2, lwd = 3)

# shape (k) comparison (loghaz)
plot(density(k_values_loghaz),
  col = rgb(1, 0, 0, 0.5),
  main = "shape (k) estimates (optim loghaz)",
  xlab = "shape (k)", ylab = "density"
)
abline(v = true_k, col = "black", lty = 2, lwd = 3)
