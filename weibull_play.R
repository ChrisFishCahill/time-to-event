# weibull density function:
#
# f(t; k, lambda) = (k / lambda) * (t / lambda)^(k - 1) * exp(-(t / lambda)^k)
#
# Where:
#   - t       = time to event (e.g., when a fish exits the river section)
#   - k       = shape parameter (controls how the risk changes over time)
#   - lambda  = scale parameter (controls the overall timing of the event)
#
# Interpretation:
#   - If k = 1: exponential (memoryless, constant hazard)
#   - If k < 1: high early risk, then decreasing hazard (fish move early)
#   - If k > 1: low early risk, then increasing hazard (fish accelerate movement)
#
# Note:
# could model log(lambda) as a function of stuff, just like a GLM

# start by visualizing some stuff
t <- seq(0, 10, length.out = 200)
k <- c(0.7, 1, 1.3)
lambda <- 2
hazards <- lapply(k, function(ki) (ki / lambda) * (t / lambda)^(ki - 1))
par(mfrow = c(1, 1))
matplot(t, do.call(cbind, hazards),
  type = "l", lty = 1, col = 2:4,
  ylab = "Hazard", xlab = "Time", main = "Weibull Hazards"
)
legend("topright", legend = paste("k =", k), col = 2:4, lty = 1)

# simulate some weibull tte data to visualize
set.seed(123)
n <- 1000
lambda <- 1 # scale parameter

k_vals <- c(0.5, 1, 1.5) # different movement styles
sim_tte <- lapply(k_vals, function(k) rweibull(n, shape = k, scale = lambda))

# function to compute hazard, survival, pdf
get_curves <- function(t, k, lambda) {
  h <- (k / lambda) * (t / lambda)^(k - 1)
  S <- exp(-(t / lambda)^k)
  f <- h * S
  list(h = h, S = S, f = f)
}

# time vector
t <- seq(0.01, 5, length.out = 200)

# compute curves for each k
curves <- lapply(k_vals, function(k) get_curves(t, k, lambda))
names(curves) <- paste0("k=", k_vals)

# plotting
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for (i in seq_along(curves)) {
  crv <- curves[[i]]
  hist(sim_tte[[i]],
    breaks = 30, probability = TRUE,
    main = paste("Weibull PDF (k =", k_vals[i], ")"), xlab = "t"
  )
  lines(t, crv$f, col = "blue", lwd = 2)
  plot(t, crv$S,
    type = "l", main = paste("Survival S(t)"),
    ylab = "", xlab = "t"
  )
  plot(t, crv$h,
    type = "l", main = paste("Hazard h(t)"),
    ylab = "", xlab = "t"
  )
}

#----------------------
# simulate some data and estimate using survreg and optim
set.seed(765)
library(survival)

# --- simulate Data ---
n_fish <- 30
lambda <- 1
true_k <- c(Control = 1, Treatment = 1.5)

tte_control <- rweibull(n_fish, shape = true_k["Control"], scale = lambda)
tte_treat <- rweibull(n_fish, shape = true_k["Treatment"], scale = lambda)

df <- data.frame(
  tte = c(tte_control, tte_treat),
  group = rep(c("Control", "Treatment"), each = n_fish)
)

# --- survreg() models (separate fits each group) ---
fits <- by(df, df$group, function(sub) {
  survreg(Surv(tte) ~ 1, data = sub, dist = "weibull")
})

extract_ci <- function(fit) {
  scale_hat <- fit$scale
  se <- fit$scale * fit$var[1]^0.5
  scale_ci <- scale_hat + c(-1.96, 1.96) * se
  shape_est <- 1 / scale_hat
  shape_ci <- 1 / rev(scale_ci)
  c(est = shape_est, lower = shape_ci[1], upper = shape_ci[2])
}
params <- t(sapply(fits, extract_ci))
rownames(params) <- names(true_k)

# --- optim() separate fits each group---
nll_opt <- function(par, t) {
  log_k <- par[1]
  log_lambda <- par[2]
  k <- exp(log_k)
  lambda <- exp(log_lambda)
  # note could write out likelihood as
  # nll <- -sum(log(k) - log(lambda) + (k - 1) * log(tte_control / lambda) -
  #  (tte_control / lambda)^k)
  nll <- -sum(dweibull(t, shape = k, scale = lambda, log = TRUE))
  nll
}
fit_ctrl <- optim(c(0, 0), nll_opt, t = tte_control, hessian = TRUE)
fit_trt <- optim(c(0, 0), nll_opt, t = tte_treat, hessian = TRUE)

k_ctrl <- exp(fit_ctrl$par[1])
k_trt <- exp(fit_trt$par[1])

# --- optim() joint model for both groups ---
nll_joint <- function(par, t, g) {
  log_k_ctrl <- par[1]
  log_k_trt <- par[2]
  log_lambda_ctrl <- par[3]
  log_lambda_trt <- par[4]
  k <- ifelse(g == "Control", exp(log_k_ctrl), exp(log_k_trt))
  # note could write out likelihood as
  # nll <- -sum(log(k) - log(lambda) + (k - 1) * log(tte_control / lambda) -
  #  (tte_control / lambda)^k)
  lambda <- ifelse(g == "Control", exp(log_lambda_ctrl), exp(log_lambda_trt))
  nll <- -sum(dweibull(t, shape = k, scale = lambda, log = TRUE))
  nll
}
init <- rep(log(1), 4)
fit_joint <- optim(
  par = init,
  fn = function(par) nll_joint(par, t = df$tte, g = df$group),
  hessian = TRUE
)
k_joint <- exp(fit_joint$par[1:2])
names(k_joint) <- c("Control", "Treatment")

# --- optim() group-specific log-hazard model ---
# Model:
# if Control:    log h(t) = a0 + b0 * log t
# if Treatment:  log h(t) = a1 + b1 * log t
# 
# recognize that 
# group-specific log-hazard model:
#   log h(t) = a_g + b_g * log(t)
#
# this is equivalent to the Weibull hazard:
#   h(t) = (k / lambda) * (t / lambda)^(k - 1)
#
# taking logs:
#   log h(t) = log(k) - log(lambda) + (k - 1) * log(t)
#
# match coefficients:
#   b_g     = k_g - 1 --> k_g = b_g + 1
#   a_g     = log(k_g) - log(lambda_g)
#
# solve for lambda_g:
#   lambda_g = exp(log(k_g) - a_g)

nll_loghaz_groups <- function(par) {
  a0 <- par[1] # intercept for Control
  b0 <- par[2] # slope for Control
  a1 <- par[3] # intercept for Treatment
  b1 <- par[4] # slope for Treatment
  logt <- log(df$tte)
  is_trt <- df$group == "Treatment"
  a <- ifelse(is_trt, a1, a0)
  b <- ifelse(is_trt, b1, b0)
  eta <- a + b * logt
  h <- exp(eta)
  H <- exp(eta - b * logt) * df$tte^(b + 1) / (b + 1)
  nll <- -sum(log(h) - H)
  nll
}

init <- rep(0, 4)
fit_loghaz_groups <- optim(init, nll_loghaz_groups, hessian = TRUE)

# extract
pars <- fit_loghaz_groups$par
a0 <- pars[1]
b0 <- pars[2]
a1 <- pars[3]
b1 <- pars[4]

# recover Weibull k and lambda
k0 <- b0 + 1
lambda0 <- exp(log(k0) - a0)

k1 <- b1 + 1
lambda1 <- exp(log(k1) - a1)

# --- plotting ---
ks_all <- data.frame(
  Control = c(
    params["Control", "est"],
    k_ctrl, k_joint["Control"], k0, true_k["Control"]
  ),
  Treatment = c(
    params["Treatment", "est"],
    k_trt, k_joint["Treatment"], k1, true_k["Treatment"]
  )
)
rownames(ks_all) <- c(
  "survreg", "optim_sep",
  "optim_joint", "loghaz_grp_optim", "truth"
)
par(mfrow=c(1,1))
plot(NULL,
  xlim = c(0.75, 2.25), ylim = c(0.5, 2.4),
  xlab = "", ylab = "Shape (k)", xaxt = "n",
  main = "Weibull Shape Parameter Estimates"
)
axis(1, at = 1:2, labels = c("Control", "Treatment"))

# CI bars from survreg
arrows(1:2, params[, "lower"], 1:2, params[, "upper"],
  angle = 90, code = 0, length = 0.1
)

# points for all models
pch_vec <- c(19, 17, 15, 3, 4) # updated symbols
col_vec <- c("blue", "darkgreen", "purple", "black", "red")

for (i in 1:nrow(ks_all)) {
  points(1:2, ks_all[i, ], pch = pch_vec[i], col = col_vec[i])
}

legend("topleft",
  legend = rownames(ks_all),
  pch = pch_vec,
  col = col_vec,
  title = "Estimator"
)
