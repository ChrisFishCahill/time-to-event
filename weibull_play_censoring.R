#----------------------
# simulate some data with censoring and estimate using survreg and log-hazard
set.seed(765)
library(survival)

# --- simulate Data ---
n_fish <- 30
lambda <- 1
true_k <- c(Control = 1, Treatment = 1.5)

# Simulate TTE data
tte_control <- rweibull(n_fish,
  shape = true_k["Control"],
  scale = lambda
)
tte_treat <- rweibull(n_fish,
  shape = true_k["Treatment"],
  scale = lambda
)
tte <- c(tte_control, tte_treat)

# Simulate censoring
censor_time <- 2
event <- as.integer(tte <= censor_time)
tte_obs <- pmin(tte, censor_time)

df <- data.frame(
  tte = tte_obs,
  event = event,
  group = rep(c("Control", "Treatment"), each = n_fish)
)

# --- survreg() models (separate fits each group) ---
fits <- by(df, df$group, function(sub) {
  survreg(Surv(tte, event) ~ 1, data = sub, dist = "weibull")
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

# --- log-hazard group-specific model with censoring ---
nll_loghaz_groups <- function(par) {
  a0 <- par[1]
  b0 <- par[2]
  a1 <- par[3]
  b1 <- par[4]
  logt <- log(df$tte)
  is_trt <- df$group == "Treatment"
  event <- df$event
  a <- ifelse(is_trt, a1, a0)
  b <- ifelse(is_trt, b1, b0)
  eta <- a + b * logt
  h <- exp(eta)
  H <- exp(eta - b * logt) * df$tte^(b + 1) / (b + 1)
  nll <- -sum(event * log(h) - H)
  nll
}

fit_loghaz_groups <- optim(rep(0, 4), nll_loghaz_groups, hessian = TRUE)
pars <- fit_loghaz_groups$par

# recover Weibull k and lambda
k0 <- pars[2] + 1
lambda0 <- exp(log(k0) - pars[1])
k1 <- pars[4] + 1
lambda1 <- exp(log(k1) - pars[3])

# --- plotting ---
ks_all <- data.frame(
  Control = c(
    params["Control", "est"], k0, true_k["Control"]
  ),
  Treatment = c(
    params["Treatment", "est"], k1, true_k["Treatment"]
  )
)
rownames(ks_all) <- c("survreg", "loghaz_grp_optim", "truth")

par(mfrow = c(1, 1))
plot(NULL,
  xlim = c(0.75, 2.25), ylim = c(0.5, 2.4),
  xlab = "", ylab = "Shape (k)", xaxt = "n",
  main = "Weibull Shape Parameter Estimates with Censoring"
)
axis(1, at = 1:2, labels = c("Control", "Treatment"))

arrows(1:2, params[, "lower"], 1:2, params[, "upper"],
  angle = 90, code = 0, length = 0.1
)

pch_vec <- c(19, 3, 4)
col_vec <- c("blue", "black", "red")

for (i in 1:nrow(ks_all)) {
  points(1:2, ks_all[i, ], pch = pch_vec[i], col = col_vec[i])
}

legend("topleft",
  legend = rownames(ks_all),
  pch = pch_vec,
  col = col_vec,
  title = "Estimator"
)
