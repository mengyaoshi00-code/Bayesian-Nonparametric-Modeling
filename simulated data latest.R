## =========================
## Weibull–AFT 
## =========================
#library
library(splines)
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(LaplacesDemon)  # rdirichlet

normalize <- function(x) x / sum(x)

weibull_S_from_mu <- function(t, mu, shape) {
  # Survival S(t | mu, k) with parameterization E[T] = mu
  # Using lambda = mu / Gamma(1 + 1/k)
  exp(- ((t * gamma(1 + 1/shape)) / mu)^shape )
}

true_tau_pattern <- function(K) {
  g_seq <- (1:K) / K * 2 * pi
  0.5 * sin(g_seq) + 0.1 * cos(2 * g_seq)
}

# --- main simulator ---
simulate_aft_weibull <- function(n = 2000, K = 8, seed = 123,
                                 shape_event = 1.5, shape_cens = 1.2) {
  set.seed(seed)
  
  ## 1) Subgroups G (unbalanced)
  probs <- normalize(exp(-0.3 * (0:(K-1))))
  G <- sample.int(K, n, replace = TRUE, prob = probs)
  
  ## 2) Covariates L1, L2 with group-specific patterns (then scale)
  L1 <- numeric(n); L2 <- numeric(n)
  split_idx <- split(seq_len(n), G)
  for (g in seq_len(K)) {  # LOOP: fill L1/L2 pattern for each subgroup g
    idx <- split_idx[[g]]
    m  <- length(idx)
    if (m == 0) next
    freq  <- 1 + 0.3 * (g - 1)
    Xg    <- seq(0, 2*pi*freq, length.out = m)
    X01   <- (Xg - min(Xg)) / (max(Xg) - min(Xg))
    osc   <- exp(-0.5 * Xg) * cos(2*pi*Xg)
    L1[idx] <- (X01)*osc + rev(X01)*rev(osc) + rnorm(m, 0, 0.1)
    phase   <- pi/4 * (g - 1)
    L2[idx] <- sin(2*pi*Xg + phase) * exp(-0.3*Xg) + rnorm(m, 0, 0.1)
  }
  L1 <- as.numeric(scale(L1))
  L2 <- as.numeric(scale(L2))
  
  ## 3) Binary covariate L3 with mild trend
  L3_prob <- 0.5 + 0.3 * sin(seq(0, 2*pi, length.out = n))
  L3 <- rbinom(n, 1, pmin(0.9, pmax(0.1, L3_prob)))
  
  ## 4) Treatment A (confounded by L)
  propensity <- plogis(0.15*L1 + 0.15*L2 + 0.1*L3)
  A <- rbinom(n, 1, propensity)
  
  ## 5) True subgroup effects tau_g
  tau_g <- true_tau_pattern(K)
  
  ## 6) Event times T via Weibull–AFT (E[T] = exp(eta))
  eta  <- 1.5 + 0.3*L1 + 0.4*L2 + 0.3*L3 + tau_g[G]*A
  mu   <- exp(eta)
  scale_event <- mu / gamma(1 + 1/shape_event)
  T_event <- rweibull(n, shape = shape_event, scale = scale_event)
  
  ## 7) Informative censoring C (depends on A and L)
  eta_c <- 2.2 + 0.15*L1 + 0.15*L2 - 0.2*A + 0.1*L3
  mu_c  <- exp(eta_c)
  scale_cens <- mu_c / gamma(1 + 1/shape_cens)
  C <- rweibull(n, shape = shape_cens, scale = scale_cens)
  
  ## 8) Observed outcomes (right-censoring)
  time  <- pmin(T_event, C)
  delta <- as.integer(T_event <= C)
  
  ## 9) 1-year binary endpoints
  t1 <- 1.0
  Y1_true <- as.integer(T_event <= t1)              # standard
  Y1_obs  <- rep(NA_integer_, n)                    # observed with missingness
  Y1_obs[delta == 1 & time <= t1] <- 1
  Y1_obs[time  >  t1]               <- 0            # known not failed by 1y
  
  ## 10) Return
  data <- data.frame(
    id = seq_len(n), time, delta, A, L1, L2, L3, G,
    true_tau = tau_g[G], Y1_true, Y1_obs
  )
  attr(data, "tau_g")        <- tau_g
  attr(data, "shape_event")  <- shape_event
  attr(data, "shape_cens")   <- shape_cens
  data
}


# For each subgroup g, it calculates the true risk difference at t =1 year

compute_true_RD1y <- function(dat, t1 = 1.0) {
  shape <- attr(dat, "shape_event")
  tau_g <- attr(dat, "tau_g")
  stopifnot(!is.null(shape), !is.null(tau_g))
  
  eta0 <- 1.5 + 0.3*dat$L1 + 0.4*dat$L2 + 0.3*dat$L3           # control
  mu0  <- exp(eta0)
  mu1  <- exp(eta0 + tau_g[dat$G])                             # treated
  
  S0 <- weibull_S_from_mu(t1, mu0, shape)
  S1 <- weibull_S_from_mu(t1, mu1, shape)
  
  by_g <- split(seq_len(nrow(dat)), dat$G)
  
  # LOOP: average 1y risk diff within each subgroup g
  RD1y <- sapply(by_g, function(idx) mean(S0[idx]) - mean(S1[idx]))
  data.frame(Subgroup = as.integer(names(by_g)), RD_1y_true = as.numeric(RD1y))
}

dat <- simulate_aft_weibull(n = 2000, K = 8, seed = 123)


# True subgroup risk differences at 1 year (from the same AFT parameters)
rd_table <- compute_true_RD1y(dat, t1 = 1.0)
rd_table


# =========================================================
# Stan (parametric AFT–Weibull with informative censoring)
# =========================================================
# --- bases for smooth f1(L1), f2(L2) ---

prep_np_data <- function(dat, K_mix = 3, df_spline = 6) {
  X1 <- bs(dat$L1, df = df_spline, intercept = FALSE)
  X2 <- bs(dat$L2, df = df_spline, intercept = FALSE)
  
  list(
    N = nrow(dat),
    K = length(unique(dat$G)),
    J1 = ncol(X1),
    J2 = ncol(X2),
    time = as.vector(dat$time),
    delta = as.integer(dat$delta),
    L3 = as.vector(dat$L3),
    A  = as.vector(dat$A),
    G  = as.integer(dat$G),
    X1 = X1, X2 = X2,
    K_mix = K_mix,                 # mixture components for log-time error
    tau_prior_sd = 1.0
  )
}
stan_data_np <- prep_np_data(dat, K_mix = 3, df_spline = 6)

stan_parametric_ic <- "
data {
  int<lower=1> N;                           // obs
  int<lower=1> K;                           // subgroups
  vector[N] time;                            // observed = min(T, C), must be > 0
  array[N] int<lower=0,upper=1> delta;       // 1=event observed, 0=censored

  vector[N] L1;
  vector[N] L2;
  vector[N] L3;
  vector[N] A;
  array[N] int<lower=1,upper=K> G;

  real<lower=0> tau_prior_sd;                // kept for interface compatibility
}
transformed data {
  int Ne = 0;
  int Nc = 0;
  array[N] int event_ind;
  array[N] int cens_ind;

  // Centering for stability 
  real mean_L1 = mean(L1);
  real mean_L2 = mean(L2);
  real mean_L3 = mean(L3);
  real mean_A  = mean(A);

  vector[N] L1c = L1 - mean_L1;
  vector[N] L2c = L2 - mean_L2;
  vector[N] L3c = L3 - mean_L3;
  vector[N] Ac  = A  - mean_A;

  for (i in 1:N) {
  // LOOP: split indices into event vs censored
    if (delta[i] == 1) { Ne += 1; event_ind[Ne] = i; }
    else               { Nc += 1;  cens_ind[Nc]  = i; }
  }
}
parameters {
  // Event model (AFT–Weibull on mean-time scale)
  real Intercept;
  real bL1;
  real bL2;
  real bL3;

  // subgroup baselines with shrinkage
  vector[K] bG_raw;
  real<lower=0> sigma_G;

  // subgroup HTE tau_g with shrinkage
  real mu_tau;
  real<lower=0> sigma_tau;
  vector[K] tau_g_raw;

  // Weibull shape for event (bounded away from 0)
  real<lower=0.2> shape_T;

  // Censoring model (Weibull on mean-time scale)
  real c_Intercept;
  real cL1;
  real cL2;
  real cL3;
  real cA;
  real<lower=0.2> shape_C;
}
transformed parameters {
  vector[K] bG   = sigma_G * bG_raw;                  // subgroup baselines
  vector[K] tau_g = mu_tau + sigma_tau * tau_g_raw;   // subgroup HTE

  vector[N] eta_T;
  vector[N] eta_C;
  for (i in 1:N) {
  // LOOP: build linear predictors for event and censoring models
    eta_T[i] = Intercept
               + bL1 * L1c[i] + bL2 * L2c[i] + bL3 * L3c[i]
               + bG[G[i]]
               + tau_g[G[i]] * Ac[i];   // subgroup-specific treatment effect
    eta_C[i] = c_Intercept
               + cL1 * L1c[i] + cL2 * L2c[i] + cL3 * L3c[i]
               + cA  * Ac[i];          // censoring depends on A,L; no G
  }
}
model {
  // Priors (weakly-informative)
  Intercept ~ normal(1.5, 2.5);
  bL1 ~ normal(0, 3);
  bL2 ~ normal(0, 3);
  bL3 ~ normal(0, 3);

  sigma_G ~ normal(0, 0.7);             // half-normal via <lower=0>
  bG_raw  ~ std_normal();

  mu_tau   ~ normal(0, 1);
  sigma_tau ~ normal(0, 0.7);           // half-normal
  tau_g_raw ~ std_normal();

  // Shapes as lognormal around reasonable centers
  shape_T ~ lognormal(log(1.5), 0.3);
  shape_C ~ lognormal(log(1.2), 0.3);

  // Censoring coefficients
  c_Intercept ~ normal(2.0, 2.5);
  cL1 ~ normal(0, 2);
  cL2 ~ normal(0, 2);
  cL3 ~ normal(0, 2);
  cA  ~ normal(0, 2);

  // --- Joint likelihood ---
  // Stan uses Weibull(shape, scale). Our AFT mean-time is mu = exp(eta).
  // scale = mu / Gamma(1 + 1/shape)  -> use lgamma for stability.
  {
    vector[Ne] t_e   = time[event_ind[1:Ne]];
    vector[Ne] eta_eT = eta_T[event_ind[1:Ne]];
    vector[Ne] eta_eC = eta_C[event_ind[1:Ne]];

    vector[Nc] t_c   = time[cens_ind[1:Nc]];
    vector[Nc] eta_cT = eta_T[cens_ind[1:Nc]];
    vector[Nc] eta_cC = eta_C[cens_ind[1:Nc]];

    // Events: f_T(t) * S_C(t)
    target += weibull_lpdf( t_e
                            | shape_T,
                              exp(eta_eT - lgamma(1 + 1/shape_T)) );
    target += weibull_lccdf( t_e
                             | shape_C,
                               exp(eta_eC - lgamma(1 + 1/shape_C)) );

    // Censored: f_C(t) * S_T(t)
    target += weibull_lpdf( t_c
                            | shape_C,
                              exp(eta_cC - lgamma(1 + 1/shape_C)) );
    target += weibull_lccdf( t_c
                             | shape_T,
                               exp(eta_cT - lgamma(1 + 1/shape_T)) );
  }
}
generated quantities {
  // de-centered intercept (undo centering)
  real b_Intercept =
    Intercept - (mean_L1*bL1 + mean_L2*bL2 + mean_L3*bL3 + mean_A*mu_tau);

  // subgroup τ_g on AFT (log-time) scale
  vector[K] subgroup_cate = tau_g;
}
"



# ---- Data ----
prepare_stan_data <- function(dat) {
  lst <- list(
    N = nrow(dat),
    K = length(unique(dat$G)),
    time  = as.vector(dat$time),
    delta = as.integer(dat$delta),
    L1 = as.vector(dat$L1),
    L2 = as.vector(dat$L2),
    L3 = as.vector(dat$L3),
    A  = as.vector(dat$A),
    G  = as.integer(dat$G),
    tau_prior_sd = 1.0
  )
  # clamp times away from zero for numerical stability
  lst$time <- pmax(lst$time, 1e-8)
  lst
}

stan_data <- prepare_stan_data(dat)


init_fun <- function() list(
  Intercept = 1.5,
  bL1 = 0, bL2 = 0, bL3 = 0,
  bG_raw = rep(0, stan_data$K),
  sigma_G = 0.3,
  
  mu_tau = 0,
  sigma_tau = 0.3,
  tau_g_raw = rep(0, stan_data$K),
  
  shape_T = 1.5,
  
  c_Intercept = 2.0,
  cL1 = 0, cL2 = 0, cL3 = 0, cA = 0,
  shape_C = 1.2
)

cat("\n=== Fitting Parametric Model (AFT–Weibull with censoring) ===\n")
fit_parametric_ic <- rstan::stan(
  model_code = stan_parametric_ic,
  data  = stan_data,
  iter  = 3000, warmup = 1500, chains = 4, seed = 123,
  init  = init_fun,
  control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.03)
)

# =========================================================
# HDP–HBB g-computation (nonparametric P_g over L)
# =========================================================


# sample one set of hierarchical Dirichlet weights for subgroup g
# returns a length-Ng vector (weights for members of g)
.sample_hdp_weights_for_g <- function(idx_g, N, alpha0 = 1, alphag = 1) {
  # global weights over all N atoms
  pi0 <- as.numeric(rdirichlet(1, rep(alpha0 / N, N)))
  # restrict and renormalize to subgroup atoms
  bg <- pi0[idx_g]
  bg <- bg / sum(bg)
  # subgroup weights
  as.numeric(rdirichlet(1, alphag * bg))
}

# HDP-HBB g-computation over p(L|G=g)
compute_subgroup_cate_hbb_hdp <- function(stan_fit, dat,
                                          timepoints = seq(0.5, 5, by = 0.5),
                                          B = 400,            # MC per person/posterior draw
                                          H = 100,            # number of HDP-HBB draws
                                          alpha0 = 1,         # global concentration
                                          alphag = 1,         # subgroup concentration
                                          subgroups_to_compute = NULL,
                                          n_draws_subsample = 1000,
                                          batch_size = 80,
                                          seed = 2026) {
  set.seed(seed)
  
  post <- rstan::extract(stan_fit)
  # Current Stan program: shape_T, Intercept, bL*, bG, subgroup_cate
  n_full <- length(post$shape_T)
  if (n_full > n_draws_subsample) {
    keep <- sample.int(n_full, n_draws_subsample)
    post <- lapply(post, function(x) if (is.matrix(x)) x[keep, , drop = FALSE] else x[keep])
  }
  n_draws <- length(post$shape_T)
  
  K  <- length(unique(dat$G))
  if (is.null(subgroups_to_compute)) subgroups_to_compute <- 1:K
  
  # match Stan centering
  mean_L1 <- mean(dat$L1); mean_L2 <- mean(dat$L2); mean_L3 <- mean(dat$L3); mean_A <- mean(dat$A)
  
  res <- vector("list", K)
  N <- nrow(dat)
  
  # LOOP: compute CATE(t) for each subgroup g via HDP–HBB
  for (g in subgroups_to_compute) {
    message(sprintf("HDP-HBB g-computation for subgroup %d/%d ...", g, K))
    dg   <- dat %>% filter(G == g)
    idxg <- which(dat$G == g)
    Ng   <- length(idxg); if (Ng == 0) next
    
    # centered covariates for subgroup rows
    L1c <- dg$L1 - mean_L1
    L2c <- dg$L2 - mean_L2
    L3c <- dg$L3 - mean_L3
    
    Inter   <- post$Intercept
    bL1     <- post$bL1
    bL2     <- post$bL2
    bL3     <- post$bL3
    bGg     <- post$bG[, g]
    tau_g   <- rstan::extract(stan_fit, pars = "subgroup_cate")$subgroup_cate[, g]
    shape_T <- post$shape_T
    
    A0c <- (0 - mean_A)
    A1c <- (1 - mean_A)
    
    nt <- length(timepoints)
    surv0 <- matrix(0, nrow = n_draws, ncol = nt)
    surv1 <- matrix(0, nrow = n_draws, ncol = nt)
    
    n_batches <- ceiling(n_draws / batch_size)
    
    # LOOP: process posterior draws in batches for memory/speed
    for (b in seq_len(n_batches)) {
      idxb <- ((b - 1) * batch_size + 1):min(b * batch_size, n_draws)
      nb   <- length(idxb)
      
      # linear predictors (Ng x nb)
      lp_base <- outer(L1c, bL1[idxb]) +
        outer(L2c, bL2[idxb]) +
        outer(L3c, bL3[idxb])
      lp_base <- sweep(lp_base, 2, Inter[idxb] + bGg[idxb], "+")
      lp0 <- sweep(lp_base, 2, tau_g[idxb] * A0c, "+")
      lp1 <- sweep(lp_base, 2, tau_g[idxb] * A1c, "+")
      
      # For each posterior draw in the batch, compute H HDP-HBB weighted survivals
      for (jl in seq_len(nb)) {
        j <- idxb[jl]; kshape <- shape_T[j]
        
        # Precompute per-person survival over the time grid (Ng x nt) under A=0 and A=1
        S0_by_i <- matrix(0, nrow = Ng, ncol = nt)
        S1_by_i <- matrix(0, nrow = Ng, ncol = nt)
        for (i in seq_len(Ng)) {
          T0 <- rweibull(B, kshape, exp(lp0[i, jl]) / gamma(1 + 1 / kshape))
          T1 <- rweibull(B, kshape, exp(lp1[i, jl]) / gamma(1 + 1 / kshape))
          S0_by_i[i, ] <- colMeans(outer(T0, timepoints, ">"))
          S1_by_i[i, ] <- colMeans(outer(T1, timepoints, ">"))
        }
        
        # H hierarchical Dirichlet draws: global -> subgroup
        S0_acc <- rep(0, nt); S1_acc <- rep(0, nt)
        
        # LOOP: average survivals with HDP–HBB weights over H resamples
        for (h in seq_len(H)) {
          w_g <- .sample_hdp_weights_for_g(idxg, N, alpha0 = alpha0, alphag = alphag)
          # weighted average across persons
          S0_acc <- S0_acc + as.numeric(crossprod(w_g, S0_by_i))
          S1_acc <- S1_acc + as.numeric(crossprod(w_g, S1_by_i))
        }
        surv0[j, ] <- surv0[j, ] + (S0_acc / H)
        surv1[j, ] <- surv1[j, ] + (S1_acc / H)
      } 
    } 
    
    ATE <- surv1 - surv0
    res[[g]] <- list(ATE = ATE, surv_0 = surv0, surv_1 = surv1,
                     timepoints = timepoints, N = Ng, H = H,
                     alpha0 = alpha0, alphag = alphag)
  }
  
  res
}




# ============== summaries/plots for tau_g =================
get_true_effects <- function(dat, K) {
  if ("true_tau" %in% names(dat)) {
    sapply(1:K, function(g) {
      vals <- unique(dat$true_tau[dat$G == g])
      vals[1]
    })
  } else {
    NULL
  }
}


# --- summarize τ_g from Stan fit ---
extract_tau_summary <- function(stan_fit, dat) {
  post <- rstan::extract(stan_fit, pars = "subgroup_cate")$subgroup_cate  # draws x K
  K <- ncol(post)
  stopifnot(K == length(unique(dat$G)))
  df <- data.frame(
    Subgroup       = 1:K,
    Posterior_Mean = colMeans(post),
    CI_Lower       = apply(post, 2, quantile, 0.025),
    CI_Upper       = apply(post, 2, quantile, 0.975),
    Subgroup_Size  = as.numeric(table(dat$G))
  )
  if ("true_tau" %in% names(dat)) {
    df$True_Effect <- sapply(1:K, function(g) unique(dat$true_tau[dat$G==g])[1])
  }
  df
}

# --- Plot: posterior mean + 95% CrI ---
plot_tau_intervals <- function(tau_sum) {
  library(ggplot2)
  p <- ggplot(tau_sum, aes(x = factor(Subgroup), y = Posterior_Mean)) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                  width = 0.15, linewidth = 0.9, color = "steelblue") +
    geom_point(size = 2.6, color = "black") +
    { if ("True_Effect" %in% names(tau_sum))
      geom_point(aes(y = True_Effect), shape = 4, size = 3.4,
                 stroke = 1.1, color = "lightblue") } +
    geom_text(aes(y = CI_Upper, label = paste0("n=", Subgroup_Size)),
              vjust = -0.8, size = 3, color = "black") +
    geom_hline(yintercept = mean(tau_sum$Posterior_Mean),
               linetype = "dotted", color = "black") +
    labs(title = "Subgroup Treatment Effects (τ_g) with 95% Credible Intervals",
         x = "Subgroup", y = "τ_g (AFT log-time scale)") +
    theme_minimal()
  p
}


tau_sum <- extract_tau_summary(fit_parametric_ic, dat)
plot_tau_intervals(tau_sum)

# weights = empirical subgroup mix from simulated data
p_g <- as.numeric(table(dat$G)) / nrow(dat)

# posterior draws of τ_g from Stan fit
post_tau <- rstan::extract(fit_parametric_ic, pars = "subgroup_cate")$subgroup_cate  # draws × K

# overall AFT effect (one value per posterior draw)
tau_bar_draws <- as.vector(post_tau %*% p_g)

# summarize
overall_AFT <- c(mean = mean(tau_bar_draws),
                 q025 = quantile(tau_bar_draws, .025),
                 q975 = quantile(tau_bar_draws, .975))
overall_AFT


