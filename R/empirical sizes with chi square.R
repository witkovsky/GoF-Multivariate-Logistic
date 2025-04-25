generate_mlogistic <- function(n, d ) {
  Y <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = diag(d))
  U <- KSinv(runif(n))  # ensure this returns KS-distributed variates
  sqrt_factor <- 1 / sqrt(U)
  X <- t(sqrt_factor * t(Y))
  return(X)
}
# Empirical characteristic function
emp_cf <- function(Z, t) {
  apply(t, 1, function(u) {
    tu <- as.vector(Z %*% u)
    mean(exp(1i * tu))
  })
}
cf_logistic <- function(t) {
  norm_t <- sqrt(rowSums(t^2))
  res <- sqrt(3) * norm_t / sinh(sqrt(3) * norm_t)
  res[is.nan(res)] <- 1
  return(res)
}

compute_Rna <- function(Z, a = 3, T = 5, grid_pts = 10) {
  d <- ncol(Z)
  grid <- expand.grid(replicate(d, seq(-T, T, length.out = grid_pts), simplify = FALSE))
  t_mat <- as.matrix(grid)
  phi_n <- emp_cf(Z, t_mat)
  phi_0 <- cf_logistic(t_mat)
  diff_sq <- (Re(phi_n) - phi_0)^2 + (Im(phi_n))^2
  weights <- exp(-a * sqrt(rowSums(t_mat^2)))
  volume_element <- (2 * T / grid_pts)^d
  Rna <- nrow(Z) * sum(diff_sq * weights) * volume_element
  return(Rna)
}


simulate_size_comparison <- function(n, d, a = 3.5, alpha = 0.05, n_sim = 10000, bins = 5) {
  # 1. Compute critical value for Rna under null
  Rvals <- numeric(n_sim)
  for (i in 1:n_sim) {
    X <- generate_mlogistic(n, d)
    mu_hat <- colMeans(X)
    Sigma_hat <- cov(X)
    Z <- t(apply(X, 1, function(x) solve(chol(Sigma_hat)) %*% (x - mu_hat)))
    Rvals[i] <- compute_Rna(Z, a)
  }
  crit_val <- quantile(Rvals, 1 - alpha)
  chi2_crit <- qchisq(1 - alpha, df = bins - 1)
 
  # 2. Run size simulations
  reject_rna <- reject_ks <- reject_energy <- reject_chi2 <- 0
 
  for (i in 1:n_sim) {
    X <- generate_mlogistic(n, d)
   
    # --- Rn,a ---
    mu_hat <- colMeans(X)
    Sigma_hat <- cov(X)
    Z <- t(apply(X, 1, function(x) solve(chol(Sigma_hat)) %*% (x - mu_hat)))
    R_stat <- compute_Rna(Z, a)
    if (R_stat > crit_val) reject_rna <- reject_rna + 1
   
    # --- KS ---
    S2 <- generate_mlogistic(n, d)
    ks_p <- fasano.franceschini.test(X, S2)$p.value
    if (ks_p < alpha) reject_ks <- reject_ks + 1
   
    # --- Energy Test ---
    Y <- generate_mlogistic(n, d)
    e_p <- eqdist.etest(rbind(X, Y), sizes = c(n, n), R = 500)$p.value
    if (e_p < alpha) reject_energy <- reject_energy + 1
   
    # --- Chi-squared Test ---
    # Mahalanobis distances
    dists <- mahalanobis(X, center = mu_hat, cov = Sigma_hat)
    edges <- qchisq(seq(0, 1, length.out = bins + 1), df = d)
    O <- hist(dists, breaks = edges, plot = FALSE)$counts
    probs <- diff(pchisq(edges, df = d))
    E <- probs * n
   
    # Check validity
    if (all(E >= 5)) {
      chi2_stat <- sum((O - E)^2 / E)
      if (chi2_stat > chi2_crit) reject_chi2 <- reject_chi2 + 1
    }
  }
 
  return(data.frame(
    SampleSize = n,
    Dimension = d,
    CriticalValue_Rna = round(crit_val, 4),
    EmpiricalSize_Rna = round(reject_rna / n_sim, 4),
    EmpiricalSize_KS = round(reject_ks / n_sim, 4),
    EmpiricalSize_Energy = round(reject_energy / n_sim, 4),
    EmpiricalSize_Chi2 = round(reject_chi2 / n_sim, 4)
  ))
}
results <- data.frame()
settings <- expand.grid(n = c(50, 100), d = c(2, 3))

for (i in 1:nrow(settings)) {
  n <- settings$n[i]
  d <- settings$d[i]
  cat(sprintf("Simulating size comparison for n = %d, d = %d...\n", n, d))
  res <- simulate_size_comparison(n = n, d = d)
  results <- rbind(results, res)
}


results