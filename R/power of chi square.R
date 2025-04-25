# Required libraries
library(MASS)       # For mvrnorm
library(sn)         # For skew-normal
install.packages("sn")  # If you haven't already
library(sn)

# Generic function for power calculation
chi_squared_power <- function(n = 100, d = 2, B = 10000, bins = 5, alpha = 0.05,
                              alt_dist, seed = 123) {
  set.seed(seed)
  rejections <- 0
  critical_value <- qchisq(1 - alpha, df = bins - 1)

  for (b in 1:B) {
    X <- alt_dist(n, d)

    mu_hat <- colMeans(X)
    Sigma_hat <- cov(X)
    dists <- mahalanobis(X, center = mu_hat, cov = Sigma_hat)

    edges <- qchisq(seq(0, 1, length.out = bins + 1), df = d)
    O <- hist(dists, breaks = edges, plot = FALSE)$counts
    probs <- diff(pchisq(edges, df = d))
    E <- probs * n

    if (any(E < 5)) next  # Optional: skip invalid chi^2 cells

    chi2 <- sum((O - E)^2 / E)
    if (chi2 > critical_value) rejections <- rejections + 1
  }

  power <- rejections / B
  return(power)
}

# -----------------------
# Alternative distributions
# -----------------------

# Standard Normal
alt_normal <- function(n, d) mvrnorm(n, mu = rep(0, d), Sigma = diag(d))

# Laplace-Normal Mixture
alt_ln <- function(p) {
  function(n, d) {
    Z1 <- mvrnorm(n, rep(0, d), diag(d))
    E <- matrix(rexp(n * d), n, d)
    Z2 <- sqrt(E / 2) * matrix(rnorm(n * d), n, d)
    mask <- rbinom(n, 1, p)
    return(mask * Z2 + (1 - mask) * Z1)
  }
}

# Skew-Normal
alt_skewnorm <- function(alpha) {
  function(n, d) {
    xi <- rep(0, d)                        # Location vector
    Omega <- diag(d)                      # Covariance matrix
    alpha_vec <- rep(alpha, d)            # Skewness vector
    rmsn(n = n, xi = xi, Omega = Omega, alpha = alpha_vec)
  }
}


# t-distribution
alt_t <- function(df) {
  function(n, d) {
    X <- matrix(rt(n * d, df), n, d)
    X / sqrt(df / (df - 2))
  }
}

# Normal Mixture
alt_norm_mixture <- function(n, d) {
  X1 <- mvrnorm(n, mu = rep(0, d), Sigma = diag(d))
  X2 <- mvrnorm(n, mu = rep(3, d), Sigma = diag(d))
  return((X1 + X2) / 2)
}

# Uniform on unit cube
alt_uniform <- function(n, d) matrix(runif(n * d), n, d)

# Pearson Type II
alt_p2 <- function(n, d) {
  X <- matrix(rnorm(n * d), n, d)
  norms <- sqrt(rowSums(X^2))
  mask <- norms < 1
  while (sum(mask) < n) {
    X_new <- matrix(rnorm(n * d), n, d)
    norms_new <- sqrt(rowSums(X_new^2))
    X <- rbind(X[mask, ], X_new[norms_new < 1, ])
    mask <- rep(TRUE, nrow(X))
  }
  return(X[1:n, ])
}

# Multivariate Logistic Mixture (MLmix)
alt_ml_mix <- function(n, d) {
  M1 <- mvrnorm(n, mu = rep(0, d), Sigma = diag(d))
  M2 <- mvrnorm(n, mu = rep(3, d), Sigma = diag(d))
  return((M1 + M2) / 2)
}

# -----------------------
# Run the power study
# -----------------------

alternatives <- list(
  "Normal"         = alt_normal,
  "LN(0.25)"       = alt_ln(0.25),
  "LN(0.5)"        = alt_ln(0.5),
  "LN(0.75)"       = alt_ln(0.75),
  "SN(0.25)"       = alt_skewnorm(2),
  "SN(0.5)"        = alt_skewnorm(5),
  "SN(0.75)"       = alt_skewnorm(10),
  "t(1)"           = alt_t(1),
  "t(10)"          = alt_t(10),
  "N1"             = alt_norm_mixture,
  "U"              = alt_uniform,
  "P2"             = alt_p2,
  "MLmix"          = alt_ml_mix
)

dims <- c(2, 3)
n_vals <- c(50, 100)

results <- data.frame()

for (d in dims) {
  for (n in n_vals) {
    for (alt_name in names(alternatives)) {
      cat("Running:", alt_name, "| d =", d, "| n =", n, "\n")
      alt_func <- alternatives[[alt_name]]
      power <- chi_squared_power(n = n, d = d, B = 1000, bins = 5, alt_dist = alt_func)
      results <- rbind(results, data.frame(
        Alternative = alt_name,
        Dimension = d,
        SampleSize = n,
        Chi2Power = round(power, 4)
      ))
    }
  }
}

# View final result
print(results)
