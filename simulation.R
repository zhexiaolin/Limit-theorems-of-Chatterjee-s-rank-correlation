# Load required library
library(MASS)
library(xtable)
library(parallel)

XInvarcalculate = function(xvec, yvec){
  n = length(xvec)
  xrank = rank(xvec, ties.method = "random")
  yrank = rank(yvec, ties.method = "random")
  ord = order(xrank)
  yrank = yrank[ord]
  yrank1 = yrank[c(2:n,n)]
  yrank2 = yrank[c(3:n,n,n)]
  yrank3 = yrank[c(4:n,n,n,n)]
  term1 = pmin(yrank,yrank1)
  term2 = pmin(yrank,yrank2)
  term3 = pmin(yrank2,yrank3)
  term4 = sapply(1:n, function(i){sum(yrank[i]<=(term1[-i]))})
  term5 = pmin(yrank1,yrank2)
  sum1 = mean((term1/n)^2)
  sum2 = mean(term1*term2/n^2)
  sum3 = mean(term1*term3/n^2)
  sum4 = mean(term4*term1/(n*(n-1)))
  sum5 = mean(term4*term5/(n*(n-1)))
  sum6 = mean(sapply(1:n, function(i){sum(pmin(term1[i],term1[-i]))})/(n*(n-1)))
  sum7 = (mean(term1/n))^2
  return(max(0,36*(sum1+2*sum2-2*sum3+4*sum4-2*sum5+sum6-4*sum7)))
}

calculateXI <- function (xvec, yvec, simple = TRUE) 
{
  n <- length(xvec)
  PI <- rank(xvec, ties.method = "random")
  fr <- rank(yvec, ties.method = "max")/n
  gr <- rank((-yvec), ties.method = "max")/n
  ord <- order(PI)
  fr <- fr[ord]
  A1 <- sum(abs(fr[1:(n - 1)] - fr[2:n]))/(2 * n)
  CU <- mean(gr * (1 - gr))
  xi <- 1 - A1/CU
  if (simple == TRUE) 
    return(xi)
  else return(list(xi = xi, fr = fr, CU = CU))
}

set.seed(666)
# Set parameters
rhos <- c(0, 0.3, 0.5, 0.7, 0.9)
ns <- c(100, 1000, 10000)
n_bootstrap <- 5000
nsim <- 5000
n_true <- 10000000
n_true_var <- 100000

# Initialize a data frame to store the results
results_df <- data.frame(
  rho = numeric(),
  n = integer(),
  alpha = numeric(),
  variance_LH = numeric(),
  variance_B = numeric(),
  LH_005 = numeric(),
  B_005 = numeric(),
  LH_01 = numeric(),
  B_01 = numeric()
)

# Loop through all parameter combinations
for (rho in rhos) {
  
  # Mean and covariance matrix for bivariate normal distribution
  mean <- c(0, 0)
  cov_matrix <- matrix(c(1, rho, rho, 1), 2)
  
  # Generate samples to approximate the true value of xi
  true_samples <- mvrnorm(n_true, mu = mean, Sigma = cov_matrix)
  true_xi <- calculateXI(true_samples[, 1], true_samples[, 2])
  
  for (n in ns) {
    # Initialize variables
    coverage_count_B_005 <- 0
    coverage_count_LH_005 <- 0
    coverage_count_B_01 <- 0
    coverage_count_LH_01 <- 0
    variance_B <- numeric(nsim)
    variance_LH <- numeric(nsim)
    
    xi_n_list <- unlist(mclapply(rep(n, n_true_var), function(n) {
      samples <- mvrnorm(n, mu = mean, Sigma = cov_matrix)
      calculateXI(samples[, 1], samples[, 2])
    }, mc.cores = 20))
    
    true_xi_var <- n*var(xi_n_list)
    
    # Run simulations
    for (sim in 1:nsim) {
      # Generate n samples from the bivariate normal distribution
      samples <- mvrnorm(n, mu = mean, Sigma = cov_matrix)
      
      # Calculate the correlation coefficient xi_n
      xi_n <- calculateXI(samples[, 1], samples[, 2])
      
      # Bootstrap resampling using the replicate function and mclapply for parallelization
      xi_bootstrap <- unlist(mclapply(rep(n, n_bootstrap), function(n) {
        bootstrap_samples <- samples[sample(1:n, floor(sqrt(n)), replace = FALSE), ]
        calculateXI(bootstrap_samples[, 1], bootstrap_samples[, 2])
      }, mc.cores = 20))
      
      xi_diff <- xi_n - true_xi
      
      # Estimate the variance using the XInvarcalculate function from your local environment
      xi_variance <- XInvarcalculate(samples[, 1], samples[, 2])
      xi_variance_B <- (floor(sqrt(n)))^2*mean((xi_bootstrap - mean(xi_bootstrap))^2)
      
      # Calculate the normal distribution confidence interval using the estimated variance and mean zero
      ci_lower_LH <- -qnorm(1 - 0.05 / 2) * sqrt(xi_variance / n)
      ci_upper_LH <- qnorm(1 - 0.05 / 2) * sqrt(xi_variance / n)
      ci_lower_B <- -qnorm(1 - 0.05 / 2) * sqrt(xi_variance_B / n)
      ci_upper_B <- qnorm(1 - 0.05 / 2) * sqrt(xi_variance_B / n)

      # Check if xi_diff is within the bootstrap confidence intervals for both versions and the normal distribution interval
      if (xi_diff >= ci_lower_B && xi_diff <= ci_upper_B) {
        coverage_count_B_005 <- coverage_count_B_005 + 1
      }
      if (xi_diff >= ci_lower_LH && xi_diff <= ci_upper_LH) {
        coverage_count_LH_005 <- coverage_count_LH_005 + 1
      }
      
      # Calculate the normal distribution confidence interval using the estimated variance and mean zero
      ci_lower_LH <- -qnorm(1 - 0.1 / 2) * sqrt(xi_variance / n)
      ci_upper_LH <- qnorm(1 - 0.1 / 2) * sqrt(xi_variance / n)
      ci_lower_B <- -qnorm(1 - 0.1 / 2) * sqrt(xi_variance_B / n)
      ci_upper_B <- qnorm(1 - 0.1 / 2) * sqrt(xi_variance_B / n)
      
      # Check if xi_diff is within the bootstrap confidence intervals for both versions and the normal distribution interval
      if (xi_diff >= ci_lower_B && xi_diff <= ci_upper_B) {
        coverage_count_B_01 <- coverage_count_B_01 + 1
      }
      if (xi_diff >= ci_lower_LH && xi_diff <= ci_upper_LH) {
        coverage_count_LH_01 <- coverage_count_LH_01 + 1
      }
      
      variance_B[sim] <- xi_variance_B
      variance_LH[sim] <- xi_variance
    }
    
    # Calculate empirical coverage probabilities for all versions
    empirical_coverage_probability_B_005 <- coverage_count_B_005 / nsim
    empirical_coverage_probability_LH_005 <- coverage_count_LH_005 / nsim
    empirical_coverage_probability_B_01 <- coverage_count_B_01 / nsim
    empirical_coverage_probability_LH_01 <- coverage_count_LH_01 / nsim
    
    results_df <- rbind(results_df, data.frame(
      rho = rho,
      n = n,
      variance_LH = sqrt(mean((variance_LH-true_xi_var)^2)),
      variance_B = sqrt(mean((variance_B-true_xi_var)^2)),
      LH_005 = empirical_coverage_probability_LH_005,
      B_005 = empirical_coverage_probability_B_005,
      LH_01 = empirical_coverage_probability_LH_01,
      B_01 = empirical_coverage_probability_B_01
    ))
    
    # Print the progress message
    print(sprintf("Simulation done for rho = %.1f and n = %d", rho, n))
  }
}

# Save the results_df data frame as a CSV file
write.csv(results_df, file = "results_df.csv", row.names = FALSE)

# Set custom column names
colnames(results_df) <- c("rho", "n", "LH", "B", "LH", "B", "LH", "B")

# Convert the rearranged data frame to a LaTeX table
latex_table <- xtable(results_df, include.rownames = FALSE)
print.xtable(latex_table, type = "latex", include.rownames = FALSE)
