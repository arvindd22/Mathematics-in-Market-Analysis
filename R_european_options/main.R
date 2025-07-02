# Load Rcpp and source the C++ Monte Carlo file
library(Rcpp)
sourceCpp("monte_carlo_option.cpp")

# S: Current stock price
# K: Strike price
# r: Risk-free interest rate
# v: Volatility (standard deviation of asset returns)
# T: Time to maturity (in years)

black_scholes_call <- function(S, K, r, v, T) {
  d1 <- (log(S / K) + (r + 0.5 * v^2) * T) / (v * sqrt(T))
  d2 <- d1 - v * sqrt(T)
  call_price <- S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  return(call_price)
}

black_scholes_put <- function(S, K, r, v, T) {
  d1 <- (log(S / K) + (r + 0.5 * v^2) * T) / (v * sqrt(T))
  d2 <- d1 - v * sqrt(T)
  put_price <- K * exp(-r * T) * pnorm(-d2) - S * pnorm(-d1)
  return(put_price)
}

################################### INPUT PARAMERTERS ###################################

# Parameters
S <- 110
K <- 100
r <- 0.05
v <- 0.2
T <- 1
num_sims <- 1e6

# Monte Carlo (Rcpp)
mc_call <- monteCarloCall(num_sims, S, K, r, v, T)
mc_put  <- monteCarloPut(num_sims, S, K, r, v, T)

# Black-Scholes (Analytical)
bs_call <- black_scholes_call(S, K, r, v, T)
bs_put  <- black_scholes_put(S, K, r, v, T)

# Print comparison
cat("Monte Carlo Call  :", mc_call, "\n")
cat("Black-Scholes Call:", bs_call, "\n")
cat("Difference         :", abs(mc_call - bs_call), "\n\n")

cat("Monte Carlo Put   :", mc_put, "\n")
cat("Black-Scholes Put :", bs_put, "\n")
cat("Difference         :", abs(mc_put - bs_put), "\n")

###################################### PLOT ###############################################

# Simulations to test
simulations <- c(1e2, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6)

# Compute BS reference price
bs_call <- black_scholes_call(S, K, r, v, T)

# Loop to collect absolute differences
abs_diff <- numeric(length(simulations))

for (i in seq_along(simulations)) {
  mc_call <- monteCarloCall(as.integer(simulations[i]), S, K, r, v, T)
  abs_diff[i] <- abs(mc_call - bs_call)
}

# Plot
plot(simulations, abs_diff,
     type = "b", col = "blue", pch = 19,
     xlab = "Number of Simulations",
     ylab = "Absolute Difference (Monte Carlo vs Black-Scholes)",
     main = "Monte Carlo Convergence to Black-Scholes",
     log = "x")  # log scale on x-axis
grid()
