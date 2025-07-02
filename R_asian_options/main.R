library(Rcpp)
library(ggplot2)
library(microbenchmark)
library(dplyr)
library(gridExtra)
sourceCpp("asian_option.cpp")


# ---------------- R Implementations ----------------

price_asian_arithmetic_mc_r <- function(S0, K, r, sigma, T, n_steps, n_sims) {
  dt <- T / n_steps
  discount <- exp(-r * T)
  payoff_sum <- 0
  for (i in 1:n_sims) {
    St <- S0
    path <- numeric(n_steps)
    for (j in 1:n_steps) {
      Z <- rnorm(1)
      St <- St * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      path[j] <- St
    }
    avg_price <- mean(path)
    payoff_sum <- payoff_sum + max(avg_price - K, 0)
  }
  return(discount * payoff_sum / n_sims)
}

price_asian_geometric_mc_r <- function(S0, K, r, sigma, T, n_steps, n_sims) {
  dt <- T / n_steps
  discount <- exp(-r * T)
  payoff_sum <- 0
  for (i in 1:n_sims) {
    St <- S0
    log_sum <- 0
    for (j in 1:n_steps) {
      Z <- rnorm(1)
      St <- St * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
      log_sum <- log_sum + log(St)
    }
    geo_price <- exp(log_sum / n_steps)
    payoff_sum <- payoff_sum + max(geo_price - K, 0)
  }
  return(discount * payoff_sum / n_sims)
}

# ---------------- Benchmark Parameters ----------------

S0 <- 100; K <- 100; r <- 0.05; sigma <- 0.2; T <- 1; n_steps <- 50
sim_counts <- c(1000, 2000, 5000, 10000, 20000)
results <- data.frame()

# ---------------- Benchmark Loop ----------------

for (n_sims in sim_counts) {
  time_cpp_arith <- microbenchmark(
    price_asian_arithmetic_mc(S0, K, r, sigma, T, n_steps, n_sims),
    times = 3
  )
  time_r_arith <- microbenchmark(
    price_asian_arithmetic_mc_r(S0, K, r, sigma, T, n_steps, n_sims),
    times = 3
  )
  
  time_cpp_geo <- microbenchmark(
    price_asian_geometric_mc(S0, K, r, sigma, T, n_steps, n_sims),
    times = 3
  )
  time_r_geo <- microbenchmark(
    price_asian_geometric_mc_r(S0, K, r, sigma, T, n_steps, n_sims),
    times = 3
  )
  
  cpp_arith_time <- median(time_cpp_arith$time) / 1e6
  r_arith_time <- median(time_r_arith$time) / 1e6
  cpp_geo_time <- median(time_cpp_geo$time) / 1e6
  r_geo_time <- median(time_r_geo$time) / 1e6
  
  results <- rbind(results, data.frame(
    Type = "Arithmetic", Simulations = n_sims, Method = "C++", Time = cpp_arith_time, Speedup = NA
  ))
  results <- rbind(results, data.frame(
    Type = "Arithmetic", Simulations = n_sims, Method = "R", Time = r_arith_time,
    Speedup = round((r_arith_time - cpp_arith_time) / r_arith_time * 100, 1)
  ))
  results <- rbind(results, data.frame(
    Type = "Geometric", Simulations = n_sims, Method = "C++", Time = cpp_geo_time, Speedup = NA
  ))
  results <- rbind(results, data.frame(
    Type = "Geometric", Simulations = n_sims, Method = "R", Time = r_geo_time,
    Speedup = round((r_geo_time - cpp_geo_time) / r_geo_time * 100, 1)
  ))
}

# ---------------- Plot Setup ----------------

create_speedup_legend <- function(data, title) {
  legend_data <- filter(data, Method == "R") %>%
    mutate(Label = paste(Simulations, ": ", Speedup, "% faster", sep = ""))
  legend_text <- paste0("C++ over R\n", paste(legend_data$Label, collapse = "\n"))
  return(legend_text)
}

# Arithmetic Plot
arith_data <- filter(results, Type == "Arithmetic")
legend_arith <- create_speedup_legend(arith_data, "Arithmetic")

p1 <- ggplot(arith_data, aes(x = Simulations, y = Time, color = Method)) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "Arithmetic Asian Option", y = "Execution Time (ms)") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = legend_arith, hjust = 1, vjust = 1.5, size = 3)
plot(p1)

# Geometric Plot
geo_data <- filter(results, Type == "Geometric")
legend_geo <- create_speedup_legend(geo_data, "Geometric")

p2 <- ggplot(geo_data, aes(x = Simulations, y = Time, color = Method)) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  labs(title = "Geometric Asian Option", y = "Execution Time (ms)") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = legend_geo, hjust = 1, vjust = 1.5, size = 3)

# Combine
plot(p2)

