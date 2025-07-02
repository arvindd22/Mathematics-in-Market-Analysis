#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double price_asian_arithmetic_mc(
    double S0, double K, double r, double sigma, double T,
    int n_steps, int n_sims, bool is_call = true
) {
  double dt = T / n_steps;
  double discount = std::exp(-r * T);
  double payoff_sum = 0.0;
  
  for (int sim = 0; sim < n_sims; ++sim) {
    double St = S0;
    double sum_prices = 0.0;
    
    for (int step = 0; step < n_steps; ++step) {
      double Z = R::rnorm(0.0, 1.0);
      St *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
      sum_prices += St;
    }
    
    double avg_price = sum_prices / n_steps;
    double payoff = is_call ? std::max(avg_price - K, 0.0) : std::max(K - avg_price, 0.0);
    payoff_sum += payoff;
  }
  
  return discount * (payoff_sum / n_sims);
}

// [[Rcpp::export]]
double price_asian_geometric_mc(
    double S0, double K, double r, double sigma, double T,
    int n_steps, int n_sims, bool is_call = true
) {
  double dt = T / n_steps;
  double discount = std::exp(-r * T);
  double payoff_sum = 0.0;
  
  for (int sim = 0; sim < n_sims; ++sim) {
    double St = S0;
    double log_sum = 0.0;
    
    for (int step = 0; step < n_steps; ++step) {
      double Z = R::rnorm(0.0, 1.0);
      St *= std::exp((r - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
      log_sum += std::log(St);
    }
    
    double geo_mean = std::exp(log_sum / n_steps);
    double payoff = is_call ? std::max(geo_mean - K, 0.0) : std::max(K - geo_mean, 0.0);
    payoff_sum += payoff;
  }
  
  return discount * (payoff_sum / n_sims);
}
