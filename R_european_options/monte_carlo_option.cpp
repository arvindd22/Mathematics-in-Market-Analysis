#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double monteCarloCall(int num_sims, double S, double K, double r, double v, double T) {
  double S_adjust = S * exp(T * (r - 0.5 * v * v));
  double payoff_sum = 0.0;
  
  for (int i = 0; i < num_sims; ++i) {
    double gauss_bm = R::rnorm(0.0, 1.0);  // Draw from N(0,1)
    double S_cur = S_adjust * exp(sqrt(v * v * T) * gauss_bm);
    payoff_sum += std::max(S_cur - K, 0.0);  // Call payoff
  }
  
  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);  // Discounted mean
}

// [[Rcpp::export]]
double monteCarloPut(int num_sims, double S, double K, double r, double v, double T) {
  double S_adjust = S * exp(T * (r - 0.5 * v * v));
  double payoff_sum = 0.0;
  
  for (int i = 0; i < num_sims; ++i) {
    double gauss_bm = R::rnorm(0.0, 1.0);  // Draw from N(0,1)
    double S_cur = S_adjust * exp(sqrt(v * v * T) * gauss_bm);
    payoff_sum += std::max(K - S_cur, 0.0);  // Put payoff
  }
  
  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);  // Discounted mean
}
