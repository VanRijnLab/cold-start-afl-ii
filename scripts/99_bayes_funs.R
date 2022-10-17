## ----------------------------------------------------------------------------
##
## Functions to run the Bayesian model with normal-gamma prior and visualise
## the model predictions.
##
## Author: Maarten van der Velde
##
## Last updated: 2019-07-15
##
## ----------------------------------------------------------------------------

require(tibble)


run_bayes_model <- function(observations) {
  # Include all observations
  alpha <- observations
  
  # Update hyperparameters
  n <- length(alpha)
  alpha_mu <- mean(alpha)
  
  mu_n <- (kappa_0 * mu_0 + n * alpha_mu) / (kappa_0 + n)
  kappa_n <- kappa_0 + n
  a_n <- a_0 + n / 2
  sum_of_squares <- sum((alpha - alpha_mu)^2)
  b_n <- b_0 + 0.5 * sum_of_squares + (kappa_0 * n * (alpha_mu - mu_0)^2) / (2 * (kappa_0 + n))
  
  
  return(tibble(mu_n = mu_n, kappa_n = kappa_n, a_n = a_n, b_n = b_n))
}

run_bayes_model_incremental <- function(observations) {
  mu_n <- c()
  kappa_n <- c()
  a_n <- c()
  b_n <- c()
  
  for(i in 1:length(observations)) {
    
    # Include i observations
    alpha <- observations[1:i]
    
    # Update hyperparameters
    n <- length(alpha)
    alpha_mu <- mean(alpha)
    
    mu_n[i] <- (kappa_0 * mu_0 + n * alpha_mu) / (kappa_0 + n)
    kappa_n[i] <- kappa_0 + n
    a_n[i] <- a_0 + n / 2
    sum_of_squares <- sum((alpha - alpha_mu)^2)
    b_n[i] <- b_0 + 0.5 * sum_of_squares + (kappa_0 * n * (alpha_mu - mu_0)^2) / (2 * (kappa_0 + n))
    
  }
  
  return(tibble(mu_n = mu_n, kappa_n = kappa_n, a_n = a_n, b_n = b_n))
}


#' Calculate the t-distribution describing the posterior predictive of a Normal-Gamma distribution.
#' 
#' @param mu_n Mean of the NG distribution.
#' @param kappa_n Number of (pseudo-)observations in the NG distribution.
#' @param a_n Shape of the NG distribution.
#' @param b_n Rate of the NG distribution.
#' @param xmin,xmax Range of x-values over which the t-distribution is calculated.
#' @param len=500 Number of samples used to calculate the t-distribution. 
#' @return A tibble of x-values and corresponding y-values.
calculate_t_distr <- function(mu_n, kappa_n, a_n, b_n, xmin = -0.25, xmax = 1, len = 500) {
  require(extraDistr)
  
  t_df <- 2 * a_n
  t_loc <- mu_n
  t_precision <- (a_n * kappa_n) / (b_n * (kappa_n + 1))
  t_scale <- sqrt(1/t_precision)
  
  x <- seq(xmin, xmax, len = len)
  y <- dlst(x, df = t_df, mu = t_loc, sigma = t_scale)
  
  return(tibble(x = x, y = y))
}

calculate_point_estimate <- function(mu_n, kappa_n, a_n, b_n) {
  d <- calculate_t_distr(mu_n, kappa_n, a_n, b_n)
  
  return(d[which.max(d$y),]$x)
}

calculate_logarithmic_pool <- function(fact_post_pred, student_post_pred) {
  lop_weighted_dist_product <- fact_post_pred$y ^ 0.5 * student_post_pred$y ^ 0.5
  lop_normalisation_const <- sum(lop_weighted_dist_product * (diff(range(fact_post_pred$x)) / length(fact_post_pred$x)))
  lop_combined_post_pred <- data.frame(x = fact_post_pred$x,
                                       y = lop_weighted_dist_product / lop_normalisation_const)
  return(lop_combined_post_pred)
}

get_mode <- function(pred_dist) {
  return(pred_dist[which.max(pred_dist$y),]$x)
}

plot_posterior_predictive <- function(mu_n, kappa_n, a_n, b_n, title = NULL) {
  tdist <- calculate_t_distr(mu_n, kappa_n, a_n, b_n)
  
  with(tdist, plot(x, y, type = "l", main = title, xlab = expression(hat(alpha)), ylab = "Density"))
  abline(v = mu_n, lty = 2, col = "red")
}

plot_posterior_predictive_stacked <- function(mu_n, kappa_n, a_n, b_n, title = NULL) {
  col <- c("blue", rev(grey.colors(length(mu_n) - 2, start = 0, end = 0.75)), "red")
  
  for(i in 1:length(mu_n)) {
    tdist <- calculate_t_distr(mu_n[i], kappa_n[i], a_n[i], b_n[i])
    
    if(i == 1) {
      with(tdist, plot(x, y, type = "l", ylim = c(0, 4), main = title, xlab = expression(hat(alpha)), ylab = "Density", col = col[i], lwd = 2))
    } else {
      with(tdist, lines(x, y, col = col[i], lwd = 2))
    }
  }
  abline(v = mu_n[length(mu_n)], lty = 2, col = col[length(mu_n)], lwd = 2)
}


plot_combined_posterior_predictive <- function(fact_post_pred, student_post_pred, combined_post_pred, fact_res, student_res, measured_alpha, title = "Model prediction") {
  with(fact_post_pred, plot(x,y, type = "l", ylim = c(0, 8), col = "red", lty = 2, main = title))
  with(student_post_pred, lines(x,y, col = "blue", lty = 3))
  with(combined_post_pred, lines(x,y, col = "purple", lwd = 2))
  lines(x = rep(tail(fact_res$mu_n, 1), 2), y = c(0, max(fact_post_pred$y)), lty = 2, col = "red")
  lines(x = rep(tail(student_res$mu_n, 1), 2), y = c(0, max(student_post_pred$y)), lty = 3, col = "blue")
  lines(x = rep(combined_post_pred[which.max(combined_post_pred$y),]$x, 2), y = c(0, max(combined_post_pred$y)), col = "purple", lwd = 2)
  points(x = measured_alpha, y = 0, pch = 16)
  legend(0.7, 7.5, legend=c("Fact model", "Student model", "Combined model", "Observed alpha"),
         col=c("red", "blue", "purple", "black"), lty = c(2,3,1,0), lwd = c(1,1,2,0), pch = c(NA,NA,NA,16), cex = 0.8)
}