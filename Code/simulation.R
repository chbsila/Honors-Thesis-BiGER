## Install Packages
if (!require(BiGER)) {
  devtools::install_github("kevin931/BiGER")
}

library(ggplot2)

source("~/Desktop/Thesis 25 CB/algorithms.R")

################################
# ONE SIMULATION
################################

evaluate_one_run <- function(n_good, n_bad, n_items, intersection) {
  
  r <- sim_ranking(n_good, n_bad, n_items)
  
  if (intersection == "A") network <- intersection_model_A(n_good, n_bad)
  if (intersection == "B") network <- intersection_model_B(n_good, n_bad)
  if (intersection == "C") network <- intersection_model_C(n_good, n_bad)
  
  out <- sim_round(r, network)
  results <- out$results
  
  kendall_tau <- function(a, b) suppressWarnings(cor(a, b, method="kendall"))
  spearman_corr <- function(a, b) suppressWarnings(cor(a, b, method="spearman"))
  
  tau_vals <- numeric(n_good)
  spear_vals <- numeric(n_good)
  
  for (i in 1:n_good) {
    post <- results[[i]]$posterior_ranking
    tau_vals[i] <- kendall_tau(r$good$true_rank, post)
    spear_vals[i] <- spearman_corr(r$good$true_rank, post)
  }
  
  c(mean_tau = mean(tau_vals),
    sd_tau   = sd(tau_vals),
    mean_spear = mean(spear_vals),
    sd_spear   = sd(spear_vals))
}

################################
# EXPERIMENT LOOP
################################

run_experiment <- function(n_items = 50,
                           total_processes = 301,
                           reps = 100) {
  
  intersections <- c("A","B","C")
  results_df <- data.frame()
  
  for (int in intersections) {
    cat("Running intersection", int, "\n")
    
    for (t in 1:100) {
      cat("  Byzantine:", t, "\n")
      
      n_bad  <- t
      n_good <- total_processes - t
      
      run_stats <- replicate(reps, evaluate_one_run(n_good, n_bad, n_items, int))
      run_stats <- t(run_stats)
      
      row <- data.frame(
        intersection = int,
        n_bad = n_bad,
        mean_tau = mean(run_stats[, "mean_tau"]),
        sd_tau   = sd(run_stats[, "mean_tau"]),
        mean_spear = mean(run_stats[, "mean_spear"]),
        sd_spear   = sd(run_stats[, "sd_spear"])
      )
      
      results_df <- rbind(results_df, row)
    }
  }
  
  write.csv(results_df, "intersection_experiment_full.csv", row.names = FALSE)
  return(results_df)
}

################################
# PLOT
################################

plot_results <- function(df) {
  
  p1 <- ggplot(df, aes(x = n_bad, y = mean_tau, color = intersection)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean_tau - sd_tau,
                    ymax = mean_tau + sd_tau,
                    fill = intersection),
                alpha = 0.2, color = NA) +
    labs(title = "Kendall Tau vs Byzantine Processes",
         x = "Number of Byzantine Processes",
         y = "Mean Kendall Tau") +
    theme_minimal()
  
  p2 <- ggplot(df, aes(x = n_bad, y = mean_spear, color = intersection)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean_spear - sd_spear,
                    ymax = mean_spear + sd_spear,
                    fill = intersection),
                alpha = 0.2, color = NA) +
    labs(title = "Spearman vs Byzantine Processes",
         x = "Number of Byzantine Processes",
         y = "Mean Spearman Correlation") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

################################
# RUN
################################

df <- run_experiment(n_items = 50, reps = 100)
plot_results(df)

