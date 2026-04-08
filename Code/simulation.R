################################
## 0. Packages and Setup
################################

if (!require(BiGER)) {
  devtools::install_github("kevin931/BiGER")
}

library(BiGER)
library(ggplot2)
library(future)
library(future.apply)

# Use multicore parallelization
plan(multicore)

# Load your algorithms
source("~/Desktop/Thesis 25 CB/algorithms.R")

################################
## 1. Single Run Evaluation
################################

evaluate_one_run <- function(n_good, n_bad, n_items, intersection) {
  
  r <- sim_ranking(n_good, n_bad, n_items)
  
  network <- switch(intersection,
                    "A"  = intersection_model_A(n_good, n_bad),
                    "A'" = intersection_model_A_prime(n_good, n_bad),
                    "B"  = intersection_model_B(n_good, n_bad),
                    "B'" = intersection_model_B_prime(n_good, n_bad),
                    "C"  = intersection_model_C(n_good, n_bad),
                    "C'" = intersection_model_C_prime(n_good, n_bad))
  
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
## 2. Parallelized Experiment Loop
################################

run_experiment_parallel <- function(n_items = 50, total_processes = 300, reps = 50) {
  
  intersections <- c("A", "A'", "B", "B'", "C", "C'")
  byzantine_counts <- 1:99
  
  # Parameter grid
  param_grid <- expand.grid(intersection = intersections,
                            n_bad = byzantine_counts)
  
  # Run in parallel
  results_list <- future_lapply(1:nrow(param_grid), function(idx) {
    int <- param_grid$intersection[idx]
    n_bad <- param_grid$n_bad[idx]
    n_good <- total_processes - n_bad
    
    run_stats <- replicate(reps, evaluate_one_run(n_good, n_bad, n_items, int), simplify = "matrix")
    run_stats <- t(run_stats)
    
    data.frame(
      intersection = int,
      n_bad = n_bad,
      mean_tau = mean(run_stats[, "mean_tau"]),
      sd_tau   = sd(run_stats[, "mean_tau"]),
      mean_spear = mean(run_stats[, "mean_spear"]),
      sd_spear   = sd(run_stats[, "mean_spear"])
    )
  }, future.seed = TRUE)
  
  results_df <- do.call(rbind, results_list)
  
  # Save CSV
  write.csv(results_df, "thing_run.csv", row.names = FALSE)
  
  return(results_df)
}

################################
## 3. Plotting Functions
################################

plot_combined_metric <- function(df, metric = "tau") {
  
  metric_col <- ifelse(metric == "tau", "mean_tau", "mean_spear")
  sd_col     <- ifelse(metric == "tau", "sd_tau", "sd_spear")
  y_lab      <- ifelse(metric == "tau", "Average Kendall's Tau", "Average Spearman Correlation")
  title_text <- ifelse(metric == "tau", "Kendall's Tau vs Byzantine Processes",
                       "Spearman Correlation vs Byzantine Processes")
  
  p <- ggplot(df, aes(x = n_bad, y = !!as.name(metric_col), color = intersection)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = !!as.name(metric_col) - !!as.name(sd_col),
                    ymax = !!as.name(metric_col) + !!as.name(sd_col),
                    fill = intersection),
                alpha = 0.2, color = NA) +
    labs(title = title_text,
         x = "Number of Byzantine Processes",
         y = y_lab) +
    theme_minimal() +
    theme(text = element_text(size = 14))
  
  print(p)
  ggsave(paste0("combined_", metric, ".png"), p, width = 8, height = 5)
}

plot_individual_intersections <- function(df, metric = "tau") {
  
  metric_col <- ifelse(metric == "tau", "mean_tau", "mean_spear")
  sd_col     <- ifelse(metric == "tau", "sd_tau", "sd_spear")
  y_lab      <- ifelse(metric == "tau", "Average Kendall's Tau", "Average Spearman Correlation")
  
  for (int in unique(df$intersection)) {
    df_int <- subset(df, intersection == int)
    
    p <- ggplot(df_int, aes(x = n_bad, y = !!as.name(metric_col))) +
      geom_line(color = "blue", size = 1) +
      geom_ribbon(aes(ymin = !!as.name(metric_col) - !!as.name(sd_col),
                      ymax = !!as.name(metric_col) + !!as.name(sd_col)),
                  alpha = 0.2, fill = "blue", color = NA) +
      labs(title = paste0(int, ": ", y_lab, " vs Byzantine Processes"),
           x = "Number of Byzantine Processes",
           y = y_lab) +
      theme_minimal() +
      theme(text = element_text(size = 14))
    
    print(p)
    ggsave(paste0("intersection_", int, "_", metric, ".png"), p, width = 8, height = 5)
  }
}

################################
## 4. Run Experiment and Plot
################################

df <- run_experiment_parallel(n_items = 50, total_processes = 300, reps = 50)

# Combined plots
plot_combined_metric(df, metric = "tau")
plot_combined_metric(df, metric = "spear")

# Individual intersection plots
plot_individual_intersections(df, metric = "tau")
plot_individual_intersections(df, metric = "spear")

