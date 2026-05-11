################################
## 0. Packages and Setup
################################

if (!require(BiGER)) {
  devtools::install_github("kevin931/BiGER")
}

############################################################
# simulations.R
#
# Distributed Rank Aggregation Experiments
# ----------------------------------------------------------
#
# This script accompanies the thesis experiments for the
# BiGER-based distributed rank aggregation framework.
#
# The file is organized into three independent experimental
# sections:
#
#   PART I   : One-round performance experiments
#   PART II  : Multi-round convergence diagnostics
#   PART III : Sensitivity analysis with respect to rho
#
# Each section:
#   - runs independently,
#   - writes its own CSV output immediately,
#   - can be executed without waiting for the other sections.
#
# Parallelization is implemented using future.apply.
#
# Output Files
# ----------------------------------------------------------
#
# PART I
#   part1_intersection_results.csv
#
# PART II
#   part2_gelman_results.csv
#   part2_process_trajectories.csv
#
# PART III
#   part3_rho_results.csv
#
############################################################

############################################################
# REQUIRED LIBRARIES
############################################################

library(BiGER)
library(dplyr)
library(coda)
library(ggplot2)
library(reshape2)
library(future.apply)

############################################################
# LOAD ALGORITHMIC DEFINITIONS
############################################################

source("~/Desktop/Thesis 25 CB/algorithms.R")

############################################################
# PARALLEL BACKEND
############################################################

available_workers <- max(1, parallel::detectCores() - 1)

plan(
  multisession,
  workers = available_workers
)

############################################################
# GLOBAL PARAMETERS
############################################################

INTERSECTIONS <- c(
  "A",
  "A'",
  "B",
  "B'",
  "C",
  "C'"
)

TOTAL_PROCESSES <- 300

N_ITEMS <- 50

REPETITIONS <- 50

############################################################
# BYZANTINE MIXTURES
############################################################

# Generates:
# (299,1), (298,2), ..., (201,99)

MIXTURES <- data.frame(
  n_good = 299:201,
  n_bad  = 1:99
)

############################################################
# OUTPUT DIRECTORY
############################################################

OUTPUT_DIR <- "simulation_outputs"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

############################################################
############################################################
# PART I
#
# ONE-ROUND PERFORMANCE EXPERIMENTS
############################################################
############################################################

# Goal:
#
# For each intersection model and each Byzantine mixture,
# run one round of BiGER aggregation and compute:
#
#   - Average Kendall's tau across good processes
#   - Average Spearman correlation across good processes
#
# relative to the ground truth ranking.
#
# Each experiment is repeated 50 times.
#
# The script reports:
#
#   - Mean Kendall's tau
#   - SD Kendall's tau
#   - Mean Spearman correlation
#   - SD Spearman correlation
#
############################################################

run_part1_single <- function(intersection,
                             n_good,
                             n_bad,
                             repetition,
                             n_items = N_ITEMS) {

  set.seed(
    1000 +
      repetition +
      n_good +
      10 * n_bad
  )

  ##########################################################
  # Simulate rankings
  ##########################################################

  r <- sim_ranking(
    n_good = n_good,
    n_bad  = n_bad,
    n_items = n_items
  )

  ##########################################################
  # Construct network
  ##########################################################

  network <- switch(
    intersection,
    "A"  = intersection_model_A(n_good, n_bad),
    "A'" = intersection_model_A_prime(n_good, n_bad),
    "B"  = intersection_model_B(n_good, n_bad),
    "B'" = intersection_model_B_prime(n_good, n_bad),
    "C"  = intersection_model_C(n_good, n_bad),
    "C'" = intersection_model_C_prime(n_good, n_bad)
  )

  ##########################################################
  # Perform one round
  ##########################################################

  out <- sim_round(r, network)

  ##########################################################
  # Compute metrics
  ##########################################################

  taus <- numeric(n_good)
  spears <- numeric(n_good)

  for (i in 1:n_good) {

    ranking_i <- out$results[[i]]$posterior_ranking

    taus[i] <- kendall_tau(
      r$good$true_rank,
      ranking_i
    )

    spears[i] <- spearman_corr(
      r$good$true_rank,
      ranking_i
    )
  }

  data.frame(
    Intersection = intersection,
    Good = n_good,
    Byzantine = n_bad,
    Repetition = repetition,
    Kendall = mean(taus, na.rm = TRUE),
    Spearman = mean(spears, na.rm = TRUE)
  )
}

############################################################
# EXECUTE PART I
############################################################

cat("\n")
cat("=================================================\n")
cat("PART I : ONE-ROUND PERFORMANCE EXPERIMENTS\n")
cat("=================================================\n")

part1_jobs <- expand.grid(
  Intersection = INTERSECTIONS,
  Mix = 1:nrow(MIXTURES),
  Repetition = 1:REPETITIONS
)

part1_raw <- future_lapply(
  1:nrow(part1_jobs),
  function(idx) {

    job <- part1_jobs[idx, ]

    n_good <- MIXTURES$n_good[job$Mix]
    n_bad  <- MIXTURES$n_bad[job$Mix]

    cat(
      "PART I |",
      "Model =", job$Intersection,
      "| Good =", n_good,
      "| Byzantine =", n_bad,
      "| Rep =", job$Repetition,
      "\n"
    )

    run_part1_single(
      intersection = job$Intersection,
      n_good = n_good,
      n_bad = n_bad,
      repetition = job$Repetition
    )
  },
  future.seed = TRUE
)

part1_raw_df <- bind_rows(part1_raw)

############################################################
# Aggregate results
############################################################

part1_summary <- part1_raw_df %>%
  group_by(
    Intersection,
    Good,
    Byzantine
  ) %>%
  summarise(
    Kendall_Mean = mean(Kendall),
    Kendall_SD   = sd(Kendall),
    Spearman_Mean = mean(Spearman),
    Spearman_SD   = sd(Spearman),
    .groups = "drop"
  )

############################################################
# Save PART I results immediately
############################################################

write.csv(
  part1_summary,
  file.path(
    OUTPUT_DIR,
    "part1_intersection_results.csv"
  ),
  row.names = FALSE
)

cat("\n")
cat("PART I COMPLETE\n")
cat("Results written to:\n")
cat("part1_intersection_results.csv\n")
cat("\n")

############################################################
############################################################
# PART II
#
# MULTI-ROUND CONVERGENCE DIAGNOSTICS
############################################################
############################################################

# Goal:
#
# Under the configuration:
#
#   (Good, Byzantine) = (201, 99)
#
# perform:
#
#   - 10 rounds of multi-round BiGER
#   - 50 repetitions
#
# For each repetition:
#
#   - Treat each good process trajectory as an MCMC chain
#   - Compute Gelman-Rubin PSRF statistics
#
# Report:
#
#   - Mean PSRF
#   - SD PSRF
#   - Maximum PSRF
#   - SD of maximum PSRF
#
# Additional Output:
#
#   part2_process_trajectories.csv
#
# which stores the Kendall's tau progression of every
# good process across all rounds and repetitions.
#
############################################################

run_part2_single <- function(intersection,
                             repetition) {

  set.seed(
    5000 +
      repetition
  )

  ##########################################################
  # Simulate rankings
  ##########################################################

  r <- sim_ranking(
    n_good = 201,
    n_bad  = 99,
    n_items = N_ITEMS
  )

  ##########################################################
  # Run multi-round BiGER
  ##########################################################

  out <- sim_rounds(
    r = r,
    number_of_rounds = 10,
    intersection = intersection,
    k = 10
  )

  ##########################################################
  # Extract Kendall trajectories
  ##########################################################

  history_tau <- out$history_tau

  tau_matrix <- matrix(
    NA,
    nrow = 10,
    ncol = 201
  )

  for (round_idx in 1:10) {

    tau_matrix[round_idx, ] <- history_tau[[round_idx]]
  }

  ##########################################################
  # Process-level trajectories
  ##########################################################

  trajectory_df <- data.frame()

  for (process_idx in 1:201) {

    process_df <- data.frame(
      Intersection = intersection,
      Repetition = repetition,
      Process = process_idx,
      Round = 1:10,
      Kendall = tau_matrix[, process_idx]
    )

    trajectory_df <- rbind(
      trajectory_df,
      process_df
    )
  }

  ##########################################################
  # Build MCMC chains
  ##########################################################

  chains <- lapply(
    1:201,
    function(i) {
      mcmc(tau_matrix[, i])
    }
  )

  chains_list <- mcmc.list(chains)

  ##########################################################
  # Gelman-Rubin statistic
  ##########################################################

  psrf <- gelman.diag(
    chains_list,
    autoburnin = FALSE
  )$psrf

  summary_df <- data.frame(
    Intersection = intersection,
    Repetition = repetition,
    PSRF_Mean = mean(psrf[, 1]),
    PSRF_Max  = max(psrf[, 1])
  )

  return(list(
    summary = summary_df,
    trajectories = trajectory_df
  ))
}

############################################################
# EXECUTE PART II
############################################################

cat("\n")
cat("=================================================\n")
cat("PART II : GELMAN-RUBIN DIAGNOSTICS\n")
cat("=================================================\n")

part2_jobs <- expand.grid(
  Intersection = INTERSECTIONS,
  Repetition = 1:REPETITIONS
)

part2_raw <- future_lapply(
  1:nrow(part2_jobs),
  function(idx) {

    job <- part2_jobs[idx, ]

    cat(
      "PART II |",
      "Model =", job$Intersection,
      "| Rep =", job$Repetition,
      "\n"
    )

    run_part2_single(
      intersection = job$Intersection,
      repetition = job$Repetition
    )
  },
  future.seed = TRUE
)

############################################################
# Separate outputs
############################################################

part2_summary_raw <- bind_rows(
  lapply(part2_raw, function(x) x$summary)
)

part2_trajectory_raw <- bind_rows(
  lapply(part2_raw, function(x) x$trajectories)
)

############################################################
# Aggregate Gelman-Rubin statistics
############################################################

part2_summary <- part2_summary_raw %>%
  group_by(Intersection) %>%
  summarise(
    Mean_PSRF = mean(PSRF_Mean),
    SD_PSRF   = sd(PSRF_Mean),
    Max_PSRF  = mean(PSRF_Max),
    SD_Max_PSRF = sd(PSRF_Max),
    .groups = "drop"
  )

############################################################
# Save PART II summary statistics
############################################################

write.csv(
  part2_summary,
  file.path(
    OUTPUT_DIR,
    "part2_gelman_results.csv"
  ),
  row.names = FALSE
)

############################################################
# Save process trajectories
############################################################

write.csv(
  part2_trajectory_raw,
  file.path(
    OUTPUT_DIR,
    "part2_process_trajectories.csv"
  ),
  row.names = FALSE
)

cat("\n")
cat("PART II COMPLETE\n")
cat("Results written to:\n")
cat("  - part2_gelman_results.csv\n")
cat("  - part2_process_trajectories.csv\n")
cat("\n")

############################################################
############################################################
# PART III
#
# RHO SENSITIVITY ANALYSIS
############################################################
############################################################

# Goal:
#
# Investigate the effect of varying rho on:
#
#   - Initial ranking quality
#   - Post-BiGER ranking quality
#
# under:
#
#   - One round of aggregation
#   - Intersection model A'
#   - 50 repetitions
#
# The experiment compares:
#
#   - Kendall's tau before aggregation
#   - Kendall's tau after aggregation
#
############################################################

RHO_VALUES <- c(
  0.1,
  0.3,
  0.5,
  0.8
)

run_part3_single <- function(rho,
                             repetition) {

  set.seed(
    9000 +
      repetition +
      round(100 * rho)
  )

  ##########################################################
  # Simulate good rankings
  ##########################################################

  good <- sim(
    n_items = N_ITEMS,
    n_processes = 201,
    rho = rho
  )

  ##########################################################
  # Initial Kendall's tau
  ##########################################################

  initial_tau <- numeric(201)

  for (i in 1:201) {

    initial_tau[i] <- kendall_tau(
      good$true_rank,
      good$r[, i]
    )
  }

  initial_mean_tau <- mean(initial_tau)

  ##########################################################
  # Simulate Byzantine rankings
  ##########################################################

  bad_rankings <- replicate(
    99,
    sample(1:N_ITEMS)
  )

  r_full <- list(
    good = good,
    bad = list(r = bad_rankings)
  )

  ##########################################################
  # One round under A'
  ##########################################################

  network <- intersection_model_A_prime(
    n_good = 201,
    n_bad  = 99
  )

  out <- sim_round(
    r_full,
    network
  )

  ##########################################################
  # Post-round Kendall's tau
  ##########################################################

  final_tau <- numeric(201)

  for (i in 1:201) {

    final_tau[i] <- kendall_tau(
      good$true_rank,
      out$results[[i]]$posterior_ranking
    )
  }

  final_mean_tau <- mean(final_tau)

  data.frame(
    Rho = rho,
    Repetition = repetition,
    Initial_Kendall = initial_mean_tau,
    Final_Kendall = final_mean_tau
  )
}

############################################################
# EXECUTE PART III
############################################################

cat("\n")
cat("=================================================\n")
cat("PART III : RHO SENSITIVITY ANALYSIS\n")
cat("=================================================\n")

part3_jobs <- expand.grid(
  Rho = RHO_VALUES,
  Repetition = 1:REPETITIONS
)

part3_raw <- future_lapply(
  1:nrow(part3_jobs),
  function(idx) {

    job <- part3_jobs[idx, ]

    cat(
      "PART III |",
      "rho =", job$Rho,
      "| Rep =", job$Repetition,
      "\n"
    )

    run_part3_single(
      rho = job$Rho,
      repetition = job$Repetition
    )
  },
  future.seed = TRUE
)

part3_raw_df <- bind_rows(part3_raw)

############################################################
# Aggregate statistics
############################################################

part3_summary <- part3_raw_df %>%
  group_by(Rho) %>%
  summarise(
    Initial_Kendall_Mean = mean(Initial_Kendall),
    Initial_Kendall_SD   = sd(Initial_Kendall),
    Final_Kendall_Mean   = mean(Final_Kendall),
    Final_Kendall_SD     = sd(Final_Kendall),
    .groups = "drop"
  )

############################################################
# Save PART III results immediately
############################################################

write.csv(
  part3_summary,
  file.path(
    OUTPUT_DIR,
    "part3_rho_results.csv"
  ),
  row.names = FALSE
)

cat("\n")
cat("PART III COMPLETE\n")
cat("Results written to:\n")
cat("part3_rho_results.csv\n")
cat("\n")
