library(BiGER)
library(ggplot2)
library(reshape2)
library(dplyr)

# Load your algorithm definitions
source("~/Desktop/Thesis 25 CB/algorithms.R")

run_experiment <- function(intersections = c("A", "A'", "B", "C", "C'"), 
                           k_values = 2:10,
                           n_good = 201,
                           n_bad  = 100,
                           n_items = 50,
                           max_rounds = 50,
                           repetitions = 20,
                           metric = "kendall",
                           seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  results_df <- data.frame()
  
  for (intersection in intersections) {
    for (k in k_values) {
      for (rep in 1:repetitions) {
        
        cat("Running Model", intersection,
            "| k =", k,
            "| Rep =", rep, "\n")
        
        # Generate fresh ranking data
        r <- sim_ranking(n_good, n_bad, n_items)
        
        # Run multi-round simulation
        out <- sim_rounds(
          r = r,
          number_of_rounds = max_rounds,
          intersection = intersection,
          k = k
        )
        
        # Choose metric
        metric_values <- if (metric == "kendall") {
          out$tau
        } else {
          out$spearman
        }
        
        df_temp <- data.frame(
          Round = 1:max_rounds,
          Value = metric_values,
          Intersection = intersection,
          k = k,
          Rep = rep
        )
        
        results_df <- rbind(results_df, df_temp)
      }
    }
  }
  
  return(results_df)
}

############################################################
# AGGREGATE RESULTS (Mean + Confidence Intervals)
############################################################

aggregate_results <- function(results_df) {
  
  summary_df <- results_df %>%
    group_by(Round, Intersection, k) %>%
    summarise(
      Mean = mean(Value),
      SD   = sd(Value),
      N    = n(),
      SE   = SD / sqrt(N),
      Lower = Mean - 1.96 * SE,
      Upper = Mean + 1.96 * SE,
      .groups = "drop"
    )
  
  return(summary_df)
}

############################################################
# PLOTTING FUNCTION
############################################################

plot_results <- function(summary_df,
                         metric = "kendall",
                         min_round = 2,
                         show_ci = FALSE,
                         facet = TRUE) {
  
  summary_df <- subset(summary_df, Round >= min_round)
  
  p <- ggplot(summary_df,
              aes(x = Round,
                  y = Mean,
                  color = as.factor(k),
                  group = interaction(k, Intersection)))
  
  p <- p + geom_line(size = 1)
  
  if (show_ci) {
    p <- p +
      geom_ribbon(aes(ymin = Lower,
                      ymax = Upper,
                      fill = as.factor(k)),
                  alpha = 0.2,
                  color = NA)
  }
  
  if (facet) {
    p <- p + facet_wrap(~ Intersection)
  }
  
  p <- p +
    labs(title = paste("Average", metric,
                       "vs Number of Rounds"),
         x = "Number of Rounds",
         y = paste("Average", metric),
         color = "k",
         fill = "k") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

############################################################
# EXECUTION BLOCK
############################################################

set.seed(1)

results <- run_experiment(
  intersections = c("A, B, C, A', C'"), # add B, C
  k_values = 10:10,
  n_good = 201,
  n_bad  = 100,
  n_items = 50,
  max_rounds = 50,
  repetitions = 30,     
  metric = "kendall",
  seed = 1
)

summary_df <- aggregate_results(results)

p <- plot_results(summary_df,
                  metric = "kendall",
                  min_round = 2,
                  show_ci = FALSE,
                  facet = TRUE)

print(p) 
