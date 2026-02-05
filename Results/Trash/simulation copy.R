## Install Packages
if (!require(BiGER)) {
  devtools::install_github("kevin931/BiGER")
}

# Source all functions
source("~/Desktop/Thesis 25 CB/algorithms copy.R")

n_bad <- 100
n_items <- 50
n_good_values <- c(101:109, seq(110, 200, by = 10))

results_df <- data.frame(
  n_good = integer(),
  n_bad = integer(),
  spearman = numeric(),
  kendall = numeric()
)

# Simulation where I increase the number of good processes 
for (n_good in n_good_values) {
  
  # simulate
  r <- sim_ranking(n_good = n_good, n_bad = n_bad, n_items = n_items)
  
  history <- sim_rounds(r, n_good_received = 101, number_of_rounds = 10)
  
  # Aggregated ranking
  #good_mu <- history[[10]]$results[[1]]$ra 
  final_results <- history[[10]]$results[[1]]$ra #number_of_rounds = 10
  aggregated_rank <- rank(-final_results$mu)
  
  # Good Process Truth 
  true_rank <- r$good$true_rank
  
  spearman <- cor(true_rank, aggregated_rank, method = "spearman")
  kendall  <- cor(true_rank, aggregated_rank, method = "kendall")
  
  results_df <- rbind(results_df, data.frame(
    n_good = n_good,
    n_bad = n_bad,
    spearman = spearman,
    kendall = kendall
  ))
}

write.csv(results_df, "~/Desktop/Thesis 25 CB/simulation_results.csv", row.names = FALSE)

#r <- sim_ranking(n_good = 200, n_bad = 100, n_items = 50)
#history <- sim_rounds(r, n_good_received = 101, number_of_rounds = 20)

#heatmap(history[[9]]$trust_matrix)

#n_good <- 200 

#row_colors <- ifelse(1:nrow(tm) <= n_good, "blue", "red")
#col_colors <- ifelse(1:ncol(tm) <= n_good, "blue", "red")

#heatmap(
  #history[[19]]$trust_matrix,
  #Rowv = NA,  
  #Colv = NA,  
  #RowSideColors = row_colors,
  #ColSideColors = col_colors,
  #scale = "none",
  #col = colorRampPalette(c("darkred", "orange", "yellow"))(100),
  #margins = c(5, 5)
#)

# sigmas <- c(0.1, 0.5, 2, 10)
# alpha <- 3
# beta <- 2
# pgamma(1 / sigmas, shape = alpha, rate = beta)




