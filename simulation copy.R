## Install Packages
if (!require(BiGER)) {
  devtools::install_github("kevin931/BiGER")
}

# Source all functions
source("~/Desktop/Thesis 25 CB/algorithms copy.R")

r <- sim_ranking(n_good = 200, n_bad = 100, n_items = 50)
history <- sim_rounds(r, n_good_received = 101, number_of_rounds = 20)

#heatmap(history[[9]]$trust_matrix)

n_good <- 200 

row_colors <- ifelse(1:nrow(tm) <= n_good, "blue", "red")
col_colors <- ifelse(1:ncol(tm) <= n_good, "blue", "red")

heatmap(
  history[[19]]$trust_matrix,
  Rowv = NA,  
  Colv = NA,  
  RowSideColors = row_colors,
  ColSideColors = col_colors,
  scale = "none",
  col = colorRampPalette(c("darkred", "orange", "yellow"))(100),
  margins = c(5, 5)
)

# sigmas <- c(0.1, 0.5, 2, 10)
# alpha <- 3
# beta <- 2
# pgamma(1 / sigmas, shape = alpha, rate = beta)



