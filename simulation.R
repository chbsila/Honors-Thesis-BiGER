## Install Packages
if (!require(BiGER)) {
  devtools::install_github("kevin931/BiGER")
}

# Source all functions
source("~/Desktop/Thesis 25 CB/algorithms.R")

r <- sim_ranking(n_good = 200, n_bad = 100, n_items = 50)

# history <- sim_rounds(r, n_good_received = 101, number_of_rounds = 10)
# history[[1]]$trust_matrix   
# history[[7]]$trust_matrix  
# heatmap(history[[5]]$trust_matrix)

table   <- intersection_model_B(200, 100)
round  <- sim_round(r, table)$results
taus <- sapply(1:200, function(i) cor(round[[i]]$posterior_ranking, r$good$true_rank, method="kendall"))
plot(taus, pch=19, col="blue", main="Kendall's Tau", xlab="Good Process ID", ylab="Tau", ylim=c(0,1))
abline(h=mean(taus), col="red", lty=2)

