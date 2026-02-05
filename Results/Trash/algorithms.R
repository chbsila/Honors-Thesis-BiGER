library(BiGER)

rank_aggregation <- function(r) {
  
  return(BiGER::BiGER(r,
                      n_r = rep(nrow(r), ncol(r)),
                      n_u = rep(0, ncol(r)))
  )
}

sim <- function(n_items,
                n_processes,
                rho=0.5){
  # ```
  # This function simulates ranked lists from the latent variable model proposed
  # by Wang et al. (2025). For now, let's focus on fully ranked lists first.
  
  # (int) num_items: Total number of items being ranked.
  # (int) num_processes: Total number of processes.
  # (vec) rho: The correlation between a gene's local importance and global importance.
  
  # Truth
  sigma_s2=rho^(-2)-1 #study variance
  mu_i <- rnorm(n_items)
  true_rank=rank(-mu_i)
  
  ## simulate data for each process
  w  <- matrix(NA,n_items,n_processes)
  r  <- matrix(NA,n_items,n_processes)
  
  for (p in 1:n_processes) {
    w[,p] <- mu_i + rnorm(n_items, mean=0, sd=sqrt(sigma_s2))
    r[,p] <- rank(-w[,p])
  }
  
  return(list(r=r,true_rank=true_rank))
}


sim_ranking <- function(n_good, n_bad, n_items) {
  # ```
  # This function simulates a batch of good ranking and a batch of
  # byzantine ranking. All the byzantine rankings collude with each other.
  
  # (int) n_good: Total number of good processes
  # (int) n_bad: Total number of bad processes.
  # (int) n_items: Total number of items being ranked.
  
  good <- sim(n_items, n_good)
  bad <- sim(n_items, n_bad) 
  return(list("good" = good, "bad" = bad))
}

# Double check the intersection model is implemented correctly 

intersection_model_B <- function(n_good, n_bad){
  
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  
  # table[i, j] = 1 means good receiver i will receive a message from j
  table <- matrix(0, nrow = n_good, ncol = n)
  
  # Each process P_i samples n - t total including itself
  
  sample_size <- (n - t - 1) # other than itself

  for (i in good_ids) {
    
    table[i, i] <- 1  # M = {i}
    
    candidates <- setdiff(1:n, i) 
    
    S_i <- sample(candidates, sample_size) # uniformly at random (Or should it be the same processes?)
    
    for (s in S_i) {
    
      table[i, s] <- 1  # M = M \union {P_s}
      good_candidates <- setdiff(good_ids, i)
      
      table[i, sample(good_ids, t)] <- 1   # M = M \union {t good processes at random}
    }
  }
  return(table)
}

sim_round <- function(r, table) {
  # ```
  # This function simulates a round of communication using intersection model B.
  # (list) r: Input rankings globally
  # (list) table: Pre-computed communication network/intersection
  
  
  results <- list()
  
  # Good process: They do legit rank aggregation on messages received
  for (i in 1:ncol(r$good$r)) { 
    r_all <- cbind(r$good$r, r$bad$r)  
    ids_used <- which(table[i, ] == 1)  # Pre-selected communication intersection B
    r_heard <- r_all[, ids_used, drop=FALSE]  
    results[[i]] <- list()
    results[[i]][["id_good"]] <- i
    results[[i]][["ids_used"]] <- ids_used      
    results[[i]][["ra"]] <- rank_aggregation(r_heard)
    results[[i]][["posterior_ranking"]] <- rank(-results[[i]]$ra$mu) 
  }
  
  # Bad process: They generate garbage 
  for (j in 1:ncol(r$bad$r)) {
    results[[ncol(r$good$r) + j]] <- list()
    results[[ncol(r$good$r) + j]][["id_good"]] <- NULL
    results[[ncol(r$good$r) + j]][["ids_used"]] <- NULL  
    results[[ncol(r$good$r) + j]][["ra"]] <- NULL
    results[[ncol(r$good$r) + j]][["posterior_ranking"]] <- sample(1:nrow(r$bad$r)) #random
  }
  
  return(list(results = results))
}

# This method simulates multiple rounds with the new algorithm
# NOT FINISHED

sim_rounds <- function(r, number_of_rounds, trust_vector) {
  
  history <- list()
  
  n_good <- ncol(r$good$r)
  n_bad  <- ncol(r$bad$r)
  n <- n_good + n_bad
  
  r_all <- cbind(r$good$r, r$bad$r)
  
  # Reservoir of each good good process i, store ids of processes whose rankings are kept
  reservoir <- vector("list", n_good)
  
  for (i in 1:n_good){
  reservoir[[i]] <- i 
  } 
  
  # Trust vector: Using a point system, a process gets a point every time its opinion is trusted by another process
  if (is.null(trust_vector)) {
    trust_vector = rep(0,n)
  }
  
  for (j in 1:number_of_rounds) {
    
    # New communication intersection/network each round
    table <- intersection_model_B(n_good, n_bad)
    results <- sim_round(r, table)
    
    # TO DO: Update reservoirs
    for (i in 1:n_good) {
      # TO DO: use sigmas to determine which ones to add or to ignore for this round
      # TO DO: update trust vector
      # TO DO : Recompute aggregation using the full reservoir
      # TO DO: Check if change in "information" is significant
    }
    
    history[[j]] <- list(results = results, trust_vector = trust_vector, table = table, reservoir = reservoir)
  }
  
  return(history)
}


