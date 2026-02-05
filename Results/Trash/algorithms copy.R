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

apply_cuts <- function(results, trust_matrix, r, max_cuts_per_process, punishment = 0.5, recovery_rate = 0.05) {
  for (i in 1:ncol(r$good$r)) {
    ids <- results[[i]]$ids_used # record actual ids used in RA
    if (is.null(ids)) next # skip 
    sigmas <- results[[i]]$ra$sigma2
    if (is.null(sigmas)) next # skip 
    
    names(sigmas) <- ids
    
    # Naive Approach: top_ids <- as.numeric(names(sort(sigmas, decreasing = TRUE))[1:max_cuts_per_process])
    sigmas <- (sigmas-min(sigmas))/(max(sigmas)-min(sigmas))
    
    # Inverse Gamma Approach
    
    #alpha <- 5

    #beta <- 5
    
    #trust_probabilities <- pgamma(1 / sigmas, shape = alpha, rate = 1 / beta)
    #names(trust_probabilities) <- ids
  
    # Stochastic and Dynamic with second chances
    #for (k in 1:ncol(trust_matrix)) {
      #if (k %in% names(trust_probabilities)) {
        # Punishment approach: trust_matrix[i, k] <- trust_matrix[i, k] * punishment
        #trust_matrix[i,k] <- trust_probabilities[as.character(k)] # Needs to be a little more sophisticated 
      #} 
    for (k in 1:ncol(trust_matrix)) {
      if (k %in% names(sigmas)) {
      trust_matrix[i,k] <- 1 - sigmas[as.character(k)] }
      else {
        trust_matrix[i, k] <- min(1, trust_matrix[i, k] + recovery_rate) 
      } 
    }
    
  }
  
  return(trust_matrix)
}

# TO DO: different intersection model + building trust up 

sim_round <- function(r, n_good_received = 101, trust_matrix) {
  # ```
  # This function simulates a round of communication. All good processes
  # get all byzantine inputs, and they do the hard work of rank aggregation.
  # They report back results.
  
  # (list) r: Input rankings globally
  # (int) n_good_received: Total number of received signals from good processes.
  
  # This implementation degrades trust with some forgetting/forgiving 
  
  # First Round
  if (missing(trust_matrix)) {
    trust_matrix <- matrix(1, ncol(r$good$r) + ncol(r$bad$r),
                           ncol(r$good$r) + ncol(r$bad$r))
  }
  
  results <- list()
  
  # Good process: They do legit rank aggregation on trusted peers
  for (i in 1:ncol(r$good$r)) {  # TO DO: Randomize order
    
    # Which processes does process i (other than itself) trust?    
    
    trusted_peers <- rbinom(length(trust_matrix[i, ]), 1, trust_matrix[i, ]) # Sampling a vector of independent Bernoulli outcomes
    
    trusted_good <- setdiff(which(trusted_peers[1:ncol(r$good$r)] == 1), i)
    
    trusted_bad  <- which(trusted_peers[(ncol(r$good$r)+1):(ncol(r$good$r)+ncol(r$bad$r))] == 1) + ncol(r$good$r)
    
    if (length(trusted_good) >= n_good_received) {
      id_good <- sample(trusted_good, n_good_received)
    } 
    else {
      id_good <- trusted_good
    }
    
    ids_used <- unique(c(id_good, trusted_bad, i)) 
    
    r_all <- cbind(r$good$r, r$bad$r)  
    r_heard <- r_all[, ids_used, drop=FALSE]  
    results[[i]] <- list()
    results[[i]][["id_good"]] <- id_good
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
    results[[ncol(r$good$r) + j]][["posterior_ranking"]] <- sample(1:nrow(r$bad$r))
  }
  
  return(list(results = results, trust_matrix = trust_matrix))
}

# This method simulates multiple rounds 
# TO DO: A notion of convergence? Kendall's Tau? Enough processes have colluded?

sim_rounds <- function(r, n_good_received = 101, number_of_rounds = 5,
                       max_cuts_per_process = 3, trust_matrix = NULL,
                       punishment = 0.5, recovery_rate = 0.05) {
  
  history <- list()
  n_all <- ncol(r$good$r) + ncol(r$bad$r)
  
  if (is.null(trust_matrix)) {
    trust_matrix <- matrix(1, ncol(r$good$r), n_all)
  }
  
  for (j in 1:number_of_rounds) { # Note for further improvement: run BiGER in parallel to make this scale faster. Has Dr. Wang done this?
    
    round_out <- sim_round(r, n_good_received = n_good_received, trust_matrix = trust_matrix)
    
    results <- round_out$results
    
    trust_matrix <- round_out$trust_matrix
    
    # Apply cutting based on current round results
    trust_matrix <- apply_cuts(results, trust_matrix, r, max_cuts_per_process, punishment, recovery_rate)
    history[[j]] <- list(results = results, trust_matrix = trust_matrix)
  }
  
  return(history)
}
