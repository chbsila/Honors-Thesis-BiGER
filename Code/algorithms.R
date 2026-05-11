library(BiGER)

rank_aggregation <- function(r) {
  # Wrapper for BiGER
  # r: matrix of rankings
  return(BiGER::BiGER(
    r,
    n_r = rep(nrow(r), ncol(r)),
    n_u = rep(0, ncol(r))
  ))
}

sim <- function(n_items, n_processes, rho = 0.5) {
  # Simulates fully ranked lists using latent variable model (Wang et al., 2025)
  # n_items: total number of items
  # n_processes: total number of processes
  # rho: correlation between local and global importance
  
  # Ground truth
  sigma_s2 <- rho^(-2) - 1
  mu_i <- rnorm(n_items)
  true_rank <- rank(-mu_i)
  
  # Simulate rankings for each process
  w <- matrix(NA, n_items, n_processes)
  r <- matrix(NA, n_items, n_processes)
  
  for (p in 1:n_processes) {
    w[, p] <- mu_i + rnorm(n_items, mean = 0, sd = sqrt(sigma_s2))
    r[, p] <- rank(-w[, p])
  }
  
  return(list(r = r, mu_i = mu_i, true_rank = true_rank))
}

sim_ranking <- function(n_good, n_bad, n_items) {
  # Simulate a batch of good and Byzantine rankings
  # n_good: number of good processes
  # n_bad: number of Byzantine processes
  # n_items: number of items being ranked
  
  good <- sim(n_items, n_good)
  bad <- sim(n_items, n_bad)
  return(list(good = good, bad = bad))
}

# ---------------------------
# Intersection models
# ---------------------------
intersection_model_A <- function(n_good, n_bad) {
  # Model A: each good process hears from itself and (n-t-1) others randomly
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  all_ids <- 1:n
  
  network <- matrix(0, n, n)
  
  for (i in good_ids) {
    network[i, i] <- 1
    network[i, sample(setdiff(all_ids, i), n - t - 1)] <- 1
  }
  
  return(network)
}

intersection_model_A_prime <- function(n_good, n_bad) {
  # Model A': each good process hears from itself, all Byzantine, and (n-2t-1) other good processes
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  byzantine_ids <- (n_good + 1):n
  network <- matrix(0, n, n)
  
  for (i in good_ids) {
    network[i, i] <- 1
    network[i, sample(setdiff(good_ids, i), n - 2 * t - 1)] <- 1
    network[i, byzantine_ids] <- 1
  }
  
  return(network)
}

intersection_model_B <- function(n_good, n_bad) {
  # Model B: gossip-like structure
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  all_ids <- 1:n
  
  network <- matrix(0, n, n)
  sample_size <- n - t - 1
  
  # Initial random sampling
  for (i in all_ids) {
    network[i, i] <- 1
    S_i <- sample(setdiff(all_ids, i), sample_size)
    network[i, S_i] <- 1
  }
  
  initial_network <- network
  
  # Gossip step: incorporate messages from random other processes
  for (i in good_ids) {
    S_i <- sample(setdiff(all_ids, i), sample_size)
    for (s in S_i) {
      network[i, ] <- pmax(network[i, ], initial_network[s, ])
    }
  }
  
  return(network)
}

intersection_model_B_prime <- function(n_good, n_bad) {
  # Model B': initialized using A'
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  all_ids <- 1:n
  sample_size <- n - t - 1
  
  network <- intersection_model_A_prime(n_good, n_bad)
  initial_network <- network
  
  for (i in good_ids) {
    S_i <- sample(setdiff(all_ids, i), sample_size)
    for (s in S_i) {
      network[i, ] <- pmax(network[i, ], initial_network[s, ])
    }
  }
  
  return(network)
}

intersection_model_C <- function(n_good, n_bad) {
  # Model C: nested supersets, tight intersections
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  all_ids <- 1:n
  network <- matrix(0, n_good, n)
  
  A_sets <- vector("list", t + 1)
  A_sets[[1]] <- sample(all_ids, n - t)
  remainder <- sample(setdiff(all_ids, A_sets[[1]]), t)
  
  for (k in 2:(t + 1)) {
    A_sets[[k]] <- c(A_sets[[k - 1]], remainder[k - 1])
  }
  
  membership <- sample(1:(t + 1), n_good, replace = TRUE)
  for (i in 1:n_good) {
    network[i, A_sets[[membership[i]]]] <- 1
  }
  
  return(network)
}

intersection_model_C_prime <- function(n_good, n_bad) {
  # Model C': A_1 contains all Byzantine processes
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  bad_ids <- (n_good + 1):n
  all_ids <- 1:n
  network <- matrix(0, n_good, n)
  
  A_sets <- vector("list", t + 1)
  A_sets[[1]] <- c(sample(1:n_good, n - 2 * t), bad_ids)
  remainder <- sample(setdiff(all_ids, A_sets[[1]]), t)
  
  for (k in 2:(t + 1)) {
    A_sets[[k]] <- c(A_sets[[k - 1]], remainder[k - 1])
  }
  
  membership <- sample(1:(t + 1), n_good, replace = TRUE)
  for (i in 1:n_good) {
    network[i, A_sets[[membership[i]]]] <- 1
  }
  
  return(network)
}

sim_round <- function(r, network) {
  n_good <- ncol(r$good$r)
  n_bad <- ncol(r$bad$r)
  r_all <- cbind(r$good$r, r$bad$r)
  
  results <- vector("list", n_good + n_bad)
  trust_matrix <- matrix(0, n_good, n_good + n_bad)
  
  # Good processes
  for (i in 1:n_good) {
    ids_used <- which(network[i, ] == 1)
    r_heard <- r_all[, ids_used, drop = FALSE]
    
    ra <- rank_aggregation(r_heard)
    post_now <- rank(-ra$mu)
    
    # Normalize sigmas to trust
    sigmas <- ra$sigma2
    sigmas <- (sigmas - min(sigmas)) / (max(sigmas) - min(sigmas))
    trust_matrix[i, ids_used] <- 1 - sigmas
    
    results[[i]] <- list(
      id_good = i,
      ids_used = ids_used,
      ra = ra,
      posterior_ranking = post_now
    )
  }
  
  # Bad processes: random rankings
  for (j in 1:n_bad) {
    idx <- n_good + j
    results[[idx]] <- list(posterior_ranking = sample(1:nrow(r_all)))
  }
  
  return(list(results = results, trust_matrix = trust_matrix))
}


kendall_tau <- function(true_rank, estimated_rank) {
  cor(true_rank, estimated_rank, method = "kendall")
}

spearman_corr <- function(true_rank, estimated_rank) {
  cor(true_rank, estimated_rank, method = "spearman")
}

sim_rounds <- function(r,
                       number_of_rounds,
                       intersection,
                       k) {
  
  n_good <- ncol(r$good$r)
  n_bad  <- ncol(r$bad$r)
  n_processes <- n_good + n_bad
  
  r_original <- cbind(r$good$r, r$bad$r)
  
  # Start evolving rankings
  current_rankings <- r_original
  r_all <- current_rankings
  
  col_process_id <- 1:n_processes
  
  duplication_count <- rep(0, n_processes)
  duplication_queue <- vector("list", n_processes)
  for (pid in 1:n_processes) {
    duplication_queue[[pid]] <- which(col_process_id == pid)
  }
  
  avg_tau_true <- numeric(number_of_rounds)
  avg_spear_true <- numeric(number_of_rounds)
  tau_convergence <- vector("list", number_of_rounds)
  spearman_convergence <- vector("list", number_of_rounds)
  
  for (round in 1:number_of_rounds) {
    
    if (intersection == "A") network <- intersection_model_A(n_good, n_bad)
    if (intersection == "A'") network <- intersection_model_A_prime(n_good, n_bad)
    if (intersection == "B") network <- intersection_model_B(n_good, n_bad)
    # if (intersection == "B'") network <- intersection_model_B_prime(n_good, n_bad)
    if (intersection == "C") network <- intersection_model_C(n_good, n_bad)
    if (intersection == "C'") network <- intersection_model_C_prime(n_good, n_bad)
    
    tau_true <- numeric(n_good)
    spear_true <- numeric(n_good)
    
    new_good_rankings <- matrix(NA, nrow = nrow(r_original), ncol = n_good)
    
    for (i in 1:n_good) {
      
      ids_used <- which(network[i, ] == 1)
      r_new <- current_rankings[, ids_used, drop = FALSE]
      
      # Run BiGER on everything so far + new suggestions
      r_combined <- cbind(r_all, r_new)
      
      ra <- rank_aggregation(r_combined)
      post_now <- rank(-ra$mu)
      sigmas <- ra$sigma2
      
      # Normalize trust
      sigmas_norm <- (sigmas - min(sigmas)) / (max(sigmas) - min(sigmas))
      trust_cols <- 1 - sigmas_norm
      
      # Extract trust only for newly suggested rankings
      start_new <- ncol(r_all) + 1
      end_new   <- ncol(r_combined)
      trust_new <- trust_cols[start_new:end_new]
      
      for (j in seq_along(ids_used)) {
        
        pid <- ids_used[j]
        trust_value <- trust_new[j]
        
        if (runif(1) < trust_value) {
          
          # If under cap, just add
          if (duplication_count[pid] < k) {
            
            r_all <- cbind(r_all, current_rankings[, pid, drop = FALSE])
            col_process_id <- c(col_process_id, pid)
            
            duplication_queue[[pid]] <- c(duplication_queue[[pid]], ncol(r_all))
            
            duplication_count[pid] <- duplication_count[pid] + 1
            
          } else {
            
            oldest_col <- duplication_queue[[pid]][1]
            
            # Remove column
            r_all <- r_all[, -oldest_col, drop = FALSE]
            col_process_id <- col_process_id[-oldest_col]
            
            # Adjust all stored indices after column removal
            for (p in 1:n_processes) {
              duplication_queue[[p]] <- duplication_queue[[p]][duplication_queue[[p]] != oldest_col]
              duplication_queue[[p]][duplication_queue[[p]] > oldest_col] <-
                duplication_queue[[p]][duplication_queue[[p]] > oldest_col] - 1
            }
            
            r_all <- cbind(r_all, current_rankings[, pid, drop = FALSE])
            col_process_id <- c(col_process_id, pid)
            
            new_col_index <- ncol(r_all)
            duplication_queue[[pid]] <- c(duplication_queue[[pid]][-1], new_col_index)
          }
        }
      }
      
      # Store aggregated ranking
      new_good_rankings[, i] <- post_now
      
      tau_true[i] <- kendall_tau(r$good$true_rank, post_now)
      spear_true[i] <- spearman_corr(r$good$true_rank, post_now)
    }
    
    # Update good processes with their aggregated rankings
    current_rankings[, 1:n_good] <- new_good_rankings
    
    # Bad processes remain random each round (simplest model for adversary behavior)
    for (j in 1:n_bad) {
      current_rankings[, n_good + j] <- sample(1:nrow(r_original))
    }
    
    tau_convergence[round] <- tau_true
    spearman_convergence[round] <- spear_true
    avg_tau_true[round] <- mean(tau_true, na.rm = TRUE)
    avg_spear_true[round] <- mean(spear_true, na.rm = TRUE)
  }
  
  return(list(
    tau = avg_tau_true,
    spearman = avg_spear_true,
    history_tau = tau_convergence,
    history_spearman = spearman_convergence
  ))
}
