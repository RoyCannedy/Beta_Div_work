print_node_score <- function(data, node, parents = character(0),
                             beta = 0, alpha = 1, verbose = TRUE) {
  stopifnot(node %in% names(data))
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents)) stopifnot(all(parents %in% names(data)))
  
  node_values <- levels(data[[node]])
  K <- length(node_values)
  alpha_bar <- alpha * K
  
  if (verbose) {
    cat("=== DEBUG for node:", node, 
        if (length(parents)) paste0(" | parents: ", paste(parents, collapse=", ")) else " | no parents",
        "| beta =", beta, "| alpha =", alpha, "===\n")
    cat("Node has", K, "levels:", paste(node_values, collapse=", "), "\n")
    cat("alpha_bar =", alpha_bar, "\n\n")
  }
  
  total_score <- 0
  
  if (length(parents) == 0) {
    # No parents: sequential through time
    for (t in seq_len(nrow(data))) {
      obs <- data[[node]][t]
      
      # Prior counts from observations 1:(t-1)
      prior <- if (t == 1) {
        setNames(rep(0, K), node_values)
      } else {
        tab <- table(factor(data[[node]][1:(t-1)], levels = node_values))
        setNames(as.integer(tab), node_values)
      }
      
      y_j <- prior[as.character(obs)]
      y_bar <- sum(prior)
      
      if (y_bar + alpha_bar > 0) {
        p_j <- (y_j + alpha) / (y_bar + alpha_bar)
        
        if (beta == 0) {
          # KL divergence limit (log-likelihood)
          score_t <- log(p_j)
        } else {
          # Beta divergence
          # First term: (1/beta) * p_j^beta
          first_term <- (1 / beta) * (p_j^beta)
          
          # Second term: sum over all k of p_k^(beta+1) / (beta+1)
          second_term <- sum(sapply(node_values, function(val) {
            p_k <- (prior[val] + alpha) / (y_bar + alpha_bar)
            p_k^(beta + 1)
          })) / (beta + 1)
          
          score_t <- first_term - second_term
        }
        
        total_score <- total_score + score_t
        
        if (verbose && t <= 5) {  # Print first 5 observations
          cat(sprintf("t=%d: obs='%s', y_j=%d, y_bar=%d, p_j=%.4f, score_t=%.6f\n",
                      t, obs, y_j, y_bar, p_j, score_t))
        }
      }
    }
    
  } else {
    # With parents: condition on parent configuration
    parent_configs <- do.call(paste, c(data[parents], sep = "_"))
    unique_configs <- unique(parent_configs)
    
    if (verbose) {
      cat("Found", length(unique_configs), "unique parent configurations\n\n")
    }
    
    for (t in seq_len(nrow(data))) {
      obs <- data[[node]][t]
      config <- parent_configs[t]
      
      # Prior counts from observations 1:(t-1) with SAME parent config
      prior <- if (t == 1) {
        setNames(rep(0, K), node_values)
      } else {
        idx <- parent_configs[1:(t-1)] == config
        if (any(idx)) {
          tab <- table(factor(data[[node]][1:(t-1)][idx], levels = node_values))
          setNames(as.integer(tab), node_values)
        } else {
          setNames(rep(0, K), node_values)
        }
      }
      
      y_j <- prior[as.character(obs)]
      y_bar <- sum(prior)
      
      if (y_bar + alpha_bar > 0) {
        p_j <- (y_j + alpha) / (y_bar + alpha_bar)
        
        if (beta == 0) {
          score_t <- log(p_j)
        } else {
          first_term <- (1 / beta) * (p_j^beta)
          second_term <- sum(sapply(node_values, function(val) {
            p_k <- (prior[val] + alpha) / (y_bar + alpha_bar)
            p_k^(beta + 1)
          })) / (beta + 1)
          score_t <- first_term - second_term
        }
        
        total_score <- total_score + score_t
        
        if (verbose && t <= 5) {
          cat(sprintf("t=%d: obs='%s', config='%s', y_j=%d, y_bar=%d, p_j=%.4f, score_t=%.6f\n",
                      t, obs, config, y_j, y_bar, p_j, score_t))
        }
      }
    }
  }
  
  if (verbose) {
    cat("\n")
    cat(sprintf("TOTAL SCORE: %.6f\n", total_score))
    cat(sprintf("(averaged per observation: %.6f)\n", total_score / nrow(data)))
    cat("===\n\n")
  }
  
  return(invisible(total_score))
}

# --- Example usage ---
# Using the alarm_small data from your original code
library(bnlearn)
set.seed(2025)
data(alarm)
sample_size <- 200
num_vars <- 10
rows_sel <- sample(nrow(alarm), min(sample_size, nrow(alarm)))
cols_sel <- sample(names(alarm), num_vars)
alarm_small <- alarm[rows_sel, cols_sel]

my_node <- names(alarm_small)[1]
my_parents <- names(alarm_small)[2:3]

cat("Example 1: No parents, beta=0 (KL)\n")
score1 <- print_node_score(alarm_small, node = my_node, parents = character(0), 
                           beta = 0, alpha = 1)

cat("\nExample 2: No parents, beta=0.5\n")
score2 <- print_node_score(alarm_small, node = my_node, parents = character(0), 
                           beta = 0.5, alpha = 1)

cat("\nExample 3: With parents, beta=0 (KL)\n")
score3 <- print_node_score(alarm_small, node = my_node, parents = my_parents, 
                           beta = 0, alpha = 1)

cat("\nExample 4: With parents, beta=0.3\n")
score4 <- print_node_score(alarm_small, node = my_node, parents = my_parents, 
                           beta = 0.3, alpha = 1)
