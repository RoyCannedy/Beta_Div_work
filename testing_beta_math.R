print_node_score <- function(data, node, parents = character(0),
                             beta = 0, alpha = 1) {
  stopifnot(node %in% names(data))
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents)) stopifnot(all(parents %in% names(data)))
  
  K <- length(levels(data[[node]]))
  A <- alpha
  
  cat("=== DEBUG for node:", node, 
      if (length(parents)) paste0(" | parents: ", paste(parents, collapse=", ")) else " | no parents",
      "| beta =", beta, "| alpha =", alpha, "===\n")
  
  if (beta == 0) {
    # Dirichlet-multinomial marginal log-likelihood path (your ll_dirichlet)
    ll_dirichlet <- function(counts, alpha_vec) {
      A <- sum(alpha_vec); N <- sum(counts)
      lgamma(A) - lgamma(A + N) + sum(lgamma(alpha_vec + counts) - lgamma(alpha_vec))
    }
    alpha_vec <- rep(alpha / K, K)
    
    if (!length(parents)) {
      tab <- table(data[[node]])
      y <- setNames(rep(0L, K), levels(data[[node]])); y[names(tab)] <- as.integer(tab)
      cat("Counts y:", paste(y, collapse=" "), "| N =", sum(y), "| K =", K, "\n")
      sc <- ll_dirichlet(as.integer(y), alpha_vec)
      cat("Dirichlet-MLL:", sc, "\n\n")
      return(invisible(sc))
    } else {
      f <- as.formula(paste("~", paste(c(parents, node), collapse = "+")))
      ct <- xtabs(f, data = data)
      margins <- seq_len(length(dim(ct)) - 1)
      parts <- apply(ct, margins, function(slice) {
        y <- as.integer(slice)
        ll_dirichlet(y, alpha_vec)
      })
      cat("Slices:", length(parts), "| per-slice Dirichlet-MLL head:", head(parts), "\n")
      sc <- sum(parts)
      cat("Sum Dirichlet-MLL:", sc, "\n\n")
      return(invisible(sc))
    }
  } else {
    # β-divergence-like branch (Tsallis-form expression)
    tsallis_piece <- function(y) {
      N <- sum(y)
      p_hat <- (y + A / K) / (N + A)
      t1 <- sum(p_hat^beta)/beta
      t2 <- sum(p_hat^(beta+1))/(beta+1)
      list(score = t1 - t2, N=N, p=p_hat, t1=t1, t2=t2)
    }
    
    if (!length(parents)) {
      tab <- table(data[[node]])
      y <- setNames(rep(0L, K), levels(data[[node]])); y[names(tab)] <- as.integer(tab)
      cat("Counts y:", paste(y, collapse=" "), "| N =", sum(y), "| K =", K, "\n")
      out <- tsallis_piece(as.integer(y))
      cat("p_hat:", paste(round(out$p, 6), collapse=" "),
          "| sum(p^β)/β:", round(out$t1, 6),
          "| sum(p^(β+1))/(β+1):", round(out$t2, 6), "\n")
      cat("β-score:", out$score, "\n\n")
      return(invisible(out$score))
    } else {
      f <- as.formula(paste("~", paste(c(parents, node), collapse = "+")))
      ct <- xtabs(f, data = data)
      margins <- seq_len(length(dim(ct)) - 1)
      parts <- apply(ct, margins, function(slice) {
        y <- as.integer(slice)
        tsallis_piece(y)$score
      })
      cat("Slices:", length(parts), "| per-slice β-scores head:", head(parts), "\n")
      sc <- sum(parts)
      cat("Sum β-score:", sc, "\n\n")
      return(invisible(sc))
    }
  }
}

print_node_score(chds_orig, node = "Admission",
                 parents = c("Social","Economic","Events"),
                 beta = 0, alpha = 1)

print_node_score(chds_orig, node = "Admission",
                 parents = c("Social","Economic","Events"),
                 beta = 0.3, alpha = 1)
