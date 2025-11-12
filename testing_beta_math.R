
library(bnlearn)
set.seed(42)

beta_div_score <- function(node, parents, data, beta = 0, alpha = 1.0, ...) {
  stopifnot(is.character(node), length(node) == 1, node %in% names(data))
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents) == 0) parents <- character(0)
  
  ll_dirichlet <- function(counts, alpha_vec) {
    A <- sum(alpha_vec); N <- sum(counts)
    lgamma(A) - lgamma(A + N) + sum(lgamma(alpha_vec + counts) - lgamma(alpha_vec))
  }
  
  if (isTRUE(all.equal(beta, 0))) {
    lev <- levels(data[[node]]); K <- length(lev)
    alpha_vec <- rep(alpha, K)  # per-category alpha
    
    if (length(parents) == 0) {
      tab <- table(data[[node]])
      y <- setNames(rep(0L, K), lev); y[names(tab)] <- as.integer(tab)
      return(ll_dirichlet(y, alpha_vec))
    } else {
      f <- as.formula(paste("~", paste(c(parents, node), collapse = "+")))
      ct <- xtabs(f, data = data)
      margins <- seq_len(length(dim(ct)) - 1)
      if (length(margins) == 0) {
        y <- as.integer(ct)
        return(ll_dirichlet(y, alpha_vec))
      } else {
        return(sum(apply(ct, margins, function(slice) {
          y <- as.integer(slice)
          ll_dirichlet(y, alpha_vec)
        })))
      }
    }
  }
  
  # beta > 0: per-category alpha smoothing p_hat = (y + alpha) / (N + K*alpha)
  if (length(parents) == 0) {
    lev <- levels(data[[node]]); K <- length(lev)
    tab <- table(data[[node]])
    y <- setNames(rep(0L, K), lev); y[names(tab)] <- as.integer(tab)
    N <- sum(y)
    p_hat <- (y + alpha) / (N + K * alpha)
    return(sum(p_hat^beta)/beta - sum(p_hat^(beta + 1))/(beta + 1))
  } else {
    lev <- levels(data[[node]]); K <- length(lev)
    f <- as.formula(paste("~", paste(c(parents, node), collapse = "+")))
    ct <- xtabs(f, data = data)
    margins <- seq_len(length(dim(ct)) - 1)
    return(sum(apply(ct, margins, function(slice) {
      y <- as.integer(slice); N <- sum(y)
      p_hat <- (y + alpha) / (N + K * alpha)
      sum(p_hat^beta)/beta - sum(p_hat^(beta + 1))/(beta + 1)
    })))
  }
}

# ---------- Helper: pretty print a single node score ----------
print_node_score <- function(data, node, parents = character(0), beta = 0, alpha = 1) {
  stopifnot(node %in% names(data))
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents)) stopifnot(all(parents %in% names(data)))
  
  K <- length(levels(data[[node]]))
  cat("=== DEBUG for node:", node, 
      if (length(parents)) paste0(" | parents: ", paste(parents, collapse=", ")) else " | no parents",
      "| beta =", beta, "| alpha =", alpha, "===\n")
  
  if (beta == 0) {
    ll_dirichlet <- function(counts, alpha_vec) {
      A <- sum(alpha_vec); N <- sum(counts)
      lgamma(A) - lgamma(A + N) + sum(lgamma(alpha_vec + counts) - lgamma(alpha_vec))
    }
    alpha_vec <- rep(alpha, K)
    
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
    tsallis_piece <- function(y, alpha, K, beta) {
      N <- sum(y)
      p_hat <- (y + alpha) / (N + K * alpha)
      t1 <- sum(p_hat^beta)/beta
      t2 <- sum(p_hat^(beta+1))/(beta+1)
      list(score = t1 - t2, N=N, p=p_hat, t1=t1, t2=t2)
    }
    
    if (!length(parents)) {
      tab <- table(data[[node]])
      y <- setNames(rep(0L, K), levels(data[[node]])); y[names(tab)] <- as.integer(tab)
      cat("Counts y:", paste(y, collapse=" "), "| N =", sum(y), "| K =", K, "\n")
      out <- tsallis_piece(as.integer(y), alpha, K, beta)
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
        tsallis_piece(y, alpha, K, beta)$score
      })
      cat("Slices:", length(parts), "| per-slice β-scores head:", head(parts), "\n")
      sc <- sum(parts)
      cat("Sum β-score:", sc, "\n\n")
      return(invisible(sc))
    }
  }
}

# ---------- KL vs beta=0 equivalence check ----------
kl_score_dirichlet <- function(y, alpha, K) {
  alpha_vec <- rep(alpha, K)
  lgamma(sum(alpha_vec)) - lgamma(sum(alpha_vec) + sum(y)) +
    sum(lgamma(alpha_vec + y) - lgamma(alpha_vec))
}

check_beta0_vs_kl <- function(data, node, parents = character(0), alpha = 1) {
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  lev <- levels(data[[node]]); K <- length(lev)
  
  if (length(parents) == 0) {
    tab <- table(data[[node]])
    y <- setNames(rep(0L, K), lev); y[names(tab)] <- as.integer(tab)
    kl <- kl_score_dirichlet(as.integer(y), alpha, K)
  } else {
    f <- as.formula(paste("~", paste(c(parents, node), collapse = "+")))
    ct <- xtabs(f, data = data)
    margins <- seq_len(length(dim(ct)) - 1)
    kl <- if (length(margins) == 0) {
      kl_score_dirichlet(as.integer(ct), alpha, K)
    } else {
      sum(apply(ct, margins, function(slice) kl_score_dirichlet(as.integer(slice), alpha, K)))
    }
  }
  
  b0 <- beta_div_score(node, parents, data, beta = 0, alpha = alpha)
  cat(sprintf("KL: %.6f | beta=0: %.6f | |diff|=%.6g\n", kl, b0, abs(kl - b0)))
  invisible(c(KL = kl, beta0 = b0))
}

# ---------- Wrapper factory ----------
make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}

# ---------- Runner across betas for a dataset ----------
run_beta_sweep <- function(dat, betas = c(0, 0.1, 0.3, 0.5, 0.7, 1.0),
                           alpha = 1, maxp = 2, restarts = 5) {
  # coerce all columns to factor
  dat[] <- lapply(dat, function(x) if (is.factor(x)) x else factor(x))
  
  res <- list()
  ref_dag <- NULL
  for (b in betas) {
    cat(sprintf("\n[hc] beta=%.3f | alpha=%.2f | maxp=%d | restarts=%d\n", b, alpha, maxp, restarts))
    wrap <- make_beta_wrapper(beta = b, alpha = alpha)
    dag  <- hc(dat, score = "custom-score", fun = wrap, maxp = maxp, restart = restarts)
    
    narcs <- narcs(dag)
    if (is.null(ref_dag)) {
      ref_dag <- dag
      ham <- 0
    } else {
      ham <- hamming(dag, ref_dag)
    }
    res[[as.character(b)]] <- list(dag = dag, narcs = narcs, hamming_vs_b0 = ham)
    cat("  arcs:", narcs, "| Hamming vs β=0:", ham, "\n")
  }
  res
}

# =========================
# Example datasets and tests
# =========================
library(bnlearn)

# --- Load your CHDS dataset ---
base_dir      <- "Desktop/betadiv_roy"
chds_path     <- file.path(base_dir, "CHDS.latentexample1.csv")

if (!file.exists(chds_path)) {
  stop("CHDS.latentexample1.csv not found at: ", chds_path)
}

chds_orig <- read.csv(chds_path, stringsAsFactors = TRUE)
chds_orig[] <- lapply(chds_orig, function(x) if (is.factor(x)) x else factor(x))

# --- Pick node and parents for testing ---
target_node <- "Admission"
target_parents <- c("Social", "Economic", "Events")

# --- Run print_node_score for KL and beta-div versions ---
cat("\n=== CHDS | β = 0 (Dirichlet-MLL / KL limit) ===\n")
print_node_score(chds_orig, node = target_node,
                 parents = target_parents,
                 beta = 0, alpha = 1)

cat("\n=== CHDS | β = 0.3 (Beta-divergence) ===\n")
print_node_score(chds_orig, node = target_node,
                 parents = target_parents,
                 beta = 0.3, alpha = 1)

