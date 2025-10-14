library(bnlearn)
set.seed(42)

# --- Paths that make sense together ---
base_dir      <- "Desktop/betadiv_roy/"
chds_path     <- file.path(base_dir, "CHDS.latentexample1.csv")
corrupted_dir <- file.path(base_dir, "CHDS_corrupted")  # where CHDS_corrupted_XXperc.csv live

chds_orig <- read.csv(chds_path, stringsAsFactors = TRUE)

#=beta_div_score
beta_div_score <- function(node, parents, data, beta = 0, alpha = 1.0, ...) {
  stopifnot(is.character(node), length(node) == 1, node %in% names(data))
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents) == 0) parents <- character(0)
  
  ll_dirichlet <- function(counts, alpha_vec) {
    A <- sum(alpha_vec); N <- sum(counts)
    lgamma(A) - lgamma(A + N) + sum(lgamma(alpha_vec + counts) - lgamma(alpha_vec))
  }
  
  if (isTRUE(all.equal(beta, 0))) {
    if (length(parents) == 0) {
      lev <- levels(data[[node]]); K <- length(lev)
      alpha_vec <- rep(alpha / K, K)
      tab <- table(data[[node]])
      y <- setNames(rep(0L, K), lev); y[names(tab)] <- as.integer(tab)
      return(ll_dirichlet(y, alpha_vec))
    } else {
      lev <- levels(data[[node]]); K <- length(lev)
      alpha_vec <- rep(alpha / K, K)
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
  
  if (length(parents) == 0) {
    lev <- levels(data[[node]]); K <- length(lev)
    tab <- table(data[[node]])
    y <- setNames(rep(0L, K), lev); y[names(tab)] <- as.integer(tab)
    N <- sum(y)
    p_hat <- (y + alpha / K) / (N + alpha)
    return(sum(p_hat^beta)/beta - sum(p_hat^(beta+1))/(beta+1))
  } else {
    lev <- levels(data[[node]]); K <- length(lev)
    f <- as.formula(paste("~", paste(c(parents, node), collapse = "+")))
    ct <- xtabs(f, data = data)
    margins <- seq_len(length(dim(ct)) - 1)
    A <- alpha
    return(sum(apply(ct, margins, function(slice) {
      y <- as.integer(slice); N <- sum(y)
      p_hat <- (y + A / K) / (N + A)
      sum(p_hat^beta)/beta - sum(p_hat^(beta+1))/(beta+1)
    })))
  }
}

# ---- Wrapper factory ----
make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}

wrapper0  <- make_beta_wrapper(beta = 0.0, alpha = 1)
wrapper01 <- make_beta_wrapper(beta = 0.1, alpha = 1)
net0  <- hc(chds_orig, score = "custom-score", fun = wrapper0)
plot(net0,  main = "ORIGINAL | beta = 0.0")
net01 <- hc(chds_orig, score = "custom-score", fun = wrapper01)
plot(net01, main = "ORIGINAL | beta = 0.1")

corrupted_dir <- "CHDS_corrupted" 

wrapper0  <- make_beta_wrapper(beta = 0.0, alpha = 1)
wrapper01 <- make_beta_wrapper(beta = 0.1, alpha = 1)

for (p in seq(0, 100, by = 10)) {
  chds_corrupted_path <- file.path(corrupted_dir, sprintf("CHDS_corrupted_%dperc.csv", p))
  chds_corrupted <- read.csv(chds_corrupted_path, stringsAsFactors = TRUE)
  
  # beta = 0.0
  net0c <- hc(chds_corrupted, score = "custom-score", fun = wrapper0)
  plot(net0c,  main = sprintf("CHDS %d%% corrupted | beta = 0.0", p))
  
  # beta = 0.1
  net01c <- hc(chds_corrupted, score = "custom-score", fun = wrapper01)
  plot(net01c, main = sprintf("CHDS %d%% corrupted | beta = 0.1", p))
}



