
library(bnlearn)
set.seed(42)


chds_path <- "Desktop/betaDiv/CHDS.latentexample1.csv"   

#  beta_div_score 
beta_div_score <- function(node, parents, data, beta = 0, alpha = 1.0, ...) {
  stopifnot(is.character(node), length(node) == 1, node %in% names(data))
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents) == 0) parents <- character(0)
  
  ll_dirichlet <- function(counts, alpha_vec) {
    A <- sum(alpha_vec); N <- sum(counts)
    lgamma(A) - lgamma(A + N) + sum(lgamma(alpha_vec + counts) - lgamma(alpha_vec))
  }
  
  # beta = 0 -> log-marginal likelihood
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
  
  # beta != 0 -> beta-divergence using posterior predictive p_hat
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

# ---- Your wrapper factory ----
make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}


# beta = 0.0
wrapper0 <- make_beta_wrapper(beta = 0.0, alpha = 1)
net0 <- hc(chds_orig, score = "custom-score", fun = wrapper0)
plot(net0, main = "beta = 0.0")

# beta = 0.1
wrapper01 <- make_beta_wrapper(beta = 0.1, alpha = 1)
net01 <- hc(chds_orig, score = "custom-score", fun = wrapper01)
plot(net01, main = "beta = 0.1")

# beta = 0.2
wrapper02 <- make_beta_wrapper(beta = 0.2, alpha = 1)
net02 <- hc(chds_orig, score = "custom-score", fun = wrapper02)
plot(net02, main = "beta = 0.2")

# beta = 0.3
wrapper03 <- make_beta_wrapper(beta = 0.3, alpha = 1)
net03 <- hc(chds_orig, score = "custom-score", fun = wrapper03)
plot(net03, main = "beta = 0.3")

# beta = 0.4
wrapper04 <- make_beta_wrapper(beta = 0.4, alpha = 1)
net04 <- hc(chds_orig, score = "custom-score", fun = wrapper04)
plot(net04, main = "beta = 0.4")

# beta = 0.5
wrapper05 <- make_beta_wrapper(beta = 0.5, alpha = 1)
net05 <- hc(chds_orig, score = "custom-score", fun = wrapper05)
plot(net05, main = "beta = 0.5")

# beta = 0.6
wrapper06 <- make_beta_wrapper(beta = 0.6, alpha = 1)
net06 <- hc(chds_orig, score = "custom-score", fun = wrapper06)
plot(net06, main = "beta = 0.6")

# beta = 0.7
wrapper07 <- make_beta_wrapper(beta = 0.7, alpha = 1)
net07 <- hc(chds_orig, score = "custom-score", fun = wrapper07)
plot(net07, main = "beta = 0.7")

# beta = 0.8
wrapper08 <- make_beta_wrapper(beta = 0.8, alpha = 1)
net08 <- hc(chds_orig, score = "custom-score", fun = wrapper08)
plot(net08, main = "beta = 0.8")

# beta = 0.9
wrapper09 <- make_beta_wrapper(beta = 0.9, alpha = 1)
net09 <- hc(chds_orig, score = "custom-score", fun = wrapper09)
plot(net09, main = "beta = 0.9")

library(bnlearn)

# --- Load the 90% corrupted dataset ---
chds_90_path <- "Desktop/betadiv_roy//CHDS_corrupted_boosted.csv"
chds_corrupted <- read.csv(chds_90_path, stringsAsFactors = TRUE)

# --- Run hc() individually for beta = 0.0 → 0.9 ---

# beta = 0.0
wrapper0 <- make_beta_wrapper(beta = 0.0, alpha = 1)
net0c <- hc(chds_corrupted, score = "custom-score", fun = wrapper0)
plot(net0c, main = "CHDS 90% corrupted | beta = 0.0")

# beta = 0.1
wrapper01 <- make_beta_wrapper(beta = 0.1, alpha = 1)
net01c <- hc(chds_corrupted, score = "custom-score", fun = wrapper01)
plot(net01c, main = "CHDS 90% corrupted | beta = 0.1")

# beta = 0.2
wrapper02 <- make_beta_wrapper(beta = 0.2, alpha = 1)
net02c <- hc(chds_corrupted, score = "custom-score", fun = wrapper02)
plot(net02c, main = "CHDS 90% corrupted | beta = 0.2")

# beta = 0.3
wrapper03 <- make_beta_wrapper(beta = 0.3, alpha = 1)
net03c <- hc(chds_corrupted, score = "custom-score", fun = wrapper03)
plot(net03c, main = "CHDS 90% corrupted | beta = 0.3")

# beta = 0.4
wrapper04 <- make_beta_wrapper(beta = 0.4, alpha = 1)
net04c <- hc(chds_corrupted, score = "custom-score", fun = wrapper04)
plot(net04c, main = "CHDS 90% corrupted | beta = 0.4")

# beta = 0.5
wrapper05 <- make_beta_wrapper(beta = 0.5, alpha = 1)
net05c <- hc(chds_corrupted, score = "custom-score", fun = wrapper05)
plot(net05c, main = "CHDS 90% corrupted | beta = 0.5")

# beta = 0.6
wrapper06 <- make_beta_wrapper(beta = 0.6, alpha = 1)
net06c <- hc(chds_corrupted, score = "custom-score", fun = wrapper06)
plot(net06c, main = "CHDS 90% corrupted | beta = 0.6")

# beta = 0.7
wrapper07 <- make_beta_wrapper(beta = 0.7, alpha = 1)
net07c <- hc(chds_corrupted, score = "custom-score", fun = wrapper07)
plot(net07c, main = "CHDS 90% corrupted | beta = 0.7")

# beta = 0.8
wrapper08 <- make_beta_wrapper(beta = 0.8, alpha = 1)
net08c <- hc(chds_corrupted, score = "custom-score", fun = wrapper08)
plot(net08c, main = "CHDS 90% corrupted | beta = 0.8")

# beta = 0.9
wrapper09 <- make_beta_wrapper(beta = 0.9, alpha = 1)
net09c <- hc(chds_corrupted, score = "custom-score", fun = wrapper09)
plot(net09c, main = "CHDS 90% corrupted | beta = 0.9")



