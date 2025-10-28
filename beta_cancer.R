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

make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}

get_cancer_data <- function(n = 5000) {
  have_fit <- FALSE
  suppressWarnings({
    if ("package:bnlearn" %in% search() || requireNamespace("bnlearn", quietly = TRUE)) {
      # Some bnlearn builds export a bn.fit named 'cancer'
      if ("cancer" %in% data(package = "bnlearn")$results[, "Item"]) {
        data("cancer", package = "bnlearn")
        # If loaded object is in env and is a bn.fit, we can sample from it:
        if (exists("cancer") && inherits(get("cancer"), "bn.fit")) {
          have_fit <- TRUE
          return(rbn(get("cancer"), n = n))
        }
      }
    }
  })
  if (!have_fit) {
    levYN <- c("no","yes")
    Pollution <- factor(sample(c("low","high"), n, replace = TRUE, prob = c(0.9, 0.1)))
    Smoker    <- factor(sample(levYN, n, replace = TRUE, prob = c(0.5, 0.5)))
    
    pc <- mapply(function(p, s) {
      if (p == "low"  && s == "no")  0.03 else
        if (p == "low"  && s == "yes") 0.05 else
          if (p == "high" && s == "no")  0.10 else
            0.20
    }, as.character(Pollution), as.character(Smoker))
    
    Cancer <- factor(ifelse(runif(n) < pc, "yes", "no"), levels = levYN)
    
    # Xray: P(X=yes | C): 0.9 if cancer, 0.2 otherwise
    px <- ifelse(Cancer == "yes", 0.9, 0.2)
    Xray <- factor(ifelse(runif(n) < px, "yes", "no"), levels = levYN)
    
    # Dyspnoea: P(D=yes | C): 0.65 if cancer, 0.30 otherwise
    pd <- ifelse(Cancer == "yes", 0.65, 0.30)
    Dyspnoea <- factor(ifelse(runif(n) < pd, "yes", "no"), levels = levYN)
    
    data.frame(Pollution, Smoker, Cancer, Xray, Dyspnoea, check.names = FALSE)
  }
}

cancer_df <- get_cancer_data(n = 5000)
str(cancer_df)   # sanity check: all factors

wrapper0  <- make_beta_wrapper(beta = 0.0,  alpha = 1)
wrapper01 <- make_beta_wrapper(beta = 0.01, alpha = 1)
wrapper02 <- make_beta_wrapper(beta = 0.00000001, alpha = 1)

net0  <- hc(cancer_df, score = "custom-score", fun = wrapper0,  maxp = 2)
net01 <- hc(cancer_df, score = "custom-score", fun = wrapper01, maxp = 2)
net02 <- hc(cancer_df, score = "custom-score", fun = wrapper02, maxp = 2)


plot(net0,  main = "CANCER | beta = 0.0")
plot(net01, main = "CANCER | beta = 0.01")
plot(net02, main = "CANCER | beta = 0.00000001")


