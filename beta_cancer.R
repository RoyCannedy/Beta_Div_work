library(bnlearn)
set.seed(42)

ensure_all_factors <- function(df) {
  for (nm in names(df)) {
    if (!is.factor(df[[nm]])) df[[nm]] <- factor(df[[nm]])
  }
  df
}

beta_div_limit_zero <- function(node, parents, data, alpha = 1.0) {
  stopifnot(is.character(node), length(node) == 1, node %in% names(data))
  data <- ensure_all_factors(data)
  if (length(parents) == 0) parents <- character(0)
  
  node_values <- levels(data[[node]])
  n_values <- length(node_values)
  total_score <- 0
  
  if (length(parents) == 0) {
    for (t in seq_len(nrow(data))) {
      obs <- data[[node]][t]
      prior <- if (t == 1) setNames(rep(0, n_values), node_values)
      else table(factor(data[[node]][1:(t-1)], levels = node_values))
      y_bar <- sum(prior)
      alpha_bar <- alpha * n_values
      if (y_bar + alpha_bar > 0) {
        p_vec <- (as.numeric(prior) + alpha) / (y_bar + alpha_bar)
        names(p_vec) <- node_values
        total_score <- total_score + log(p_vec[as.character(obs)])
      }
    }
  } else {
    parent_configs <- do.call(paste, c(data[parents], sep = "_"))
    for (t in seq_len(nrow(data))) {
      obs <- data[[node]][t]
      config <- parent_configs[t]
      if (t == 1) {
        prior <- setNames(rep(0, n_values), node_values)
      } else {
        idx <- parent_configs[1:(t-1)] == config
        if (any(idx)) prior <- table(factor(data[[node]][1:(t-1)][idx], levels = node_values))
        else prior <- setNames(rep(0, n_values), node_values)
      }
      y_bar <- sum(prior)
      alpha_bar <- alpha * n_values
      if (y_bar + alpha_bar > 0) {
        p_vec <- (as.numeric(prior) + alpha) / (y_bar + alpha_bar)
        names(p_vec) <- node_values
        total_score <- total_score + log(p_vec[as.character(obs)])
      }
    }
  }
  return(total_score)
}


beta_div_score <- function(node, parents, data, beta = 1.0, alpha = 1.0) {
  stopifnot(is.character(node), length(node) == 1, node %in% names(data))
  data <- ensure_all_factors(data)
  if (!is.factor(data[[node]])) data[[node]] <- factor(data[[node]])
  if (length(parents) == 0) parents <- character(0)
  
  if (beta == 0 || beta < 1e-4) {
    return(beta_div_limit_zero(node, parents, data, alpha))
  }
  
  node_values <- levels(data[[node]])
  n_values <- length(node_values)
  total_score <- 0
  
  if (length(parents) == 0) {
    for (t in seq_len(nrow(data))) {
      obs <- data[[node]][t]
      prior <- if (t == 1) setNames(rep(0, n_values), node_values)
      else table(factor(data[[node]][1:(t-1)], levels = node_values))
      
      y_j <- prior[obs]
      y_bar <- sum(prior)
      alpha_bar <- alpha * n_values
      if (y_bar + alpha_bar > 0) {
        p_j <- (y_j + alpha) / (y_bar + alpha_bar)
        first_term <- (1 / beta) * (p_j^beta)
        second_term <- sum(sapply(node_values, function(val) {
          p_k <- (prior[val] + alpha) / (y_bar + alpha_bar)
          p_k^(beta + 1)
        })) / (beta + 1)
        total_score <- total_score + (first_term - second_term)
      }
    }
  } else {
    parent_configs <- do.call(paste, c(data[parents], sep = "_"))
    for (t in seq_len(nrow(data))) {
      obs <- data[[node]][t]
      config <- parent_configs[t]
      prior <- if (t == 1) setNames(rep(0, n_values), node_values)
      else {
        idx <- parent_configs[1:(t-1)] == config
        if (any(idx)) table(factor(data[[node]][1:(t-1)][idx], levels = node_values))
        else setNames(rep(0, n_values), node_values)
      }
      y_j <- prior[obs]
      y_bar <- sum(prior)
      alpha_bar <- alpha * n_values
      if (y_bar + alpha_bar > 0) {
        p_j <- (y_j + alpha) / (y_bar + alpha_bar)
        first_term <- (1 / beta) * (p_j^beta)
        second_term <- sum(sapply(node_values, function(val) {
          p_k <- (prior[val] + alpha) / (y_bar + alpha_bar)
          p_k^(beta + 1)
        })) / (beta + 1)
        total_score <- total_score + (first_term - second_term)
      }
    }
  }
  
  return(total_score)
}

# ---------- 2) Wrapper for hc ----------
make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}

# ---------- 3) Get a CANCER dataset ----------
get_cancer_data <- function(n = 5000) {
  have_fit <- FALSE
  suppressWarnings({
    if ("package:bnlearn" %in% search() || requireNamespace("bnlearn", quietly = TRUE)) {
      if ("cancer" %in% data(package = "bnlearn")$results[, "Item"]) {
        data("cancer", package = "bnlearn")
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
          if (p == "high" && s == "no")  0.10 else 0.20
    }, as.character(Pollution), as.character(Smoker))
    
    Cancer <- factor(ifelse(runif(n) < pc, "yes", "no"), levels = levYN)
    px <- ifelse(Cancer == "yes", 0.9, 0.2)
    Xray <- factor(ifelse(runif(n) < px, "yes", "no"), levels = levYN)
    pd <- ifelse(Cancer == "yes", 0.65, 0.30)
    Dyspnoea <- factor(ifelse(runif(n) < pd, "yes", "no"), levels = levYN)
    
    data.frame(Pollution, Smoker, Cancer, Xray, Dyspnoea, check.names = FALSE)
  }
}

cancer_df <- ensure_all_factors(get_cancer_data(n = 5000))
str(cancer_df)

# ---------- 4) Run HC with beta=0, 0.9, 0.01 ----------
wrapper0   <- make_beta_wrapper(beta = 0.0,  alpha = 1)
wrapper09  <- make_beta_wrapper(beta = 0.9,  alpha = 1)
wrapper001 <- make_beta_wrapper(beta = 0.01, alpha = 1)

net0   <- hc(cancer_df, score = "custom-score", fun = wrapper0,   maxp = 2)
net09  <- hc(cancer_df, score = "custom-score", fun = wrapper09,  maxp = 2)
net001 <- hc(cancer_df, score = "custom-score", fun = wrapper001, maxp = 2)

plot(net0,   main = "CANCER | beta = 0.0 (limit to KL)")
plot(net09,  main = "CANCER | beta = 0.9")
plot(net001, main = "CANCER | beta = 0.01")
