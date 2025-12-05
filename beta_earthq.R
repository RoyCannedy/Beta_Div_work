library(bnlearn)
set.seed(42)

ensure_all_factors <- function(df) {
  for (nm in names(df)) if (!is.factor(df[[nm]])) df[[nm]] <- factor(df[[nm]])
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
      cfg <- parent_configs[t]
      if (t == 1) {
        prior <- setNames(rep(0, n_values), node_values)
      } else {
        idx <- parent_configs[1:(t-1)] == cfg
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
  total_score
}

# --- Sequential β-divergence scorer (same math) ---
beta_div_score <- function(node, parents, data, beta = 0, alpha = 1.0, ...) {
  stopifnot(is.character(node), length(node) == 1, node %in% names(data))
  data <- ensure_all_factors(data)
  if (length(parents) == 0) parents <- character(0)
  
  if (beta == 0 || beta < 1e-8) {
    return(beta_div_limit_zero(node, parents, data, alpha))
  }
  
  node_values <- levels(data[[node]])
  n_values <- length(node_values)
  total_score <- 0
  
  if (length(parents) == 0) {
    # One global running prior
    for (t in seq_len(nrow(data))) {
      prior <- if (t == 1) setNames(rep(0, n_values), node_values)
      else table(factor(data[[node]][1:(t-1)], levels = node_values))
      y_bar <- sum(prior)
      alpha_bar <- alpha * n_values
      if (y_bar + alpha_bar > 0) {
        p_vec <- (as.numeric(prior) + alpha) / (y_bar + alpha_bar)
        # Tsallis / β-divergence piece
        total_score <- total_score + (sum(p_vec^beta)/beta - sum(p_vec^(beta+1))/(beta+1))
      }
    }
  } else {
    parent_configs <- do.call(paste, c(data[parents], sep = "_"))
    for (t in seq_len(nrow(data))) {
      cfg <- parent_configs[t]
      if (t == 1) {
        prior <- setNames(rep(0, n_values), node_values)
      } else {
        idx <- parent_configs[1:(t-1)] == cfg
        if (any(idx)) prior <- table(factor(data[[node]][1:(t-1)][idx], levels = node_values))
        else prior <- setNames(rep(0, n_values), node_values)
      }
      y_bar <- sum(prior)
      alpha_bar <- alpha * n_values
      if (y_bar + alpha_bar > 0) {
        p_vec <- (as.numeric(prior) + alpha) / (y_bar + alpha_bar)
        total_score <- total_score + (sum(p_vec^beta)/beta - sum(p_vec^(beta+1))/(beta+1))
      }
    }
  }
  
  total_score
}

make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}

# --- EARTHQUAKE dataset (try bnlearn, else synthesize classic BN) ---
get_earthquake_data <- function(n = 5000) {
  have_obj <- FALSE
  
  suppressWarnings({
    if ("package:bnlearn" %in% search() || requireNamespace("bnlearn", quietly = TRUE)) {
      if ("earthquake" %in% data(package = "bnlearn")$results[, "Item"]) {
        data("earthquake", package = "bnlearn")
        obj <- get("earthquake")
        have_obj <- TRUE
        
        if (is.data.frame(obj)) return(obj)
        
        if (inherits(obj, "bn.fit")) return(rbn(obj, n = n))
      }
    }
  })
  
  if (!have_obj) {
    levYN <- c("no", "yes")
    
    
    B <- rbinom(n, 1, 0.001)  # P(B = yes) = 0.001
    E <- rbinom(n, 1, 0.002)  # P(E = yes) = 0.002
    
    pA <- numeric(n)
    for (i in seq_len(n)) {
      if (B[i] == 1 && E[i] == 1)      pA[i] <- 0.95
      else if (B[i] == 1 && E[i] == 0) pA[i] <- 0.94
      else if (B[i] == 0 && E[i] == 1) pA[i] <- 0.29
      else                             pA[i] <- 0.001
    }
    A <- rbinom(n, 1, pA)
    
    pJ <- ifelse(A == 1, 0.90, 0.05)
    J  <- rbinom(n, 1, pJ)
    

    pM <- ifelse(A == 1, 0.70, 0.01)
    M  <- rbinom(n, 1, pM)
    
    data.frame(
      Burglary  = factor(ifelse(B == 1, "yes", "no"), levels = levYN),
      Earthquake = factor(ifelse(E == 1, "yes", "no"), levels = levYN),
      Alarm     = factor(ifelse(A == 1, "yes", "no"), levels = levYN),
      JohnCalls = factor(ifelse(J == 1, "yes", "no"), levels = levYN),
      MaryCalls = factor(ifelse(M == 1, "yes", "no"), levels = levYN),
      check.names = FALSE
    )
  }
}

earthquake_df <- ensure_all_factors(get_earthquake_data(n = 5000))
str(earthquake_df)

wrapper0   <- make_beta_wrapper(beta = 0.0,  alpha = 1)
wrapper009 <- make_beta_wrapper(beta = 0.09, alpha = 1)
wrapper090 <- make_beta_wrapper(beta = 0.40, alpha = 1)

net_eq_b0   <- hc(earthquake_df, score = "custom-score", fun = wrapper0,   maxp = 2)
net_eq_b009 <- hc(earthquake_df, score = "custom-score", fun = wrapper009, maxp = 2)
net_eq_b090 <- hc(earthquake_df, score = "custom-score", fun = wrapper090, maxp = 2)

# --- Plots ---
plot(net_eq_b0,   main = "EARTHQUAKE | beta = 0.00")
plot(net_eq_b009, main = "EARTHQUAKE | beta = 0.09")
plot(net_eq_b090, main = "EARTHQUAKE | beta = 0.40")
