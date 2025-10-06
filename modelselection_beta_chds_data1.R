library(bnlearn)

set.seed(42)


chds_data  <- read.csv("Desktop/betaDiv/CHDS.latentexample1.csv", stringsAsFactors = TRUE)
data1_full      <- read.csv("Desktop/betaDiv/data1.csv", stringsAsFactors = TRUE)
data1      <- data1_full[1:200, 1:6]

# Helper: create wrapper with beta
make_beta_wrapper <- function(beta, alpha = 1) {
  function(node, parents, data, args) {
    beta_div_score(node, parents, data, beta = beta, alpha = alpha)
  }
}

# Function: run model selection
run_beta_sweep <- function(data, betas, maxp = 2) {
  results <- list()
  empty_net <- empty.graph(nodes = names(data))
  
  for (b in betas) {
    cat(sprintf("Running hill-climbing for beta = %.3f...\n", b))
    wrapper <- make_beta_wrapper(beta = b)
    
    net <- hc(
      data,
      score = "custom-score",
      fun   = wrapper,
      start = empty_net,
      maxp  = maxp,
      debug = FALSE
    )
    
    results[[as.character(b)]] <- net
    cat(sprintf("Beta = %.3f: %d arcs\n", b, nrow(net$arcs)))
  }
  
  return(results)
}

# Run sweeps
beta_seq <- seq(0, 1, by = 0.1)

cat("\n=== CHDS Data Sweep ===\n")
chds_results <- run_beta_sweep(chds_data, beta_seq)

cat("\n=== data1.csv (reduced) Sweep ===\n")
data1_results <- run_beta_sweep(data1, beta_seq)

# Summarize arc counts
data1_arc_counts <- sapply(data1_results, function(net) nrow(net$arcs))



cat("\nArc counts (data1 reduced):\n")
print(data1_arc_counts)



plot(data1_results[["0.4"]], main = "data1 reduced (beta = 0.4)")
plot(data1_results[["0.5"]], main = "data1 reduced (beta = 0.5)")
