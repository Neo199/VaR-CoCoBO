Contamination <- function(x, runlength, seed) {
  # Initialize outputs
  FnVar <- NA
  FnGrad <- NA
  FnGradCov <- NA
  ConstraintGrad <- NA
  ConstraintGradCov <- NA
  
  n <- length(x)
  # Input validation
  if (length(x) != n || any(x > 1) || any(x < 0) || runlength <= 0 || seed <= 0 || seed != round(seed)) {
    cat(sprintf("x has %d elements, elements of x are binary, 
                 runlength should be positive and real, seed should be a positive integer.\n", n))
    return(list(fn = NA, constraint = NA, ConstraintCov = NA))
  } else {
    # Parameters
    nGen <- runlength         # number of independent generations
    u <- x                    # prevention binary decision variable
    X <- matrix(0, n, nGen)   # fraction contaminated at each stage for each generation
    epsilon <- rep(0.05, n)   # error probability
    p <- rep(0.1, n)          # proportion limit
    cost <- rep(1, n)         # cost for prevention at stage i
    
    # Beta parameters for initial contamination, contamination rate, restoration rate
    initialAlpha <- 1
    initialBeta <- 30
    contamAlpha <- 1
    contamBeta <- 17 / 3
    restoreAlpha <- 1
    restoreBeta <- 3 / 7
    
    # Set random seed
    set.seed(seed)
    
    # Generate initial fraction of contamination
    initialX <- rbeta(nGen, initialAlpha, initialBeta)
    
    # Generate rates of contamination and restoration
    Lambda <- matrix(rbeta(n * nGen, contamAlpha, contamBeta), n, nGen)
    Gamma <- matrix(rbeta(n * nGen, restoreAlpha, restoreBeta), n, nGen)
    
    # Determine contamination fractions
    X[1, ] <- Lambda[1, ] * (1 - u[1]) * (1 - initialX) + (1 - Gamma[1, ] * u[1]) * initialX
    for (i in 2:n) {
      X[i, ] <- Lambda[i, ] * (1 - u[i]) * (1 - X[i - 1, ]) + (1 - Gamma[i, ] * u[i]) * X[i - 1, ]
    }
    
    # Limit and cost of contamination control
    limit <- 1 - epsilon
    fn <- sum(cost * u)
    
    # Constraint checking
    con <- matrix(0, nGen, n)
    for (j in 1:nGen) {
      con[j, ] <- X[, j] <= p
    }
    le <- sum(rowSums(con) == n)
    constraint <- rep(0, n)
    for (k in 1:n) {
      constraint[k] <- (sum(con[, k]) / runlength) - limit[k]
    }
    ConstraintCov <- cov(con)
    
    return(list(fn = fn, constraint = constraint, ConstraintCov = ConstraintCov))
  }
}


# x <- c(1, 0, 1)  # Binary decision vector
# runlength <- 1000
# seed <- 42
# result <- Contamination(x, runlength, seed)
# 
# cat("Cost (fn):", result$fn, "\n")
# cat("Constraints:", result$constraint, "\n")
# cat("Constraint Covariance Matrix:\n")
# print(result$ConstraintCov)


contamination_prob <- function(x, n_samples, seed) {
  # Declare gamma factor (Lagrange constants)
  gamma <- 10
  
  # Find the total number of input samples
  num_inputs <- nrow(x)
  out <- numeric(num_inputs)
  
  # Iterate over each input sample
  for (i in 1:num_inputs) {
    # Run contamination study
    contamination_result <- Contamination(x[i, ], n_samples, seed)
    cost <- contamination_result$fn
    constraint <- contamination_result$constraint
    
    # Compute total output
    out[i] <- cost - sum(gamma * constraint)
  }
  
  return(out)
}

# 
# # Define input matrix (binary variables)
# x <- matrix(c(0, 1, 1, 0, 1, 0), nrow = 3, byrow = TRUE)
# 
# # Parameters
# n_samples <- 1000
# seed <- 5
# 
# # Compute contamination probabilities
# out <- contamination_prob(x, n_samples, seed)
# 
# print(out)
# 
