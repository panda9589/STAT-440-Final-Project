library(Matrix)
library(MASS)
library(stats4)
library(ggplot2)


# Read in data
y <- scan(file.path("final-project-data.txt"))
n <- length(y)

latent_pois_pdf <- function(yi, lambda) {
  lambda^yi * exp(-1 * lambda) / factorial(yi)
}

# here y, z are vector, beta0 and sigma are scalar
# params = c(beta0, sigma)

set.seed(20898576)
L = 500
Q <- 2 * Diagonal(n) # Sparse matrix storage
for (i in 1:(n-1)) {
  Q[i+1, i] <- -1
  Q[i, i+1] <- -1
}
Q[n, n] <- 1

generate_zs <- function(n, Q, L) {
  # Check if Q is positive-definite
  if (any(eigen(Q)$values <= 0)) {
    stop("The precision matrix Q is not positive-definite.")
  }
  
  # Cholesky decomposition of the precision matrix Q
    L_chol <- solve(chol(Q)) # Q = LL^T
  
  # Initialize a matrix to store the generated z vectors
  zs <- matrix(0, nrow = n, ncol = L)
  
  # Generate each z vector individually
  for (i in 1:L) {
    # Generate standard normal random variables
    w <- rnorm(n, 0, 1)
    
    # Calculate z = solve(L_chol) %*% w
    zs[, i] <- solve(L_chol, w)
  }
  
  return(zs)
}

# Example usage
zs <- generate_zs(n, Q, L)

# Extracting the 5th simulation of z
z_5 <- zs[, 5]
print(z_5)

# Checking the range and standard deviation
range_z1 <- range(zs[, 1])
range_z1000 <- range(zs[, L])
sd_z1 <- sd(zs[, 1])
sd_z1000 <- sd(zs[, L])

print(range_z1)
print(range_z1000)
print(sd_z1)
print(sd_z1000)


prod_pois <- function(params, z) {
  beta0 = params[1]
  sigma = params[2]
  lambda = exp(beta0 + sigma * z)
  len <- length(y)
  ans <- 1
  for (i in 1:len) {
    cur <- latent_pois_pdf(y[i], lambda[i])
    ans <- ans * cur
  }
  return(ans)
}

# create_D <- function(n) {
#   # Create an n x n identity matrix
#   mat <- diag(1, n, n)
#   # Modify elements to the left of the diagonal starting from the second row
#   for (i in 2:n) {
#     mat[i, i-1] <- -1
#   }
#   return(mat)
# }
# # Example usage for n = 5
# D <- create_D(5)
# print(D)
# t(D) %*% D

likelihood_func_pois <- function(params) {
  beta0 = params[1]
  sigmasq = params[2]
  ans = 0
  for (i in 1:L) {
    ans <- ans + sum(dpois(y, exp(beta0 + sqrt(sigmasq) * zs[, i]), log = TRUE))
    # ans <- ans + prod_pois(params, zs[, i])
  }
  ans - log(L)
}

nll <- function(params) (-likelihood_func_pois(params))



# # Define a sequence of beta0 values and sigma values
# beta0_values <- seq(-20, 20, by = 0.5)
# sigma_values <- seq(-20, 20, by = 0.5)
# 
# # Calculate NLL for varying beta0 and fixed sigma = 0.1
# nll_beta0 <- sapply(beta0_values, function(beta0) nll(c(beta0, 0.1)))
# 
# # Calculate NLL for varying sigma and fixed beta0 = 0
# nll_sigma <- sapply(sigma_values, function(sigma) nll(c(0, sigma)))
# 
# Plot NLL as a function of beta0 for fixed sigma = 0.1
# plot(beta0_values, nll_beta0, type = "l", col = "blue",
#      main = "NLL as a function of beta0 (sigma = 0.1)",
#      xlab = "beta0", ylab = "NLL")
# 
# # Plot NLL as a function of sigma for fixed beta0 = 0
# plot(sigma_values, nll_sigma, type = "l", col = "red",
#      main = "NLL as a function of sigma (beta0 = 0)",
#      xlab = "sigma", ylab = "NLL")
# 
# 
# library(car)
# 
# # Define the range of values for beta0 and sigma
# sigma_vals <- seq(-1, 1, length.out = 20)
# beta0_vals <- seq(-6, 1, length.out = 20)
# 
# # Create a grid of values
# grid <- expand.grid(beta0 = beta0_vals, sigma = sigma_vals)
# 
# # Compute nll values for each combination
# grid$nll <- apply(grid, 1, function(x) nll(c(x['beta0'], x['sigma'])))
# 
# 
# # 3D scatter plot
# scatter3d(x = grid$beta0, y = grid$sigma, z = grid$nll, 
#           xlab = "Beta0", ylab = "Sigma", zlab = "NLL", 
#           point.col = "blue", surface = FALSE) # Set surface=TRUE to see a fitted surface as well
# 


# ***** NEWTON METHOD (find beta0, sigma) ******************************************************************************************************************************


nllg <- function(params) tryCatch(numDeriv::grad(nll, params), error = function(e) return)
nllh <- function(params) tryCatch(numDeriv::hessian(nll, params, method.args = 
                                                      list(eps=1e-3, d=0.3, zero.tol=sqrt(.Machine$double.eps/7e-7), r=3, v=2, show.details=FALSE)), error = function(e) return)


# Starting values
beta0_start = 0
sigma_start = 0.1
theta_start <- c(beta0_start, sigma_start)
# temp = optim(theta_start, nll, "BFGS")

# Set up the algorithm
run_newton <- function (reflect = FALSE, halving = FALSE) {
  theta = theta_start
  history <- list(); history[[1]] <- theta
  steps <- list(); steps[[1]] <- theta
  maxstep <- 1
  itr <- 1
  eps <- 1e-02
  while (max(abs(nllg(theta))) > eps) {
    step = solve(nllh(theta), nllg(theta))
    ## reflection on boundary: ----------------------
    if (reflect) {
      step = solve(abs(nllh(theta)), nllg(theta))
      # if(theta < -1/max(lambda)){
      #   theta <- (-2/max(lambda)) -theta
      # }
    }
    steps[[itr + 1]] <- step
    # step halving ----------------------------------
    if (halving) {
      while (min(step) > maxstep) step = step / 2
    }
    # the actual iteration **************************
    theta <- theta - step
    itr <- itr + 1
    history[[itr]] <- theta
    print(theta)
  }
  return (list(theta, history))
}
# 
results_neither = run_newton(FALSE, FALSE)
# results_halving = run_newton(FALSE, TRUE)
# results_reflect = run_newton(TRUE, FALSE)
# results_both = run_newton(TRUE, TRUE)


# Extracting final theta values
# final_theta_neither <- results_neither[[1]]
# final_theta_halving <- results_halving[[1]]
# final_theta_reflect <- results_reflect[[1]]
final_theta_both <- results_both[[1]]

# Evaluate the NLL for each final theta
nll_both <- nll(final_theta_both)
# 
# # Organize NLL values and corresponding strategies into a data frame
# nll_values <- data.frame(
#   Strategy = c("Neither", "Halving", "Reflect", "Both"),
#   NLL = c(nll_neither, nll_halving, nll_reflect, nll_both),
#   Beta0 = c(final_theta_neither[1], final_theta_halving[1], final_theta_reflect[1], final_theta_both[1]),
#   Sigma = c(final_theta_neither[2], final_theta_halving[2], final_theta_reflect[2], final_theta_both[2])
# )
# 
# # Find the minimum NLL value and corresponding strategy
# min_nll <- nll_values[which.min(nll_values$NLL),]
# 
# # Print the results
# print(paste("The strategy with the minimum NLL is:", min_nll$Strategy))
# print(paste("Minimum NLL value:", min_nll$NLL))
# print(paste("Beta0 at minimum:", min_nll$Beta0))
# print(paste("Sigma at minimum:", min_nll$Sigma))
# 
# param_estimate <- c(Beta0 = min_nll$Beta0, Sigma = min_nll$Sigma)
# 
param_estimate <- final_theta_both
H <- nllh(param_estimate)

usd <- sqrt(diag(solve(H))) # solve(H) is nxn but NOT SPARSE. Bad!

# Calculate the inverse of the Hessian to get the covariance matrix
cov_matrix <- solve(H)

# Extract variances from the diagonal of the covariance matrix
var_beta0 <- cov_matrix[1, 1]
var_sigma <- cov_matrix[2, 2]

# Calculate standard deviations
sd_beta0 <- sqrt(var_beta0)
sd_sigma <- sqrt(var_sigma)

# Calculate confidence intervals (CIs) using z-value approx 1.98
ci_beta0 <- c(param_estimate[1] - 2 * sd_beta0, param_estimate[1] + 2 * sd_beta0)
ci_sigma <- c(param_estimate[2] - 2 * sd_sigma, param_estimate[2] + 2 * sd_sigma)

# Print standard deviations
cat("Standard deviation for beta0:", sd_beta0, "\n")
cat("Standard deviation for sigma:", sd_sigma, "\n")

# Print confidence intervals
cat("95% Confidence interval for beta0:", paste(round(ci_beta0, 2), collapse = " to "), "\n")
cat("95% Confidence interval for sigma:", paste(round(ci_sigma, 2), collapse = " to "), "\n")





beta0 <- param_estimate[1]
sigmasq <- param_estimate[2]^2
param = c(beta0, sigmasq^0.5)

nll_u <- function(u) 
  -sum(dpois(y, exp(u + beta0), log = TRUE)) + (sum(diff(u)^2) + (u[1] - 0)^2) / (-2 * sigmasq)
Q <- 2 * Diagonal(n) # Sparse matrix storage
for (i in 1:(n-1)) {
  Q[i+1, i] <- -1
  Q[i, i+1] <- -1
}
# x <- seq(-10, 10, by = 0.1)
# nll_values <- sapply(x, function(u) nll(u))
# 
# plot(x, nll_values, type = 'l', col = 'blue', main = 'NLL vs. x', xlab = 'x', ylab = 'NLL')


Q[n, n] <- 1
Q[1:5, 1:5] # Show the structure of Q
nllg <- function(u) {
  val <- -y + exp(u + beta0) + (1/sigmasq) * as.numeric(Q %*% u)
  val[1] <- val[1] - beta0 / sigmasq
  val
}
# Hessian
nllh <- function(u) (1/sigmasq) * Q + Diagonal(n, x = exp(u + beta0))
# 
# nll_u_V <- Vectorize(nll_u)
# 
# nllg <- function(u) tryCatch(numDeriv::grad(nll_u_V, u, method = "simple"), error = function(e) return)
# nllh <- function(u) tryCatch(numDeriv::hessian(nll_u_V, u, method = "complex"), error = function(e) return)

# Newton
u <- rep(beta0, n) # Starting value equal to the marginal mean
tmpgrad <- nllg(u)
eps <- 1e-03 # Convergence tolerance
maxitr <- 1000
itr <- 1
maxstep <- 2
gradvals <- numeric(maxitr)

# temp_u = optim(u, nll, "L-BFGS-B")

# Run the algorithm
while(max(abs(tmpgrad)) > eps & itr < maxitr) {
  # Newton step
  tmpgrad <- nllg(u) # Save the gradient for diagnostics
  print("here")
  tmphess <- nllh(u)
  print("heree")
  step <- as.numeric(-solve(tmphess, tmpgrad))
  print("hereese")
  
  while (max(abs(step)) > maxstep) step <- step / 2
  u <- u + step
  gradvals[itr] <- max(abs(tmpgrad))
  itr <- itr + 1
}
# Check convergence
itr # If itr == maxitr, algorithm was stopped prematurely
max(abs(nllg(u))) # If max(abs(tmpgrad)) > eps, algorithm failed to converge




# Confidence intervals
H <- nllh(u) # Sparse


L <- chol(H)  # Cholesky decomposition, L is lower triangular
# Solve L * y = e_i for each unit vector e_i to get columns of L^{-1}
inv_diag <- rep(NA, ncol(H))
for (i in 1:ncol(H)) {
  e_i <- rep(0, ncol(H))
  e_i[i] <- 1
  v <- solve(L, e_i, sparse = TRUE)  # Solving L * y = e_i
  # Since H^-1 = (L^T)^-1 * L^-1 and we need diag(H^-1), we use y directly
  inv_diag[i] <- sum(v^2)  # The i-th diagonal element of H^-1
}

# Standard deviations are the square roots of the diagonal entries
usd <- sqrt(inv_diag)





# usd <- sqrt(diag(solve(H))) # solve(H) is nxn but NOT SPARSE. Bad!
# Make the plot
pdf("stat 440 final plot.pdf", width = 10, height = 8)
plot(1:n, y, 
     main = "Observed counts and true vs predicted means (my code)", 
     xlab = "Day", ylab = "Count/mean")
lines(1:n, exp(u + beta0), lwd = 2, col = 'red')
lines(1:n, exp(u + beta0 - 2 * usd), lwd = 2, col = 'red', lty = "dashed")
lines(1:n, exp(u + beta0 + 2 * usd), lwd = 2, col = 'red', lty = "dashed")
dev.off()



# finally, check if it is a sparse matrix

# List all objects in the environment
objects <- ls()

# Initialize a variable to track the presence of dense square matrices matching dimension n
dense_square_found <- FALSE


# Loop through each object and check if it is a matrix and if it is sparse
for (obj_name in objects) {
  obj <- get(obj_name)  # Get the object by name
  if (is(obj, "Matrix") || is.matrix(obj)) {  # Check if object is a matrix or a Matrix object
    if (inherits(obj, "CsparseMatrix") || inherits(obj, "TsparseMatrix")) {
      cat(obj_name, "is a sparse matrix.\n")
    } else {
      # Check if the matrix dimensions are both exactly n
      if (is.matrix(obj) && dim(obj)[1] == n && dim(obj)[2] == n) {
        cat(obj_name, "is a dense square matrix of dimension n x n. Dimensions:", n, "x", n, "\n")
        dense_square_found <- TRUE
      } else {
        cat(obj_name, "is not a sparse matrix, but its dimensions do not qualify it as dense (n x n). Dimensions:", dim(obj)[1], "x", dim(obj)[2], "\n")
      }
    }
  }
}

# Final message based on the presence of dense matrices of size n x n
if (dense_square_found) {
  cat("Frown :(. There are dense square matrices of size n x n in your environment.\n")
} else {
  cat("Celebrate! There are no dense square matrices of size n x n; all others are sparse or do not meet the dense criteria.\n")
}
