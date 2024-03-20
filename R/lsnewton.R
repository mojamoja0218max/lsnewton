#' least squares regression based on Newton's algorithm
#'
#' This package implements the Newton's algorithm for solving the least square regression
#'
#' @param X design matrix (n x p)
#' @param y outcome matrix (n x 1)
#' @param tolerance tolerance
#' @param max_iterations maximum number of iterations allowed
#' @return beta vector of length p
#' @export
#' @examples
#' # example code
#' set.seed(123)
#' n <- 100  # Number of data points
#' p <- 3    # Number of predictor variables
#' X <- matrix(rnorm(n * p), ncol = p)  # Create a design matrix with p columns
#' beta_true <- c(2, 3, 1)  # True coefficients
#' y <- X %*% beta_true + rnorm(n)  # Generate response variable
#' ls_newton(y, X)
ls_newton <- function(y, X, tolerance = 1e-06, max_iterations = 100) {
  ## assuming X contains the intercept
  p <- ncol(X)
  beta <- rep(0, p)
  for (iteration in 1:max_iterations) {
    r <-  y - X %*% beta
    step <- solve(t(X) %*% X) %*% t(X) %*% r
    beta <- beta + step
    if (max(abs(step)) < tolerance) {
      break
    }
  }
  ## Assign the lsnewton class and return
  class(beta) <- "lsnewton"

  ## Custom print function for lsnewton objects
  print.lsnewton <- function(obj){
    cat("Estimated beta coefficients:\n", obj)
  }

  print(beta)
}

