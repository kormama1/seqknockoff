#' Simulate Gaussian and binary covariate predictors
#'
#' simulate Gaussian predictors with mean zero
#' and covariance structure determined by "cov_type" argument. Then p_b
#' randomly selected columns are dichotomized using the sign function
#'
#' @param n number of observations (rows of X)
#' @param p total number of covariates (columns of X) both continuous and binary
#' @param p_b number of binary covariates (0 <= p_b <= p)
#' @param cov_type character string specifying the covariance function. Can be one of
#' "cov_diag" (independent columns), "cov_equi" (equi-correlated columns), or "cov_ar1" (ar1-correlated columns).
#' The columns are shuffled during simulation
#' @param rho correlation parameter; input to the cov_type function
#'
#' @return the simulated data.frame with n rows and p columns (p_b of which are binary and p-p_b of which are gaussian).
#' Each column is either of class "numeric" or "factor".
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' # all columns are continuous:
#' X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' round(cor(X), 2)
#'
#' # two of the six columns are dichotomized (and set to class factor):
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # The class of each column:
#' unlist(lapply(X, class))
generate_X <- function(n, p, p_b, cov_type, rho=0.5) {

  p_c <- p - p_b

  sigma_z <- do.call(get(cov_type), args=list(n=n, p=p, rho=rho))

  Z <- data.frame(mvtnorm::rmvnorm(n, rep(0, p), sigma_z, "chol"))
  X <- Z

  # Random indices for continuous and binary columns of X
  inds_c <- sample(1:p, size=p_c, replace=FALSE)
  inds_b <- setdiff(1:p, inds_c)

  # Deterministic indices for continuous and binary columns of Z
  seq.c <- 1:p_c
  seq.b <- setdiff(1:p, seq.c)

  # Data frame:
  X[, inds_c] <- identity(Z[, seq.c])
  X[, inds_b] <- lapply((1+sign(Z[, seq.b]))/2, as.factor)

  return(X)

}

# Covariance matrices scaled to be approximately 1/n on the diagonal

cov_diag <- function(n, p, rho=NULL) {
  # Diagonal covariance
  s <- diag(p) / n
  return(s)
}

cov_equi <- function(n, p, rho = 0.5) {
  # Equicorrelated covariance
  s <- (diag(1 - rho, p, p) + rho) / n
  return(s)
}

cov_ar1 <- function(n, p, rho = 0.5) {
  # AR(1) covariance
  s <- toeplitz(rho^(0:(p - 1))) / n
  return(s)
}


#' Simulate Gaussian response from a sparse regression model
#'
#' @param x matrix corresponding to the regression design matrix
#' @param p_nn number of non-null covariate predictors. The regression coefficients corresponding to
#' columns 1:p_nn of x will be non-zero, all other are set to zero.
#' @param a amplitude of non-null regression coefficients
#'
#' @return
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' # Simulate 4 Gaussian and 2 binary covariate predictors:
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # Obtain the corresponding model matrix:
#' x <- model.matrix(~., data=X)[,-1]
#'
#' # Scale the binary dummy-variables so column-wise variance of x is same:
#' which.factor <- as.numeric(which(sapply(X, is.factor)))
#' x[,which.factor] <- x[,which.factor]/(0.5*sqrt(nrow(x)))
#'
#' # Simulate response from model y = 2*x[,1] + 2*x[,2] + epsilon, where epsilon ~ N(0,1)
#' y <- generate_y(x, p_nn=2, a=2)
generate_y <- function(x, p_nn, a) {

  n <- nrow(x)
  p <- ncol(x)

  beta <- rep(c(a, 0), c(p_nn, p - p_nn))

  mu <- as.numeric(x %*% beta)

  y <- mu + rnorm(n)

  return(y)

}
