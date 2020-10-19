#' Sequential knockoffs for continuous and categorical variables
#'
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param seq_simulator function that simulates sequential knockoffs. Default is the function \code{sim_EN}, which simulates response from an estimated elastic-net model
#' @param ... other parameters passed to the function seq_simulator. For the default (elastic-net sequential simulator, \code{seq_simulator = sim_EN})
#' these other parameters are passed to cv.glmnet.
#'
#' @details \code{seqknockoff} performs sequential knockoff simulation.
#' @return sequential knockoff copy of X. A data.frame or tibble of same type and dimensions as X.
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # knockoffs based on sequential elastic-net regression with penalty alpha:
#' Xk <- seqknockoff(X, alpha=0.5)
seqknockoff <- function(X, seq_simulator = sim_EN, ...) {

  if (ncol(X)<=2) stop("The number of columns, ncol(X) needs to be > 2")

  knockoffs <- X

  # add ".tilde" to column names:
  names(knockoffs) <- paste0(names(knockoffs),".tilde")

  # Randomly shuffle column indicies of X:
  shf <- sample(ncol(X))

  # Loop through the columns of input data (in random order)
  loop.count <- 1
  for (i in shf) {

    y <- X[[i]] # i-th column serves as response
    x.mat <- X[,-i] # columns[-i] serve as covariates

    if (loop.count > 1) x.mat <- cbind(knockoffs[,shf[1:(loop.count-1)]], x.mat)

    x.mat <- model.matrix(~., data = x.mat)[,-1]

    knockoffs[[i]] <- seq_simulator(y = y, x = x.mat, ...)

    loop.count <- loop.count + 1

  }

  # remove ".tilde" from column names:
  names(knockoffs) <- gsub(".tilde","", names(knockoffs))

  return(knockoffs)

}


#' Simulate from elastic-net regression model
#'
#' @param y response vector (either "numeric" or "factor") that gets passed to cv.glmnet
#' @param x input matrix that gets passed to cv.glmnet
#' @param ... other parameters passed to the function cv.glmnet
#'
#' @return
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = x[,1] + rnorm(100)
#'
#' # simulate from elastic-net regression with elastic-net penalty alpha:
#' ysim = sim_EN(y=y, x=x, alpha=0.5)
#'
#' # simulated versus input response:
#' plot(y, ysim)
sim_EN <- function(y=y, x=x, ...) {

  if (is.factor(y)) {

    classes <- levels(y)

    K <- length(classes)

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="multinomial", intercept=TRUE, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.min")[[2]])[-1]

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.min")

    mat.multinom <- apply(mu, 1, function(prob) rmultinom(n=1, size=1, prob=prob))

    y.sim <- classes[apply((1:K)*mat.multinom, 2, max)]

    y.sim <- factor(y.sim, levels=classes)

    rmse <- NULL

  } else {

    if(!is.numeric(y)) stop("class(y) needs to be either 'numeric' or 'factor'")

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="gaussian", intercept=TRUE, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.min"))[-1]

    # columns of predictor matrix corresponding to non-zero beta.coefs:
    non.zero.cols <- which(beta.coefs != 0)

    # Total number of non-zero parameters (including intercept, hence + 1)
    s.lambda = length(non.zero.cols) + 1

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.min")

    rmse = sqrt(sum((y-mu)^2)/(length(y) - s.lambda))

    y.sim <- rnorm(n=length(y), mean=mu, sd=rmse)

  }

  return(y.sim)

}
