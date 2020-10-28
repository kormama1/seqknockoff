## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE------------------------------------------------
library(seqknockoff)

## ---- eval=FALSE---------------------------------------------------------
#  # Generate a 2000 x 30 Gaussian data.frame under AR1(rho=0.5) covariance structure,
#  # with 10 of the columns  dichotomized
#  X <- generate_X(n=2000, p=30, p_b=10, cov_type = "cov_ar1", rho=0.5)

## ---- eval=FALSE---------------------------------------------------------
#  # Generate y ~ N(X%*%beta,I_n) where first 10 beta-coefficients are = 3.5, all other = 0.
#  y <- generate_y(X = X, p_nn=10, a=3.5)

## ---- warning=FALSE, eval=FALSE------------------------------------------
#  # Simulate sequential knockoff of X:
#  Xk <- knockoffs_seq(X)

## ---- eval=FALSE---------------------------------------------------------
#  # Simulate sequential knockoff of X based on sequential elastic-net with penalty parameter alpha.
#  Xk <- knockoffs_seq(X, alpha=0.5)

## ---- eval=FALSE---------------------------------------------------------
#  # Generate a 2000 x 30 Gaussian data.frame under AR1(rho=0.5) covariance structure,
#  X <- generate_X(n=2000, p=30, p_b=0, cov_type = "cov_ar1", rho=0.5)
#  # Simulate second order multivariate Gaussian MX knockoff:
#  Xk <- knockoffs_mx(X)

## ---- message=FALSE, warning=FALSE---------------------------------------
set.seed(123)
X <- generate_X(n=2000, p=30, p_b=10, cov_type = "cov_ar1", rho=0.5)
y <- generate_y(X = X, p_nn=10, a=3.5)
S <- knockoff_filter(X, y, fdr=c(0.05, 0.1, 0.2))
S

## ---- message=FALSE------------------------------------------------------
library(clustermq)

## ---- message=FALSE------------------------------------------------------
set.seed(123)
X <- generate_X(n=2000, p=30, p_b=0, cov_type = "cov_ar1", rho=0.5)
y <- generate_y(X = X, p_nn=10, a=3.5)
S <- Q(knockoff_filter, fdr = rep(0.2, 50), const = list(X=X, y=y, knockoffs=knockoffs_mx, 
       statistic=stat_glmnet), pkgs="seqknockoff", n_jobs=50, seed=1)

## ---- eval=FALSE---------------------------------------------------------
#  S <- list()
#  for (count in 1:50) {
#    S[[count]] <- knockoff_filter(X, y, fdr=0.2, knockoffs=knockoffs_mx, statistic=stat_glmnet)
#  }

## ------------------------------------------------------------------------
# Obtain a single set of selected variables from the list of knockoff selections:
multi_select(S, p=30) # p is the total number of variables

## ---- fig.width=7.15, fig.height=4.5, fig.fullwidth=TRUE-----------------
plot_heatmap(S, labels=paste0("x",1:30))

## ------------------------------------------------------------------------
sim_experiment <- function(n, p, cov_type="cov_ar1", rho, p_nn, a) {
  # Generate Gaussian design matrix:
  X <- generate_X(n=n, p=p, p_b=0, cov_type="cov_ar1", rho=rho)
  # Simulate y given X:
  y <- generate_y(X=X, p_nn=p_nn, a=a)
  # Apply MX and sequential knockoff filters:
  S.mx <- knockoff_filter(X=X, y=y, fdr=0.2, knockoffs=knockoffs_mx)
  S.seq <- knockoff_filter(X=X, y=y, fdr=0.2, knockoffs=knockoffs_seq)
  # Evaluate fdp and tpp for each method:
  results <- data.frame(method = c("MX Knockoff", "Seq Knockoff"),
                        fdp = c(eval_fdp(S.mx, negatives=(p_nn+1):p), 
                                eval_fdp(S.seq, negatives=(p_nn+1):p)),
                        tpp = c(eval_tpp(S.mx, positives=1:p_nn), 
                                eval_tpp(S.seq, positives=1:p_nn)),
                        a = a)
  
  return(results)
  
}

## ---- message=FALSE, warning=FALSE---------------------------------------
# Recommended to run only on HPC:
results_list <- Q(sim_experiment, a = rep(1:5, each=100), const = list(n=2000, p=200, 
                  cov_type="cov_ar1", rho=0.5, p_nn=50), pkgs="seqknockoff", n_jobs=500, seed=1)
results <- do.call("rbind", results_list)

## ---- fig.width=7.15, fig.height=4.5, fig.fullwidth=TRUE-----------------
library(ggplot2)
ggplot(mapping=aes(x=factor(a), y=fdp), data=results) + 
  geom_boxplot(aes(fill=method)) + 
  stat_summary(aes(group=method, colour=method), fun.y="mean", geom="line") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), minor_breaks = NULL) +
  xlab("Amplitude") +
  ggtitle("False Discovery Rate Curves") +
  theme(legend.position="bottom")

## ---- fig.width=7.15, fig.height=4.5, fig.fullwidth=TRUE-----------------
ggplot(mapping=aes(x=factor(a), y=tpp), data=results) + 
  geom_boxplot(aes(fill=method)) + 
  stat_summary(aes(group=method, colour=method), fun.y="mean", geom="line") +
  xlab("Amplitude") +
  ggtitle("Power Curves") +
  theme(legend.position="bottom")

