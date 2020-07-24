#'@title Coefficient estimation
#'@description Estimates elements of the inverse covariance matrix column-wise
#'by solving a sequence of lasso sub-problems for detection of \code{type="graph"}
#'. Each subproblem is solved using the library \code{glmnet}
#'\insertCite{friedman2010regularization}{cpgraph} in which the regularization
#'path is computed over a given range of values of the regularization parameter.
#'This function implements the 'fit' part of the split, fit and minimize procedure.
#'@details The parameters are obtained using the \code{glmnet} regularization
#'path with parameters specified in \code{lambda_vec}. If \code{lambda_vec} is
#'omitted the path is set to
#'\deqn{seq(lam.sep*(sqrt(log(max(p,n))/n),l.end,length.out=length.lam)}
#'If \code{par=TRUE} the coefficient estimation is done in parallel using the
#'library \code{parallel} with the number of CPU cores set to \code{NumCores}.
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@param x an n by p matrix of observations.
#'@param lambda_vec a grid of values of the regularization parameter to be used
#' for cross-validation.
#'@param cv number of cross-validation folds.
#'@param NumCores number of processor cores to be used
#'@param par If TRUE parallelization of the lasso estimation is performed.
#'@param length_lam length of grid.
#'@param lam_sep numeric value. Used in the regularizer grid if not provided by
#' the user.
#'@param l_end end value of the lambda sequence
#'@importFrom MASS mvrnorm
#'@return
#'\item{coef_mat}{Matrix of coefficients arranged columnwise}
#'\item{la_min_vec}{Vector of selected lambda values for each sub-problem.}
#'@export
#'@references
#'\insertAllCited
#'@examples
#' n = 300; p = 50
#' sigma = t(toeplitz(1:p)) %*% toeplitz(1:p)
#' x = mvrnorm(n, rep(0, p), sigma)
#' mod = coef_estimate(x)
#' mod$la_min_vec
#' mod$coef_mat[ , 1]

coef_estimate = function(x, lambda_vec = NULL, cv = 5, NumCores = 4, par = F,
                         length_lam = 50, lam_sep = 1.75, l_end = 2) {
  n = dim(x)[1]; p = dim(x)[2]; pvn = max(p, n)
  if (is.null(lambda_vec) == T) {
    gsp = lam_sep * ( sqrt(log(pvn)  /n) )
    lambda_vec = seq(gsp, l_end,length.out = length_lam)
    }
  coef_mat = matrix(0, (p-1),p);
  # For parallel run (par=T)
  if (par == T) {
    if (is.null(NumCores) == TRUE) { NumCores = max(1, detectCores()) }
    cl <- makeCluster(NumCores);  registerDoSNOW(cl)
    pb = txtProgressBar(max = p, style = 3);
    progress <- function(r) setTxtProgressBar(pb,r)
    opts <- list(progress = progress)
    f = foreach (j = 1:p, .packages = c("glmnet"),
                 .options.snow = opts) %dopar% {
    mod = cv.glmnet(x[ ,-j], x[ ,j],  nfolds = cv, lambda = lambda_vec,
                    intercept = FALSE)
    lam_min = mod$lambda.min;
    hb = coef(mod, s = lam_min)[-1]
    out = list(coef = hb, la_min = lam_min)
    }
    close(pb);  stopCluster(cl); la_min_vec=c()
    for (j in 1:p) {
      coef_mat[ , j] = f[[j]]$coef
      la_min_vec[j] = f[[j]]$la_min
      }
  }
  # For sequantial run (par=F)
  if (par==F) {
    la_min_vec = c();
    for(j in 1:p){mod = cv.glmnet(x[,-j], x[,j],  nfolds = cv,
                                  lambda = lambda_vec, intercept = FALSE)
    lam_min = mod$lambda.min;
    hb = coef(mod, s = lam_min)[-1]
    coef_mat[ ,j] = hb
    la_min_vec[j] = lam_min}
    }
  list( coef_mat = coef_mat,
        la_min_vec = la_min_vec)
  }



