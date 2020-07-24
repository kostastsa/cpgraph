#'@title Split and fit function.
#'@description Split the data at a given point and fit to obtain parameters
#' for each segment.
#'@details This function separates the data into segments at a point
#' \code{n*t_split} and estimates parameters. If  \code{type="mean"} the
#' parameters are the soft-thresholded sample means of the segments.
#' If \code{type="graph"} the parameters are estimated using the procedure of
#' \insertCite{yuan2010high}{cpgraph}. In this case the parameters are obtained
#' using the \code{glmnet} regularization path with parameters specified in
#'\code{lambda_vec}. If \code{lambda_vec} is omitted the path is set to
#'\deqn{seq(lam_sep*(sqrt(log(max(p,n))/n),l_end,length.out=length_lam)}
#'If \code{par=TRUE} the coefficient estimation is done in parallel using
#'the library \code{parallel} with the number of CPU cores set to \code{NumCores}.
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@param x an n by p matrix of observations.
#'@param t_split a split value between 0 and 1.
#'@param type a string of characters denoting the detection type : "mean" or
#' "graph"
#'@param lambda_vec a grid of values of the regularization parameter
#' to be used for cross-validation.
#'@param length_lam length of grid.
#'@param lam_sep numeric value. Used in the regularizer grid if not provided by
#' the user.
#'@param par If TRUE parallelization of the lasso estimation is performed.
#'@param cv number of cross-validation folds.
#'@param bic_const penalty used in the BIC evaluation.
#'@param NumCores number of processor cores to be used
#'@param start_grid grip of values for the initial split location
#'@param l_end end value of the lambda sequence
#'@return
#'\item{param_pre}{Parameters for the first segment}
#'\item{param_post}{Parameters for the second segment}
#'@export
#' @references
#' \insertAllCited
#'@examples
#' p=50;n=225;t0=0.2
#' params=gen.param(p=p, reg_param = TRUE)
#' sig1=params$pre.cov; sig2=params$post.cov
#' x=gen.data(n=n,p=p,t0=t0,sig1,sig2)$data
#' # split in the middle of the dataset
#' snf = split_and_fit(x, t_split = 0.5)
#' pre = snf$pre_mat
#' post = snf$post_mat

split_and_fit = function(x, t_split, type = "mean", lambda_vec = NULL, cv = 5,
                         NumCores = 4, par = F, length_lam = 50, lam_sep = 1.75,
                         l_end = 2, k_seq = seq(0, 0.5, length.out = 50)) {
  n = dim(x)[1]
  p = dim(x)[2]
  cut = floor(n * t_split); pre = (1:cut); post = ((cut+1):n)
  k = bic_k(x, t_split = t_split)
  if ( type == "mean" ) {
    if ( t_split == 1 ) {
      param_pre = soft_thresh(mean(x), k)
      param_post = NULL
    } else {
      param_pre = soft_thresh(mean(x[pre, ]), k)
      param_post = soft_thresh(mean(x[post, ]), k)
    }
    return(list( param_pre = param_pre,
                 param_post = param_post ))
  }
  if ( type == "graph" ) {
    pvn = max(p, n)
    if (is.null(lambda_vec) == T) {
      gsp = lam_sep * ( sqrt(log(pvn) / n) )
      lambda_vec = seq(gsp, l_end, length.out = length_lam)
    }
    if ( t_split == 1 ) {
      param_pre = coef_estimate(x, lambda_vec, cv, NumCores, par,
                                length_lam, lam_sep, l_end)$coef_mat
      param_post = NULL
    } else {
      param_pre = coef_estimate(x[pre, ], lambda_vec, cv, NumCores, par,
                                length_lam, lam_sep, l_end)$coef_mat
      param_post = coef_estimate(x[post, ], lambda_vec, cv, NumCores, par,
                                 length_lam, lam_sep, l_end)$coef_mat
      }
    return(list( param_pre = param_pre,
                 param_post = param_post ))
  }
}
