#'@title Least squares evaluation
#'@description Computes the least squares cost of a data matrix at a changepoint
#'or for all values of the changepoint, given matrices of parameters.
#'@details If \code{tau="all"} the evaluation happens for all values of the
#'changepoint and a seqeunce of real values is returned. Otherwise if \code{tau}
#'is in (0,1) a single value at that point is returned.
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@param x an n by p matrix of observations.
#'@param tau a real parameter between (0,1) for point evaluation or 'all'.
#'@param mat_pre an p by p-1 of coefficients for the first segment
#'@param mat_pre an p by p-1 of coefficients for the second segment
#'@return
#'\item{}{A real value if \code{tau} is in (0,1)}
#'\item{}{A vector of reals if \code{tau="all"}}
#'@examples
#'p = 50; n = 225; t0 = 0.3
#'params = gen_param(p = p, reg_param = T)
#'sig1 = params$pre_cov; sig2 = params$post_cov
#'x = gen_data(n = n, p = p, t0 = t0, sig1, sig2)$data
#'pre_mat = params$mu_mat
#'post_mat = params$gamma_mat
#' # Evaluation at a point
#'ls_eval(x, 0.5, type = "graph", pre_mat, post_mat)
#' # Evaluation everywhere
#'ls_eval(x, "all",  type = "graph", pre_mat, post_mat)
ls_eval = function(x, tau, type = "mean", param_pre, param_post = NULL) {
  n = dim(x)[1]
  if ( type == "mean" ) {
    if (tau == "all"){
      mat1 = matrix(0, n, n)
      mat1[upper.tri(mat1,diag=T)] = 1
      mat2 = 1-mat1
      loss_mat = rowSums((x - param_pre) ^ 2) * mat1 +
        rowSums((x - param_post) ^ 2) * mat2
      return(( colSums(loss_mat) ) / n)
    } else if (tau>=0 & tau<1) {
        cut = floor(n * tau)
        pre = (1:cut); post = ((cut+1):n)
        dat_pre = x[pre, ]
        dat_post = x[post, ]
        return((1 / n) * ( sum((dat_pre - param_pre) ^ 2) +
                             sum((dat_post - param_post) ^ 2)))
    } else if (tau == 1) {
        return((1 / n) * ( sum((x - param_pre) ^ 2)))
    } else {
        stop("tau out of bounds")
    }
  }
  if ( type == "graph") {
    if (tau == "all"){
      mat1 = matrix(0, n, n)
      mat1[upper.tri(mat1, diag = T)] = 1
      mat2 = 1-mat1
      loss_mat = rowSums((x - pred_mat(x, param_pre)) ^ 2) * mat1 +
                 rowSums((x - pred_mat(x, param_post)) ^ 2) * mat2
      return(( colSums(loss_mat) ) / n)
    } else if (tau>=0 & tau<1) {
        cut = floor(n*tau)
        pre = (1:cut); post = ((cut+1):n)
        dat_pre = x[pre, ]
        dat_post = x[post, ]
        pred_pre = pred_mat(dat_pre, param_pre)
        pred_post = pred_mat(dat_post, param_post)
        error_pre = dat_pre - pred_pre
        error_post = dat_post - pred_post
        ( sum(error_pre ^ 2) + sum(error_post ^ 2) ) / n
    } else if (tau == 1) {
        pred = pred_mat(x, param_pre)
        return( (1 / n) * ( sum((x - pred) ^ 2)) )
    } else {
        stop("tau out of bounds")
    }
  }
}
