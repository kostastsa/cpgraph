#'@title Prediction matrix calculation.
#'@description Computes the matrix of the predicted values according to a linear
#'model for all elements of the data matrix.
#'@param x an n by p matrix of observations.
#'@param coef_mat an p by p-1 of coefficients of the linear model
#' for each component.
#'@return
#'\item{pred}{A n by p matrix of the predicted values.}
#'@examples
#' n = 300; p = 100
#' x = matrix(rnorm(n * p), n, p)
#' mod = coef_estimate(x)
#' mod = pred_mat(x, coef.mat = mod$coef.mat)
#' pred = mod$predicted
pred_mat = function(x, coef_mat) {
  n = dim(x)[1]; p = dim(x)[2];
  pred = matrix(0, n, p)
  for (j in 1:p) {
    pred[ , j] = x[ , -j] %*% coef_mat[ , j]
    }
  pred
  }


