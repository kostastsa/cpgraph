#' R6 class \code{cpgraph}
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@description
#' A cpgraph object is a data matrix along with information about the recursive
#' application of the detection method.
cpgraph <- R6Class("cpgraph",
                   public = list(
  #' @field data A n by p matrix of observations.
  data = "matrix",

  #' @field type The type of the change detection.
  type = "mean",

  #' @field run_flag A boolean flag denoting whether binary segmentation has been
  #' applied.
  run_flag = FALSE,

  #' @field cp_declared A vector of detected changepoints.
  cp_declared = NULL,

  #' @field cp_generation The generation of detected changepoints on the binary.
  #' tree
  cp_generation = NULL,

  #' @description
  #' Print method for cpgraph object.
  #' @param ... Object of class cpgraph.
    print = function(...) {
    shape = dim(self$data)
    cat("Object of class cpgraph with ", shape[1], " observations in ",
        shape[2], " dimensions.\nSearching for change of type '", self$type,
        "' ...\n", sep = "")
    if (self$run_flag == TRUE) {
      cat("Declared changepoints at locations: \n",
          sort(self$cp_declared))
      }
    },

  #'@description Finds changepoint estimates by recursively applying the split, fit
  #'and minimize procedure.
  #'
  #'The \code{bic_const} threshold parameter increases as recursion moves deeper
  #'at a rate proportional to n / (e - s) to avoid false positives due to small
  #'sample size.
  #'@param s Starting point of segment
  #'@param e Endpoint of segment
  #'@param edge Percentage of throaway data from the edges.
  #'@param lambda_vec A grid of values of the regularization parameter to be used
  #' for cross-validation.
  #'@param length_lam Length of grid.
  #'@param lam_sep Numeric value. Used in the regularizer grid if not provided by
  #' the user.
  #'@param par If TRUE parallelization of the lasso estimation is performed.
  #'@param cv Number of cross-validation folds.
  #'@param bic_const Penalty used in the BIC evaluation.
  #'@param NumCores Number of processor cores to be used
  #'@param start_grid Grid of values for the initial split location
  #'@param l_end End value of the lambda sequence
  #'@param gen A counter for the depth of recursion
  #'@export
  #'@examples
  #' p = 50; n1 = 500; n2 = 500; t01 = 0.2; t02 = 0.6
  #' params1 = gen_param(p = p)
  #' params2 = gen_param(p = p)
  #' sig1 = params1$pre_cov; sig2 = params1$post_cov
  #' sig3 = params2$pre_cov; sig4 = params2$post_cov
  #' x = gen_data(n = n1, p = p, t0 = t01, sig1, sig2)$data
  #' y = gen_data(n = n2, p = p, t0 = t02, sig3, sig4)$data
  #' data = rbind(x, y)
  #' dat = cpgraph$new()
  #' dat$data = data
  #' dat
  #' dat$bin_seg_std(1, n)
  #' dat
  changedetect = function(s = NULL, e = NULL, edge = 0.05, lambda_vec = NULL,
                          length_lam = 50, lam_sep = 1.75, par = F, cv = 5,
                          bic_const = 0.25, NumCores = 4, start_grid = 0.5,
                          l_end = 2, gen = 0) {
                self$run_flag <<- TRUE
                x = self$data
                n = dim(x)[1]
                if (is.null(s)) { s = 1 ; e = n }
                type = self$type
                gen = gen + 1
                if( e - s < 1) {}
                else {
                  est = double_update(x[s:e, ], type, edge, lambda_vec,
                                      length_lam, lam_sep, par, cv, bic_const,
                                      NumCores, start_grid, l_end)$upd2
                  if(est!=1) {
                    b_loc =  floor(est * (e - s)) + s
                    self$cp_declared <<- c(self$cp_declared, b_loc / n)
                    self$cp_generation <<- c(self$cp_generation, gen)
                    self$changedetect(s, b_loc, edge, lambda_vec,
                                     length_lam, lam_sep, par, cv,
                                     # bic_const increasing as segment size
                                     # decreases to avoid false positives
                                     bic_const = 0.25 * n / (b_loc - s),
                                     NumCores, start_grid, l_end, gen)
                    self$changedetect(b_loc + 1, e, edge, lambda_vec,
                                     length_lam, lam_sep, par, cv,
                                     bic_const = 0.25 * n / (e - b_loc),
                                     NumCores, start_grid, l_end, gen)
                  }
                }
              }
  ))

