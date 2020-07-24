#'@title Detect possibly multiple changes of certain type.
#'@description Find changepoints using the split fit and minimize method in
#'a recursive binary segmentation algorithm.
#'@details The function creates a cpgraph object and runs the method change_detect
#'on the dataset. It returns the updated object.
#'@author Konstantinos Tsampourakis <kostas.tsampourakis@@gmail.com>
#'@param data an n by p matrix of observations.
#'@param type a string of characters denoting the detection type : "mean" or
#' "graph"
#'@return
#'\item{}{Returns cpgraph object with detected changepoints}
#'@export
#'@examples
#' p = 50; n = 225; t0 = 0.2
#' params = gen_param(p = p, reg_param = TRUE)
#' sig1 = params$pre_cov; sig2 = params$post_cov
#' x = gen_data(n = n, p = p, t0 = t0, sig1, sig2)$data
#' cp = changedetect(x, type = "graph")
#' cp
changedetect = function(data, type, verbose = FALSE) {
  cp_m = cpgraph$new()
  cp_m$data = data
  cp_m$type = type
  cp_m$change_detect()
  if ( verbose ) { cp_m }
  return(cp_m)
}
