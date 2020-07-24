soft_thresh = function(x, k){
  sign(x)*apply(as.array(x-k), 1, function(x) max(x,0))
}

