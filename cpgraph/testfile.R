# # Graph dataset
# p=50; n1 = 500; n2 = 500; t01 = 0.2; t02 = 0.6
# params1 = gen_param(p = p, s1 = 2, s2 = 4)
# params2 = gen_param(p = p, s1 = 2, s2 = 3)
# sig1 = params1$pre_cov; sig2 = params1$post_cov
# sig3 = params2$pre_cov; sig4 = params2$post_cov
# x = gen_data(n = n1, p = p, t0 = t01, sig1, sig2)$data
# y = gen_data(n = n2, p = p, t0 = t02, sig3, sig4)$data
# data = rbind(x, y)



# Mean dataset
library(MASS)
x = mvrnorm(300, rep(1, p), sig1)
y = mvrnorm(700, rep(0, p), sig1)
data = rbind(x,y)

cp = changedetect(data, "mean", verbose = FALSE)
cp$cp_declared

p = 50; n = 225; t0 = 0.2
params = gen_param(p = p, reg_param = TRUE)
sig1 = params$pre_cov; sig2 = params$post_cov
x = gen_data(n = n, p = p, t0 = t0, sig1, sig2)$data
cp_list = changepoints(x, type = "graph")
cp_list



