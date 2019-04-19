rm(list = ls())

library(HDInterval)

priorz = 1
sample_size.vec = c(0, 100, 200, 400, 800, 1600)

# size_i = 200
# load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))
#
lambdas = seq(0, 1, 0.01)

getTheta_ci = function(L, U, lambda, ci = 0.95){
      # L: lower bounds for ACE
      # U: upper bounds for ACE

      # theta_lambda = lamdba * L + (1 - lambda) * U

      # ci: quantiles
      n = length(L)

      theta_lambda = lambda * L + (1 - lambda) * U

      # find the posterior density CI
      bounds = hdi(theta_lambda, ci)
      return(bounds)
}

ci.mat = matrix(NA, nrow = length(sample_size.vec), 2)

for (i in 1 : length(sample_size.vec)){

      size_i = sample_size.vec[i]
      load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))

      print(paste("sample size = ", size_i))

      tmp_ci.mat = matrix(NA, ncol = 2, nrow = length(lambdas))
      #tmp_ci.mat[, 1] = lambdas
      for (j in 1: length(lambdas)){
            lambda_j = lambdas[j]
            tmp_ci.mat[j, ] = getTheta_ci(ace.bound.mat[,1], ace.bound.mat[,2], lambda = lambda_j, ci = 0.95)
            print(paste("j=", j))
      }

      ci.mat[i,] = c(min(tmp_ci.mat[,1]), max(tmp_ci.mat[,2]))
}

save(ci.mat, file = "ace_point_coverage_2019-01-31.RData")
