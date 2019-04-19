source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R")
library(MCMCpack)

load(file = "general_true_ace_and_ace-bound.Rdata")

do.one.tranp = function(){
      count.vec = c(rmultinom(1, 25, p_observed[1:4]),
                    rmultinom(1, 25, p_observed[5:8]),
                    rmultinom(1, 25, p_observed[9:12]),
                    rmultinom(1, 25, p_observed[13:16]))

      names(count.vec) = names(p_observed)

      # Specification of prior for p(x,y|z)
      prior.z = rep(1, 4) # for each arm

      ## post
      post.z = rep(prior.z, 4) + count.vec

      num.sims <- 10000  # Number of posterior simulations for each arm

      theta.sims = t(replicate(num.sims, dirichlet(post.z[1:4])))

      colnames(theta.sims) = c("X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0")
      for(i in 1: 3){
            theta.sim.z = t(replicate(num.sims, dirichlet(post.z[(i*4+1):((i+1)*4)])))
            name = paste(c("X0,Y0", "X0,Y1", "X1,Y0", "X1,Y1"), "|Z", i, sep = "")
            colnames(theta.sim.z) = name
            theta.sims = cbind(theta.sims,theta.sim.z)
      }


      ## Remove sampled distributions that violate the IV inequalites:
      all.ivs.ok <- rep(NA,num.sims)

      for (i in 1:num.sims){
            all.ivs.ok[i] <- check.iv.ineq.arbz(theta.sims[i,])
      }


      posterior.theta.sims.iv <- theta.sims[all.ivs.ok,]

      ace.bound.mat = matrix(ncol = 2, nrow = nrow(posterior.theta.sims.iv))

      for(k in 1: nrow(posterior.theta.sims.iv)){
            ace.bound.mat[k,] = ace.bounds(posterior.theta.sims.iv[k,])
      }

      # find the largest interval
      L_star = quantile(ace.bound.mat[,1], probs = 0.95)
      U_star = quantile(ace.bound.mat[,2], probs = 0.05)
      # print(c(L_star, U_star))
      # empircal CDF for lower and upper bounds
      lower_cdf = ecdf(ace.bound.mat[,1])
      upper_cdf = ecdf(ace.bound.mat[,2])
      if(L_star >= U_star){
            return(c("NA", "NA")) #There is no valid credible Interval
      }
      k = 0
      m = 0
      while (round(lower_cdf(L_star) - upper_cdf(L_star), 3) != 0.95){
            if(k == 51){
                  return(c("NA", "NA")) #There is no valid credible Interval
            }
            L_k = quantile(ace.bound.mat[,1], probs = 0.95 + 0.001 * k)
            L_star = L_k
            k = k + 1
      }
      while (round(lower_cdf(U_star) - upper_cdf(U_star), 3) != 0.95){
            if(m == 51){
                  return(c("NA", "NA")) #There is no valid credible Interval
            }
            U_m = quantile(ace.bound.mat[,2], probs = 0.05 - 0.001 * m)
            U_star = U_m
            m = m + 1
      }
      return(c(L_star, U_star))
}

set.seed(17)
bound.mat = replicate(1000, do.one.tranp())
rownames(bound.mat) = c("LowerBounds", "UpperBounds")
save(bound.mat, file = "transparent95_bounds.Rdata")
