source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R")
library(MCMCpack)

load(file = "general_true_ace_and_ace-bound.Rdata")

do.one.tranp.v2 = function(){
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
        
        N_pair = nrow(ace.bound.mat)
        
        L.vec = ace.bound.mat[,1]
        L_sort.vec = sort(L.vec)
        U.vec = ace.bound.mat[,2]
        U_sort.vec = sort(U.vec)
        
        t.vec = sort(c(L.vec, U.vec)) #increasing order
        
        start = L_sort.vec[floor(0.95 * N_pair)]
        end = U_sort.vec[ceiling(0.05 * N_pair)]
        
        start_index = which(t.vec == start)[1]
        end_index = which(t.vec == end)[1]
        
        if(start_index > end_index){
                return(c("NA", "NA"))
        }else{
                T.vec = c()
                for(i in start_index:end_index){
                        t_i = t.vec[i]
                        if(sum(L.vec <= t_i) - sum(U.vec <= t_i) == ceiling(0.95 * N_pair)){
                                T.vec = c(T.vec, t_i)  
                        }
                }
                if(length(T.vec) == 0){
                        return(c("NA", "NA")) 
                }else{
                        return(c(min(T.vec), max(T.vec)))
                }
        }
}

set.seed(17)
bound.mat = replicate(1000, do.one.tranp.v2())
rownames(bound.mat) = c("LowerBounds", "UpperBounds")
save(bound.mat, file = "transparent95_bounds_1004.Rdata")
