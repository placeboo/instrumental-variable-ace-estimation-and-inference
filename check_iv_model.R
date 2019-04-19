check.iv.ineqs = function(pyx.z,verbose=TRUE){
      if(verbose){
            print(paste("p(y0, x0 | z0) + p(y1, x0 | z1) = ", round(pyx.z[1] + pyx.z[7],4)),quote=FALSE)
            print(paste("p(y1, x0 | z0) + p(y0, x0 | z1) = ", round(pyx.z[3] + pyx.z[5],4)),quote=FALSE)
            print(paste("p(y0, x1 | z0) + p(y1, x1 | z1) = ", round(pyx.z[2] + pyx.z[8],4)),quote=FALSE)
            print(paste("p(y1, x1 | z0) + p(y0, x1 | z1) = ", round(pyx.z[4] + pyx.z[6],4)),quote=FALSE)}
      violated <- (((pyx.z[1] + pyx.z[7])>1) | #4
                         ((pyx.z[3] + pyx.z[5])>1) | #4
                         ((pyx.z[2] + pyx.z[8])>1) | #5
                         ((pyx.z[4] + pyx.z[6])>1)   #5
      )
      if(verbose){
            if(violated){
                  print("pyx.z does not satisfy inequalities for IV Model!")}
            else{print("All IV inequalities satisfied: distribution compatible with IV model.")}}
      return((!violated))
}

## once it voilates the iv-model, we maximize the log-likehood subjecting to iv-inequalies to fix the empirical distributions
prob.iv = function(n){
      # p = (p(x0,y0|z0), p(x1,y0|z0), p(x0,y1|z0), p(x1,y1|z0), p(x0,y0|z1), p(x1,y0|z1), p(x0,y1|z1), p(x1,y1|z1))
      n0 = sum(n[1:4])
      n1 = sum(n[5:8])

      pyx.z = c(n[1:4]/n0, n[5:8]/n1)
      pyx.z_new = rep(0, 8)

      ## check whether it voilates the iv-inequalities
      is_voilate = c((pyx.z[1] + pyx.z[7])>1,
                           (pyx.z[3] + pyx.z[5])>1,
                           (pyx.z[2] + pyx.z[8])>1,
                           (pyx.z[4] + pyx.z[6])>1)
      if(any(is_voilate)){ # if it voilated the iv-inequalies
            pair_idx = which(is_voilate)
            ## index of voilating the iv-model
            index = rbind(c(1,3,2,4), c(7,5,8,6))[,pair_idx]
            idx_rest = rbind(c(1,3,2,4), c(7,5,8,6))[,-pair_idx]

            idx_z0 = index[1]
            idx_z1 = index[2]
            ## with respect to counts
            n0_voilate = n[idx_z0]
            n1_voilate = n[idx_z1]
            ## rest counts
            n0_rest = n[-idx_z0][1:3]
            n1_rest = n[-idx_z1][5:7]

            ## suppose P(y,x|z=0) + p(1-y,x|z=1)>1
            ## feasible set for P_new(y,x|z=0)
            lower.b = max(1 - n1_voilate/n1, 0)
            upper.b = min(n0_voilate/n0, 1)
            delta = upper.b - lower.b

            ## the maximizer of log likelihood L is q(y,x|z)
            ## q(y,x|z) = (n(xyz) + n(y(1-x)(1-z)) + n(yx(1-z)) + n((1-y)(1-x)(1-z)))/n

            ## check is q(y,x|z) is in the feasible set
            q = (n[idx_z0] + n1 - n[idx_z1])/sum(n)

            if(q > lower.b & q < upper.b){ # in the feasible set
                  pyx.z_new[idx_z0] = q
            }else if(q <= lower.b){
                  pyx.z_new[idx_z0] = lower.b + delta/10^4
            }else if(q >= upper.b){
                  pyx.z_new[idx_z0] = upper.b - delta/10^4
            }else{
                  stop("upper bound is smaller than lower bound")
            }

            ## P_new(1-y,x|z=1) = 1 - P_new(y,x|z=0)
            pyx.z_new[idx_z1] = 1 - pyx.z_new[idx_z0]


            for(i in idx_rest[1,]){
                  if(i == 1){ ## p(x0,y0|z0)
                        pyx.z_new[1] = n[1]/(n0 - n0_voilate) * pyx.z_new[idx_z1]
                        pyx.z_new[7] = n[7]/(n1 - n1_voilate) * pyx.z_new[idx_z0]
                  }else if(i == 2){ ## p(x1,y0|z0)
                        pyx.z_new[2] = n[2]/(n0 - n0_voilate) * pyx.z_new[idx_z1]
                        pyx.z_new[8] = n[8]/(n1 - n1_voilate) * pyx.z_new[idx_z0]
                  }else if(i == 3){ ## p(x0,y1|z0)
                        pyx.z_new[3] = n[3]/(n0 - n0_voilate) * pyx.z_new[idx_z1]
                        pyx.z_new[5] = n[5]/(n1 - n1_voilate) * pyx.z_new[idx_z0]
                  }else{ ## p(x1,y1|z1)
                        pyx.z_new[4] = n[4]/(n0 - n0_voilate) * pyx.z_new[idx_z1]
                        pyx.z_new[6] = n[6]/(n1 - n1_voilate) * pyx.z_new[idx_z0]
                  }
            }
      }else{
            pyx.z_new = pyx.z
      }
      return(pyx.z_new)
}


Lr = function(alternative_p, null_p, n){
      # p, n, are vectors: (x0,y0|z0), (x1,y0|z0), (x0,y1|z0), (x1,y1|z0),
      # (x0,y0|z1), (x1,y0|z1), (x0,y1|z1), (x1,y1|z1)
      # (p1): alternative model
      # (p2): null model
      lnl1 = sum(n * log(alternative_p))
      lnl2 = sum(n * log(null_p))

      return(2 * (lnl1 - lnl2))
}
