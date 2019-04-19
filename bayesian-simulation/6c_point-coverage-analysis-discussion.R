                                                                                                                                                               rm(list = ls())
load(file = "ace_point_coverage_2019-01-31.RData")
priorz = 1
sample_size.vec = c(0, 100, 200, 400, 800, 1600)

findEmpPro = function(Ls, Us, interval){
      tmp = 0.1 * Ls + 0.9 * Us
      count = sum((Ls >= interval[1] & Ls <= interval[2]) *
                        (Us >= interval[1] & Us <= interval[2]) *
                        (tmp >= interval[1] & tmp
                         <= interval[2]))
      return(count / length(Ls))
}

emp.vec = c()
for (i in 1 : length(sample_size.vec)){

        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))

        #emp.vec = c(emp.vec, findEmpPro(ace.bound.mat[,1], ace.bound.mat[,2], ci.mat[i,]))
        emp.vec = c(emp.vec, findEmpPro(ace.bound.mat[,1], ace.bound.mat[,2], quant.mat[i,]))
}

#----------------------------------------------------------#
# empirical point coverage
#----------------------------------------------------------#
lambdas = seq(0, 1, 0.05)
point_coverage = c()

for (i in 1 : length(sample_size.vec)){

        size_i = sample_size.vec[i]
        load(paste("bounds_sample_size_", size_i,"priorz_", priorz,".Rdata", sep = ""))

        Ls = ace.bound.mat[,1]
        Us = ace.bound.mat[,2]

        L_quant = quantile(Ls, 0.05)
        U_quant = quantile(Us, 0.95)

        count.vec = c()
        for (j in 1: length(lambdas)){
                lambda_j = lambdas[j]

                thetas = lambda_j * Ls + (1 - lambda_j) * Us

                count.vec = c(count.vec, sum(thetas >= L_quant & thetas <= U_quant) / length(thetas))
        }
        point_coverage = c(point_coverage, min(count.vec))

}
