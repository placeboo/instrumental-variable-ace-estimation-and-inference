rm(list = ls())
source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R")
library(MCMCpack)
library(xtable)
#---------------------------------------------------#
# Example 1, no confounder between treatment and outcome
# True probability distribution
# z, x, y
# P(HE) = 0.4
# P(HU) = 0.15
# P(AR) = 0.25
# P(NR) = 0.2

# P(y(x0), y(x1), x(z0), x(z1), x(z2), x(z3)) = P(y(x0), y(x1)) * P(x(z0), x(z1), x(z2), x(z3))

# P(x(z0), x(z1), x(z2), x(z3)) ~ uniform(0,1)
# sum(P(x(z0), x(z1), x(z2), x(z3))) = 1
#---------------------------------------------------#
p_he = 0.4
p_hu = 0.15
p_ar = 0.25
p_nr = 0.2

# true ACE
ace_true = p_he - p_hu
set.seed(17)
x_z.vec = sample(1:100, 4 * 4)
p_x_z.vec = x_z.vec / sum(x_z.vec) # normalization

# Vector sequence rank as following:
# 
# # 1. x(z0)=0, x(z1)=0, x(z2)=0, x(z3)=0
# # 2. x(z0)=1, x(z1)=0, x(z2)=0, x(z3)=0
# # 3. x(z0)=0, x(z1)=1, x(z2)=0, x(z3)=0
# # 4. x(z0)=0, x(z1)=0, x(z2)=1, x(z3)=0
# # 5. x(z0)=0, x(z1)=0, x(z2)=0, x(z3)=1
# # 6. x(z0)=1, x(z1)=1, x(z2)=0, x(z3)=0
# # 7. x(z0)=1, x(z1)=0, x(z2)=1, x(z3)=0
# # 8. x(z0)=1, x(z1)=0, x(z2)=0, x(z3)=1
# # 9. x(z0)=0, x(z1)=1, x(z2)=1, x(z3)=0
# # 10. x(z0)=0, x(z1)=1, x(z2)=0, x(z3)=1
# # 11. x(z0)=0, x(z1)=0, x(z2)=1, x(z3)=1
# # 12. x(z0)=1, x(z1)=1, x(z2)=1, x(z3)=0
# # 13. x(z0)=1, x(z1)=1, x(z2)=0, x(z3)=1
# # 14. x(z0)=1, x(z1)=0, x(z2)=1, x(z3)=1
# # 15. x(z0)=0, x(z1)=1, x(z2)=1, x(z3)=1
# # 16. x(z0)=1, x(z1)=1, x(z2)=1, x(z3)=1

# Observed probability
p_observed = rep(0, 4 * 4)
# Z=0
# # P(Y=0, X=0|Z=0) = P(Y(x0)=0, X(z0)=0)
p_observed[1] = (p_he + p_nr) * sum(p_x_z.vec[c(1,3:5,9:11, 15)])
# # P(Y=1, X=0|Z=0) = P(Y(x0)=1, X(z0)=0)
p_observed[2] = (p_ar + p_hu) * sum(p_x_z.vec[c(1,3:5,9:11, 15)])
# # P(Y=0, X=1|Z=0) = P(Y(x1)=0, X(z0)=1)
p_observed[3] = (p_nr + p_hu) * sum(p_x_z.vec[c(2, 6:8, 12:14, 16)])
# # P(Y=1, X=1|Z=0) = P(Y(x1)=1, X(z0)=1)
p_observed[4] = (p_ar + p_he) * sum(p_x_z.vec[c(2, 6:8, 12:14, 16)])
#
# Z=1
# # P(Y=0, X=0|Z=1) = P(Y(x0)=0, X(z1)=0)
p_observed[5] = (p_he + p_nr) * sum(p_x_z.vec[c(1,2,4,5,7,8,11,14)])
# # P(Y=1, X=0|Z=1) = P(Y(x0)=1, X(z1)=0)
p_observed[6] = (p_ar + p_hu) * sum(p_x_z.vec[c(1,2,4,5,7,8,11,14)])
# # P(Y=0, X=1|Z=1) = P(Y(x1)=0, X(z1)=1)
p_observed[7] = (p_nr + p_hu) * sum(p_x_z.vec[-c(1,2,4,5,7,8,11,14)])
# # P(Y=1, X=1|Z=1) = P(Y(x1)=1, X(z1)=1)
p_observed[8] = (p_ar + p_he) * sum(p_x_z.vec[-c(1,2,4,5,7,8,11,14)])
#
# Z=2
# # P(Y=0, X=0|Z=2) = P(Y(x0)=0, X(z2)=0)
p_observed[9] = (p_he + p_nr) * sum(p_x_z.vec[c(1:3,5,6,8,10,13)])
# # P(Y=1, X=0|Z=2) = P(Y(x0)=1, X(z2)=0)
p_observed[10] = (p_ar + p_hu) * sum(p_x_z.vec[c(1:3,5,6,8,10,13)])
# # P(Y=0, X=1|Z=2) = P(Y(x1)=0, X(z2)=1)
p_observed[11] = (p_nr + p_hu) * sum(p_x_z.vec[-c(1:3,5,6,8,10,13)])
# # P(Y=1, X=1|Z=2) = P(Y(x1)=1, X(z2)=1)
p_observed[12] = (p_ar + p_he) * sum(p_x_z.vec[-c(1:3,5,6,8,10,13)])
#
# Z=3
# # P(Y=0, X=0|Z=3) = P(Y(x0)=0, X(z3)=0)
p_observed[13] = (p_he + p_nr) * sum(p_x_z.vec[c(1:4, 6,7, 9, 12)])
# # P(Y=1, X=0|Z=3) = P(Y(x0)=1, X(z3)=0)
p_observed[14] = (p_ar + p_hu) * sum(p_x_z.vec[c(1:4, 6,7, 9, 12)])
# # P(Y=0, X=1|Z=3) = P(Y(x1)=0, X(z3)=1)
p_observed[15] = (p_nr + p_hu) * sum(p_x_z.vec[-c(1:4, 6,7, 9, 12)])
# # P(Y=1, X=1|Z=3) = P(Y(x1)=1, X(z3)=1)
p_observed[16] = (p_ar + p_he) * sum(p_x_z.vec[-c(1:4, 6,7, 9, 12)])

names(p_observed) = paste(c("X0,Y0", "X0,Y1", "X1,Y0", "X1,Y1"), "|Z",rep(0:3, each =4), sep = "")

# True ACE bounds
ace_bound_true = ace.bounds(p_observed)

save(ace_true, ace_bound_true, p_observed, file = "bayesian-simulation/true_ace_and_ace-bound.Rdata")


# table
tmp = matrix(as.character(c(0,0,0,1,1,0,1,1)), ncol =2, byrow = T)
mat = cbind(rep(0:3, each = 4), rbind(tmp, tmp, tmp, tmp))
mat = cbind(mat, round(p_observed, 3))
print(xtable(mat), include.rownames=FALSE)

#---------------------------------------------------#
# example 2. more general
# P(y(x0), y(x1), x(z0), x(z1), x(z2), x(z3)) \neq P(y(x0), y(x1)) * P(x(z0), x(z1), x(z2), x(z3))
# there are 2^6 = 64
# order as following
# y(x0)=0, y(x1)=0, 
# # 1. x(z0)=0, x(z1)=0, x(z2)=0, x(z3)=0
# # 2. x(z0)=1, x(z1)=0, x(z2)=0, x(z3)=0
# # 3. x(z0)=0, x(z1)=1, x(z2)=0, x(z3)=0
# # 4. x(z0)=0, x(z1)=0, x(z2)=1, x(z3)=0
# # 5. x(z0)=0, x(z1)=0, x(z2)=0, x(z3)=1
# # 6. x(z0)=1, x(z1)=1, x(z2)=0, x(z3)=0
# # 7. x(z0)=1, x(z1)=0, x(z2)=1, x(z3)=0
# # 8. x(z0)=1, x(z1)=0, x(z2)=0, x(z3)=1
# # 9. x(z0)=0, x(z1)=1, x(z2)=1, x(z3)=0
# # 10. x(z0)=0, x(z1)=1, x(z2)=0, x(z3)=1
# # 11. x(z0)=0, x(z1)=0, x(z2)=1, x(z3)=1
# # 12. x(z0)=1, x(z1)=1, x(z2)=1, x(z3)=0
# # 13. x(z0)=1, x(z1)=1, x(z2)=0, x(z3)=1
# # 14. x(z0)=1, x(z1)=0, x(z2)=1, x(z3)=1
# # 15. x(z0)=0, x(z1)=1, x(z2)=1, x(z3)=1
# # 16. x(z0)=1, x(z1)=1, x(z2)=1, x(z3)=1
#
# y(x0)=0, y(x1)=1
# # above 16
#
# y(x0)=1, y(x1)=0
# # above 16
#
# y(x0)=1, y(x1)=1
# # above 16
#---------------------------------------------------#
set.seed(2)

alpha = sample(x = 0.1 * c(1:10),size = 64, replace = T)
# P(y(x0), y(x1), x(z0), x(z1), x(z2), x(z3))
p.vec = rdirichlet(1, alpha)
# each row is one patient category
p.mat = matrix(p.vec, nrow = 4, byrow = T)
rownames(p.mat) = c("NR", "HE", "HU", "AR")
p_nr = sum(p.mat[1,])
p_he = sum(p.mat[2,])
p_hu = sum(p.mat[3,])
p_ar = sum(p.mat[4,])

# true ACE
ace_true = p_he - p_hu

# Observed probability
p_observed = rep(0, 4 * 4)
# Z=0
# # P(Y=0, X=0|Z=0) = P(Y(x0)=0, X(z0)=0)
p_observed[1] = sum(p.mat[c(1,2), c(1,3:5,9:11, 15)])
# # P(Y=1, X=0|Z=0) = P(Y(x0)=1, X(z0)=0)
p_observed[2] = sum(p.mat[c(3,4), c(1,3:5,9:11, 15)])
# # P(Y=0, X=1|Z=0) = P(Y(x1)=0, X(z0)=1)
p_observed[3] = sum(p.mat[c(1,3), -c(1,3:5,9:11, 15)])
# # P(Y=1, X=1|Z=0) = P(Y(x1)=1, X(z0)=1)
p_observed[4] = sum(p.mat[c(2,4), -c(1,3:5,9:11, 15)])
#
# Z=1
# # P(Y=0, X=0|Z=1) = P(Y(x0)=0, X(z1)=0)
p_observed[5] = sum(p.mat[c(1,2), c(1,2,4,5,7,8,11,14)])
# # P(Y=1, X=0|Z=1) = P(Y(x0)=1, X(z1)=0)
p_observed[6] = sum(p.mat[c(3,4), c(1,2,4,5,7,8,11,14)])
# # P(Y=0, X=1|Z=1) = P(Y(x1)=0, X(z1)=1)
p_observed[7] = sum(p.mat[c(1,3), -c(1,2,4,5,7,8,11,14)])
# # P(Y=1, X=1|Z=1) = P(Y(x1)=1, X(z1)=1)
p_observed[8] = sum(p.mat[c(2,4), -c(1,2,4,5,7,8,11,14)])
#
# Z=2
# # P(Y=0, X=0|Z=2) = P(Y(x0)=0, X(z2)=0)
p_observed[9] = sum(p.mat[c(1,2), c(1:3,5,6,8,10,13)])
# # P(Y=1, X=0|Z=2) = P(Y(x0)=1, X(z2)=0)
p_observed[10] = sum(p.mat[c(3,4), c(1:3,5,6,8,10,13)])
# # P(Y=0, X=1|Z=2) = P(Y(x1)=0, X(z2)=1)
p_observed[11] = sum(p.mat[c(1,3), -c(1:3,5,6,8,10,13)])
# # P(Y=1, X=1|Z=2) = P(Y(x1)=1, X(z2)=1)
p_observed[12] = sum(p.mat[c(2,4), -c(1:3,5,6,8,10,13)])
#
# Z=3
# # P(Y=0, X=0|Z=3) = P(Y(x0)=0, X(z3)=0)
p_observed[13] = sum(p.mat[c(1,2), c(1:4, 6,7, 9, 12)])
# # P(Y=1, X=0|Z=3) = P(Y(x0)=1, X(z3)=0)
p_observed[14] = sum(p.mat[c(3,4), c(1:4, 6,7, 9, 12)])
# # P(Y=0, X=1|Z=3) = P(Y(x1)=0, X(z3)=1)
p_observed[15] = sum(p.mat[c(1,3), -c(1:4, 6,7, 9, 12)])
# # P(Y=1, X=1|Z=3) = P(Y(x1)=1, X(z3)=1)
p_observed[16] = sum(p.mat[c(2,4), -c(1:4, 6,7, 9, 12)])

names(p_observed) = paste(c("X0,Y0", "X0,Y1", "X1,Y0", "X1,Y1"), "|Z",rep(0:3, each =4), sep = "")

# True ACE bounds
ace_bound_true = ace.bounds(p_observed)

save(ace_true, ace_bound_true, p_observed, file = "bayesian-simulation/general_true_ace_and_ace-bound.Rdata")

# table
tmp = matrix(as.character(c(0,0,0,1,1,0,1,1)), ncol =2, byrow = T)
mat = cbind(rep(0:3, each = 4), rbind(tmp, tmp, tmp, tmp))
mat = cbind(mat, round(p_observed, 3))
print(xtable(mat), include.rownames=FALSE)


