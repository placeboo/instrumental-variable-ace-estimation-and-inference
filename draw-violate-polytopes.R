rm(list = ls())
library(dplyr)
library(ggplot2)
library(rgl)
library(hexbin)
source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")
source("check-iv-z-arbitrary.R")

#---------------------------------------------------#
# Bayesion Method of IV model, arbitray level of z
#
# Notation:
# qij = P(Y=i, D=j|Z=z), p_ij = P(Y0=i, Y1=j), i,j = 0,1
# IV Inequality:
#     q11.bar + q01.bar <= 1;
#     q10.bar + q00.bar <= 1;
#     q10.bar <= q1.lowbar + q10-01.lowbar;
#     q00.bar <= q0.lowbar + q00-11.lowbar;
#     q11.bar <= q1.lowbar + q00-11.lowbar;
#     q01.bar <= q0.lowbar + q10-01.lowbar;
#     1 <= q1.lowbar + q10-01.lowbar + q0.lowbar + q00-11.lowbar
#---------------------------------------------------#

#---------------------------------------------------#
# Simulation, RCT
# generate data, with certain level of Z
# Z ~ Uniform{0, 1, 2, ..., |z|}
# V ~ Unifrom{0,1,2} covariates
# Y_{Z = i, Z=z} ~ Ber(1/2), i=0,1
# Xz|Yxz = y, V=v ~ Ber(expit(-4+y+z+5/2*z-3/2*v^2+y*(z+v)+v*z*(1+2*y)))
#---------------------------------------------------#
#n_sample = 100
#z_level = 3

n_sample = 100
z_level = 4

set.seed(1)
u = rnorm(n = n_sample, 0, 1)

z_true = sample(x = 0:(z_level-1), replace = T, size = n_sample)

#v = sample(x = 0:(z_level-1), replace = T, size = n_sample)
x_true = rbinom(n = n_sample,
                prob = exp(1 * u + 0.5 * z_true)/(1 + exp(1 * u + 0.5 * z_true)),
                size = 1)
y_true = rbinom(n = n_sample,
                prob = exp(3 * u + 0.7 * x_true)/(1 + exp(3 * u + 0.7 * x_true)),
                size = 1)


#y_miss1 = ifelse(exp(3 * u + 0.7 * (1-x_true))/(1 + exp(3 * u + 0.7 *(1- x_true)))>0.5,1,0)

y_miss = rbinom(n = n_sample,
                 prob = exp(3 * u + 0.7 *  (1-x_true))/(1 + exp(3 * u + 0.7 *  (1-x_true))),
                 size = 1)

ace_true = (sum(y_true[x_true == 1]) + sum(y_miss[x_true==0]))/n_sample - (sum(y_true[x_true == 0]) + sum(y_miss[x_true==1]))/n_sample


cor(z_true, y_true)
cor(z_true, x_true)

data_true = data.frame(z = z_true, x = x_true, y = y_true)

count_true = data_true %>% group_by(z, x, y) %>%
        summarise(n = n()) %>%
        as.data.frame()


arm_size = data_true %>% group_by(z) %>%
        summarise(n = n()) %>%
        as.data.frame()

data_true = merge(count_true, arm_size, by = "z")

#pyx.z_true = data_true[,4] / data_true[,5]

#(pyx.z_true)


# Specification of prior for p(x,y|z)
prior.z = rep(1, 4) # for each arm

## post
post.z = rep(prior.z, z_level) + count_true[, 4]

num.sims <- 10000  # Number of posterior simulations for each arm

theta.sims = t(replicate(num.sims, dirichlet(post.z[1:4])))
colnames(theta.sims) = c("X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0")
for(i in 1: (z_level-1)){
        theta.sim.z = t(replicate(num.sims, dirichlet(post.z[(i*4+1):((i+1)*4)])))
        name = paste(c("X0,Y0", "X0,Y1", "X1,Y0", "X1,Y1"), "|Z", i, sep = "")
        colnames(theta.sim.z) = name
        theta.sims = cbind(theta.sims,theta.sim.z)
}

summary(theta.sims)

## Remove sampled distributions that violate the IV inequalites:
all.ivs.ok <- rep(NA,num.sims)

for (i in 1:num.sims){
        all.ivs.ok[i] <- check.iv.ineq.arbz(theta.sims[i,])
}

summary(all.ivs.ok)

mean(all.ivs.ok) # 0.2902

posterior.theta.sims.iv <- theta.sims[all.ivs.ok,]

ace.bound.mat = matrix(ncol = 2, nrow = nrow(posterior.theta.sims.iv))
for(k in 1: nrow(posterior.theta.sims.iv)){
      ace.bound.mat[k,] = ace.bounds(posterior.theta.sims.iv[k,])
}
png("figure/z4-effec-tri-plot.png")
do.tri.plot(ace.bound.mat, title.txt = "ACE Bounds", line = ace_true)
dev.off()

h <- hexbin(ace.bound.mat)

png("figure/z4-effec-hexbin-plot.png")
plot(h, xlab = "Lower Bound on ACE", ylab = "Upper Bound on ACE")
dev.off()

# P(lower_bound < true ace, upper_bound > true ace)
N = nrow(posterior.theta.sims.iv)
length(intersect(c(1: N)[ace.bound.mat[,1] < ace_true], 
          c(1: N)[ace.bound.mat[,2] > ace_true])) / N 
# 0.62474

#---------------------------------------------------#
# voilates the IV inequ
#---------------------------------------------------#
tmp = theta.sims[!all.ivs.ok, ]
check.tmp = matrix(NA, ncol = 3, nrow = nrow(tmp))
for(i in 1: nrow(tmp)){
        check.tmp[i, ] = which.iv.ineq.violate(tmp[i,])
}
# typical violates
# (F, T, T)
# (T, F, T)
# (T, T, F)
is_FTT = apply(check.tmp, 1, function(x) all(x == c("FALSE", "TRUE", "TRUE")))
is_TFT = apply(check.tmp, 1, function(x) all(x == c("TRUE", "FALSE", "TRUE")))
is_TTF = apply(check.tmp, 1, function(x) all(x == c("TRUE", "TRUE", "FALSE")))


prob_FTT = tmp[is_FTT, ][1,]
prob_TFT = tmp[is_TFT, ][1,]
prob_TTF = tmp[is_TTF, ][1,]


# output
tmp = rbind(c(0,0,0), c(0,0,1), c(0,1,0),c(0,1,1))
tmp = apply(rbind(tmp, tmp, tmp, tmp),2, as.character)
out.mat = cbind(tmp, round(prob_FTT, 3), round(prob_TFT,3), round(prob_TTF,3))

print(xtable(out.mat), include.rownames = F)
#---------------------------------------------------#
# polytopes
#---------------------------------------------------#
switch.fnt = function(prob){
        # original "X0,Y0|Z0","X0,Y1|Z0","X1,Y0|Z0","X1,Y1|Z0"
        # switch to
        # "X0,Y0|Z0","X1,Y0|Z0","X0,Y1|Z0","X1,Y1|Z0"
        tmp = prob[2]
        prob[2] = prob[3]
        prob[3] = tmp
        return(prob)
}

# FTT
png("figure/z4-violate-two-level.png")
simp <- do.simplex(phi=30,theta=120,r=1000,main="IV Inequalities, only violate the ineq involved 2 levels of IV.") # Set up plot
do.polytope(switch.fnt(prob_FTT[1:4]), simp, red="green")
do.polytope(switch.fnt(prob_FTT[5:8]), simp, red="blue")
do.polytope(switch.fnt(prob_FTT[9:12]), simp, red="purple")
do.polytope(switch.fnt(prob_FTT[13:16]), simp, red="orange")
dev.off()

png("figure/z4-intersection-violate-two-level.png")
simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersections, only violate the ineq involved 2 levels of IV.") # Set up plot
do.polytope(c(switch.fnt(prob_FTT[1:4]), switch.fnt(prob_FTT[5:8])), simp, red="green")
do.polytope(c(switch.fnt(prob_FTT[1:4]),switch.fnt(prob_FTT[9:12])), simp, red="purple")
do.polytope(c(switch.fnt(prob_FTT[1:4]),switch.fnt(prob_FTT[13:16])), simp, red="orange")
do.polytope(c(switch.fnt(prob_FTT[13:16]), switch.fnt(prob_FTT[5:8])), simp, red="blue")
do.polytope(c(switch.fnt(prob_FTT[13:16]),switch.fnt(prob_FTT[9:12])), simp, red="cyan")
#do.polytope(c(switch.fnt(prob_FTT[5:8]),switch.fnt(prob_FTT[9:12])), simp, red="magenta") # no interaction
dev.off()

# png("figure/violate-two-level-z0z1-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A0 and A1") # Set up plot
# do.polytope(c(switch.fnt(prob_FTT[1:4]), switch.fnt(prob_FTT[5:8])), simp, red="green")
# dev.off()
#
#
# png("figure/violate-two-level-z0z2-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A0 and A2") # Set up plot
# do.polytope(c(switch.fnt(prob_FTT[1:4]),switch.fnt(prob_FTT[9:12])), simp, red="purple")
# dev.off()

## Set up interactive plot
open3d()
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1))
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)
do.polytope.rgl(switch.fnt(prob_FTT[1:4]), red="green")
do.polytope.rgl(switch.fnt(prob_FTT[5:8]), red="blue")
do.polytope.rgl(switch.fnt(prob_FTT[9:12]), red="purple")
do.polytope.rgl(switch.fnt(prob_FTT[13:16]), red="orange")
browseURL(paste("file://", writeWebGL(dir=file.path("figure/", "webGL"), width=700), sep=""))

open3d()
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1))
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)
do.polytope.rgl(c(switch.fnt(prob_FTT[1:4]), switch.fnt(prob_FTT[5:8])), red="purple")
do.polytope.rgl(c(switch.fnt(prob_FTT[1:4]),switch.fnt(prob_FTT[9:12])), red="green")
do.polytope.rgl(c(switch.fnt(prob_FTT[1:4]),switch.fnt(prob_FTT[13:16])), red="orange")
do.polytope.rgl(c(switch.fnt(prob_FTT[13:16]), switch.fnt(prob_FTT[5:8])), red="blue")
do.polytope.rgl(c(switch.fnt(prob_FTT[13:16]),switch.fnt(prob_FTT[9:12])), red="cyan")
do.polytope.rgl(c(switch.fnt(prob_FTT[5:8]),switch.fnt(prob_FTT[9:12])), red="magenta")
browseURL(paste("file://", writeWebGL(dir=file.path("figure/", "webGL"), width=700), sep=""))

########## TFT
png("figure/z4-violate-three-level.png")
simp <- do.simplex(phi=30,theta=120,r=1000,main="IV Inequalities, only violate the ineq involved 3 levels of IV.") # Set up plot
do.polytope(switch.fnt(prob_TFT[1:4]), simp, red="green")
do.polytope(switch.fnt(prob_TFT[5:8]), simp, red="blue")
do.polytope(switch.fnt(prob_TFT[9:12]), simp, red="purple")
do.polytope(switch.fnt(prob_FTT[13:16]), simp, red="orange")
dev.off()

png("figure/z4-intersection-violate-three-level.png")
simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersections, only violate the ineq involved 3 levels of IV.") # Set up plot
do.polytope(c(switch.fnt(prob_TFT[1:4]), switch.fnt(prob_TFT[5:8])), simp, red="green")
do.polytope(c(switch.fnt(prob_TFT[1:4]),switch.fnt(prob_TFT[9:12])), simp, red="purple")
do.polytope(c(switch.fnt(prob_TFT[1:4]),switch.fnt(prob_TFT[13:16])), simp, red="orange")
do.polytope(c(switch.fnt(prob_TFT[13:16]), switch.fnt(prob_TFT[5:8])), simp, red="blue")
do.polytope(c(switch.fnt(prob_TFT[13:16]),switch.fnt(prob_TFT[9:12])), simp, red="cyan")
do.polytope(c(switch.fnt(prob_TFT[5:8]),switch.fnt(prob_TFT[9:12])), simp, red="magenta")
dev.off()

# png("figure/violate-three-level-z0z1-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A0 and A1") # Set up plot
# do.polytope(c(switch.fnt(prob_TFT[1:4]), switch.fnt(prob_TFT[5:8])), simp, red="purple")
# dev.off()
#
#
# png("figure/violate-three-level-z0z2-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A0 and A2") # Set up plot
# do.polytope(c(switch.fnt(prob_TFT[1:4]),switch.fnt(prob_TFT[9:12])), simp, red="green")
# dev.off()
#
# png("figure/violate-three-level-z1z2-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A1 and A2") # Set up plot
# do.polytope(c(switch.fnt(prob_TFT[5:8]),switch.fnt(prob_TFT[9:12])), simp, red="red")
# dev.off()

#
open3d()
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1))
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)
do.polytope.rgl(switch.fnt(prob_TFT[1:4]), red="green")
do.polytope.rgl(switch.fnt(prob_TFT[5:8]), red="blue")
do.polytope.rgl(switch.fnt(prob_TFT[9:12]), red="purple")
do.polytope.rgl(switch.fnt(prob_FTT[13:16]), red="orange")
browseURL(paste("file://", writeWebGL(dir=file.path("figure/", "webGL"), width=700), sep=""))

# collection of intersection
open3d()
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1))
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)
do.polytope.rgl(c(switch.fnt(prob_TFT[1:4]), switch.fnt(prob_TFT[5:8])), red="purple")
do.polytope.rgl(c(switch.fnt(prob_TFT[1:4]),switch.fnt(prob_TFT[9:12])), red="green")
do.polytope.rgl(c(switch.fnt(prob_TFT[1:4]),switch.fnt(prob_TFT[13:16])), red="orange")
do.polytope.rgl(c(switch.fnt(prob_TFT[13:16]), switch.fnt(prob_TFT[5:8])), red="blue")
do.polytope.rgl(c(switch.fnt(prob_TFT[13:16]),switch.fnt(prob_TFT[9:12])), red="cyan")
do.polytope.rgl(c(switch.fnt(prob_TFT[5:8]),switch.fnt(prob_TFT[9:12])), red="magenta")
browseURL(paste("file://", writeWebGL(dir=file.path("figure/", "webGL"), width=700), sep=""))

########## TTF
png("figure/z4-violate-4-level.png")
simp <- do.simplex(phi=30,theta=120,r=1000,main="IV Inequalities, only violate the ineq involved 3 levels of IV.") # Set up plot
do.polytope(switch.fnt(prob_TTF[1:4]), simp, red="green")
do.polytope(switch.fnt(prob_TTF[5:8]), simp, red="blue")
do.polytope(switch.fnt(prob_TTF[9:12]), simp, red="purple")
do.polytope(switch.fnt(prob_TTF[13:16]), simp, red="orange")
dev.off()

png("figure/z4-intersection-violate-4-level.png")
simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersections, only violate the ineq involved 3 levels of IV.") # Set up plot
do.polytope(c(switch.fnt(prob_TTF[1:4]), switch.fnt(prob_TTF[5:8]),switch.fnt(prob_TTF[9:12])), simp, red="purple")
do.polytope(c(switch.fnt(prob_TTF[1:4]), switch.fnt(prob_TTF[5:8]),switch.fnt(prob_TTF[13:16])), simp, red="green")
do.polytope(c(switch.fnt(prob_TTF[1:4]),switch.fnt(prob_TTF[9:12]), switch.fnt(prob_TTF[13:16])), simp, red="orange")
do.polytope(c(switch.fnt(prob_TTF[13:16]), switch.fnt(prob_TTF[9:12]), switch.fnt(prob_TTF[5:8])), simp, red="blue")
dev.off()

# png("figure/violate-three-level-z0z1-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A0 and A1") # Set up plot
# do.polytope(c(switch.fnt(prob_TFT[1:4]), switch.fnt(prob_TFT[5:8])), simp, red="purple")
# dev.off()
#
#
# png("figure/violate-three-level-z0z2-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A0 and A2") # Set up plot
# do.polytope(c(switch.fnt(prob_TFT[1:4]),switch.fnt(prob_TFT[9:12])), simp, red="green")
# dev.off()
#
# png("figure/violate-three-level-z1z2-intersectoin.png")
# simp <- do.simplex(phi=30,theta=120,r=1000,main="Intersection of A1 and A2") # Set up plot
# do.polytope(c(switch.fnt(prob_TFT[5:8]),switch.fnt(prob_TFT[9:12])), simp, red="red")
# dev.off()

#
open3d()
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1))
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)
do.polytope.rgl(switch.fnt(prob_TTF[1:4]), red="green")
do.polytope.rgl(switch.fnt(prob_TTF[5:8]), red="blue")
do.polytope.rgl(switch.fnt(prob_TTF[9:12]), red="purple")
do.polytope.rgl(switch.fnt(prob_TTF[13:16]), red="orange")
browseURL(paste("file://", writeWebGL(dir=file.path("figure/", "webGL"), width=700), sep=""))

# collection of intersection
open3d()
# resize window
par3d(windowRect = c(200, 200, 800, 800))
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1))
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)
do.polytope.rgl(c(switch.fnt(prob_TTF[1:4]), switch.fnt(prob_TTF[5:8]),switch.fnt(prob_TTF[9:12])), red="purple")
do.polytope.rgl(c(switch.fnt(prob_TTF[1:4]), switch.fnt(prob_TTF[5:8]),switch.fnt(prob_TTF[13:16])), red="green")
do.polytope.rgl(c(switch.fnt(prob_TTF[1:4]),switch.fnt(prob_TTF[9:12]), switch.fnt(prob_TTF[13:16])), red="orange")
do.polytope.rgl(c(switch.fnt(prob_TTF[13:16]), switch.fnt(prob_TTF[9:12]), switch.fnt(prob_TTF[5:8])), red="blue")
browseURL(paste("file://", writeWebGL(dir=file.path("figure/", "webGL"), width=700), sep=""))


print(xtable(prob_FTT))
