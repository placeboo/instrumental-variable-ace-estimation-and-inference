#---------------------------------------------------#
# Notation:
# qij(z) = P(Y=i, X=j|Z=z), qi(z) = qi0(z) + qi1(z)
# q00-01 = q00(z) + q01(z)
# qij.bar = max_z(P(Y=i, X=j|Z=z))
# qij.lowbar = min_z(P(Y=i, X=j|Z=z))
# IV Inequality:
#     q11.bar + q01.bar <= 1;
#     q10.bar + q00.bar <= 1;
#     q10.bar <= q1.lowbar + q10-01.lowbar;
#     q00.bar <= q0.lowbar + q00-11.lowbar;
#     q11.bar <= q1.lowbar + q00-11.lowbar;
#     q01.bar <= q0.lowbar + q10-01.lowbar;
#     1 <= q1.lowbar + q10-01.lowbar + q0.lowbar + q00-11.lowbar
#---------------------------------------------------#
check.iv.ineq.arbz = function(pyx.z){
      # pyx.z = (P(x0,y0|z), P(x0,y1|z), P(x1,y0|z), P(x1,y1|z))
      zlevel = length(pyx.z)/4

      q00.vec = pyx.z[(0:(zlevel-1)) * 4 + 1]
      q10.vec = pyx.z[(0:(zlevel-1)) * 4 + 2]
      q01.vec = pyx.z[(0:(zlevel-1)) * 4 + 3]
      q11.vec = pyx.z[(1:zlevel) * 4]

      each.ineq = c(max(q11.vec) + max(q01.vec) <= 1,
        max(q10.vec) + max(q00.vec) <= 1,
        max(q10.vec) <= min(q10.vec + q11.vec) + min(q10.vec + q01.vec),
        max(q00.vec) <= min(q00.vec + q01.vec) + min(q00.vec + q11.vec),
        max(q11.vec) <= min(q10.vec + q11.vec) + min(q00.vec + q11.vec),
        max(q01.vec) <= min(q00.vec + q01.vec) + min(q10.vec + q01.vec),
        1 <= min(q10.vec + q11.vec) + min(q10.vec + q01.vec) + min(q00.vec + q01.vec) + min(q00.vec + q11.vec))

      return(all(each.ineq))
}

which.iv.ineq.violate = function(pyx.z){
        # pyx.z = (P(x0,y0|z), P(x0,y1|z), P(x1,y0|z), P(x1,y1|z))
        zlevel = length(pyx.z)/4

        q00.vec = pyx.z[(0:(zlevel-1)) * 4 + 1]
        q10.vec = pyx.z[(0:(zlevel-1)) * 4 + 2]
        q01.vec = pyx.z[(0:(zlevel-1)) * 4 + 3]
        q11.vec = pyx.z[(1:zlevel) * 4]

        each.ineq = c(max(q11.vec) + max(q01.vec) <= 1,
                      max(q10.vec) + max(q00.vec) <= 1,
                      max(q10.vec) <= min(q10.vec + q11.vec) + min(q10.vec + q01.vec),
                      max(q00.vec) <= min(q00.vec + q01.vec) + min(q00.vec + q11.vec),
                      max(q11.vec) <= min(q10.vec + q11.vec) + min(q00.vec + q11.vec),
                      max(q01.vec) <= min(q00.vec + q01.vec) + min(q10.vec + q01.vec),
                      1 <= min(q10.vec + q11.vec) + min(q10.vec + q01.vec) + min(q00.vec + q01.vec) + min(q00.vec + q11.vec))
        return(c(all(each.ineq[1:2]), all(each.ineq[3:6]), each.ineq[7]))
}


dirichlet <- function(alpha){
        theta <- vector(4,mode="numeric")
        for(i in 1:length(alpha)){
                if(!(alpha[i]==0))
                        theta[i] <- rgamma(1,alpha[i],1)
        }
        theta <- theta / sum(theta)
}

g_fct = function(pyx.z,i,j){
        zlevel = length(pyx.z)/4

        q00.vec = pyx.z[(0:(zlevel-1)) * 4 + 1]
        q10.vec = pyx.z[(0:(zlevel-1)) * 4 + 2]
        q01.vec = pyx.z[(0:(zlevel-1)) * 4 + 3]
        q11.vec = pyx.z[(1:zlevel) * 4]

        tmp1 = get(paste("q", j, i, ".vec", sep = ""))
        tmp2 = get(paste("q", "0", 1-i, ".vec", sep = ""))
        tmp3 = get(paste("q", "1", 1-i, ".vec", sep = ""))

        part1 = min(tmp1 + tmp2 + tmp3)

        part2 = rep(NA, zlevel)
        for(i in 1: zlevel){
                part2[i] = min(tmp1[i] + tmp2[i] + tmp1 + tmp3)
        }
        return(min(part1, min(part2)))
}


ace.bounds = function(pyx.z){
        ace.lower = 1 - g_fct(pyx.z, i = 0, j = 1) - g_fct(pyx.z, i = 1, j = 0)
        ace.upper = g_fct(pyx.z, i = 0, j = 0) + g_fct(pyx.z, i = 1, j = 1) - 1
        return(c(ace.lower, ace.upper))
}


#-------------------------------------------------#
# Thomas' codes
# small revise
# adding the location of red lines
#-------------------------------------------------#
# draw
do.tri.plot <- function(bounds,title.txt,col="#00000110", line = c(0,0)){
        plot(c(-1.001,1.001),c(-1.001,1.001)
             ,axes=FALSE,xlab=expression(paste("Lower Bound on ACE")),type="n",ylab=expression(paste("Upper Bound on ACE")),asp=1)
        axis(1,at=seq(-1,1,by=0.2),cex.axis=1)
        axis(2,at=seq(-1,1,by=0.2),cex.axis=1)
        title(title.txt,cex=3,line=1)
        points(c(-1,1,-1,-1),c(-1,1,1,-1), type="l",lty=1)
        points(bounds[,1],bounds[,2],pch=19, cex=0.3)
        #points(bounds[,1],bounds[,2],pch=".", cex=1)
        points(c(-1,line[1],line[1]),c(line[2],line[2],1),col="red",lty=1,type="l")
}


c.xyz <- function(x,y,z){
        return(1+ x + 2*y + 4*z)}

gfn <- function(pxy.z,x,y){
        num.z <- (length(pxy.z)/4)-1
        min.1 <- 1
        min.2 <- 1
        for(z in 0:num.z){
                min.1 <- min(pxy.z[c.xyz(x,y,z)]+pxy.z[c.xyz(1-x,y,z)]+pxy.z[c.xyz(1-x,1-y,z)],min.1)
        }
        if(num.z > 0){
                for(z in 1:num.z){
                        for(zs in 0:(z-1)){
                                min.2 <- min(pxy.z[c.xyz(x,y,z)]+pxy.z[c.xyz(1-x,0,z)]+pxy.z[c.xyz(x,y,zs)]+pxy.z[c.xyz(1-x,1,zs)],
                                             pxy.z[c.xyz(x,y,zs)]+pxy.z[c.xyz(1-x,0,zs)]+pxy.z[c.xyz(x,y,z)]+pxy.z[c.xyz(1-x,1,z)],
                                             min.2)
                        }
                }
        }
        return(min(min.1,min.2))
}


ace.bounds.tr <- function(pxy.z){
        lower.bnd <- 1 - gfn(pxy.z,1,0)- gfn(pxy.z,0,1)
        upper.bnd <- gfn(pxy.z,0,0) + gfn(pxy.z,1,1) - 1
        return(c(lower.bnd,upper.bnd))
}
