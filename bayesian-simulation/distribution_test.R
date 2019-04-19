# c_h = (1-h)*c_0 + h*c_1
# c_0 <= c_1, random variable
library(MCMCpack)
# case 1, both skewed
nbins = 20
N = 1000

w0 = c(-.5, -.2)
w1 = c(-.4, .2)

w0_seq = seq(from = w0[1], to = w0[2], length.out =  nbins+1)
w1_seq = seq(from = w1[1], to = w1[2], length.out =  nbins+1)


seed = 2
set.seed(seed)

N0 = as.vector(rmultinom(1,N,rdirichlet(1, rep(1, nbins))))
N1 = as.vector(rmultinom(1,N,rdirichlet(1, rep(1, nbins))))

# N0 = round(N * rdirichlet(1, rep(1, nbins)), 0)
# N1 = round(N * rdirichlet(1, rep(1, nbins)), 0)

# sumN0 = sum(N0)
# sumN1 = sum(N1)
# if (sumN0 != sumN1) {
#       diff = abs(sumN0 - sumN1)
#       idx = sample(1:nbins, 1)
#       if (sumN0 < sumN1) {
#             N0[idx] = N0[idx] + diff
#       } else {
#             N1[idx] = N1[idx] + diff
#       }
# }
#N0 = sort(N0, decreasing = T)
#N1 = sort(N1, decreasing = F)

c0 = c()
c1 = c()
for (i in 1:nbins) {
      c0 = c(c0, runif(N0[i], min = w0_seq[i], max = w0_seq[i+1]))
      c1 = c(c1, runif(N1[i], min = w1_seq[i], max = w1_seq[i+1]))
}


c0 = sort(c0,decreasing=TRUE)
c1 = sort(c1)

# c0 < c1
idx = which(c0 > c1)
tmp = c0[idx]
c0[idx] = c1[idx]
c1[idx] = tmp

plot(density(c0))
lines(density(c1))

h = seq(0, 1, 0.01)

a = 0.05
l = quantile(c0, a)
u = quantile(c1, 1 - a)
# l = -0.3
# u = 0.4
print(paste("l=",l,"u=",u))
count.vec = c()
for (i in h) {
      ch = i * c0 + (1 - i) * c1
      count.vec = c(count.vec, sum(ch >= l & ch <= u) / length(ch))
}
plot(1:length(count.vec), count.vec)
min(count.vec)
