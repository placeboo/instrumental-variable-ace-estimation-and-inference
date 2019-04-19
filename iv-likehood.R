rm(list = ls())

source("check_iv_model.R")

## input counts
## p(x0,y0|z0), p(x1,y0|z0), p(x0,y1|z0), p(x1,y1|z0),
## p(x0,y0|z1), p(x1,y0|z1), p(x0,y1|z1), p(x1,y1|z1)


## fake data, which voilates the iv-model
count = c(150, 50, 100, 200, 50, 325, 100, 25)
p_emp = c(count[1:4]/sum(count[1:4]), count[5:8]/sum(count[5:8]))

check.iv.ineqs(p_emp)
## voilate


p_new = prob.iv(n = count)
check.iv.ineqs(p_new)

lr = Lr(alternative_p = p_emp, null_p = p_new, n = count)

## sample data from p_new, multinomial distribution
n0 = sum(count[1:4])
n1 = sum(count[5:8])
## for z=0, sample counts from p_new[1:4]
set.seed(17)
n_sample = 10000
sample_z0 = rmultinom(n_sample, size = n0, prob = p_new[1:4])
sample_z1 = rmultinom(n_sample, size = n1, prob = p_new[5:8])
sample_z = rbind(sample_z0, sample_z1)

## empirical probability for alternative
sample_p_z = rbind(sample_z0/n0, sample_z1/n1)
sample_pnew_z = apply(sample_z, 2, function(x) prob.iv(x))
# ## test
test = sample(1:n_sample, size = 1)
check.iv.ineqs(sample_pnew_z[, test])
sample_p_z[,test]
sample_pnew_z[,test]

## distribution of LR(sample, null)
sim_lr = 2 * (apply(sample_z * log(sample_p_z), 2, sum) - apply(sample_z * log(sample_pnew_z), 2, sum))

## Pr(sim_lr > lr)
sum(sim_lr >= lr) / n_sample

1-pchisq(lr, df = 1)
##
## plot
plot(1:n_sample, sim_lr)
hist(sim_lr, breaks = 50)
## 0.1168
