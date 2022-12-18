## Statistics 240 midterm project sourcecode - Hongju Pae (11. 30. 2020)

# load packages --------------------------------------------------
library(MASS)
library(matrixcalc)
library(xtable)

# generate MV data (without function) -------------------------------------

# generate MV normal data
n = 100 # sample size, chosen by trial-and-error 
A_iev = mvrnorm(n, mu=mean_A, Sigma=cov_iev)
B_iev = mvrnorm(n, mu=mean_B, Sigma=cov_iev)
C_iev = mvrnorm(n, mu=mean_C, Sigma=cov_iev)
iev <- array(c(A_iev, B_iev, C_iev), dim=c(100, 3, 3))
dimnames(iev)[[3]] <- c("A", "B", "C")

A_iuv = mvrnorm(n, mu=mean_A, Sigma=cov_iuv)
B_iuv = mvrnorm(n, mu=mean_B, Sigma=cov_iuv)
C_iuv = mvrnorm(n, mu=mean_C, Sigma=cov_iuv)
iuv <- array(c(A_iuv, B_iuv, C_iuv), dim=c(100, 3, 3))
dimnames(iuv)[[3]] <- c("A", "B", "C")

A_pd = mvrnorm(n, mu=mean_A, Sigma=cov_pd)
B_pd = mvrnorm(n, mu=mean_B, Sigma=cov_pd)
C_pd = mvrnorm(n, mu=mean_C, Sigma=cov_pd)
pd <- array(c(A_pd, B_pd, C_pd), dim=c(100, 3, 3))
dimnames(pd)[[3]] <- c("A", "B", "C")

A_nd = mvrnorm(n, mu=mean_A, Sigma=cov_nd)
B_nd = mvrnorm(n, mu=mean_B, Sigma=cov_nd)
C_nd = mvrnorm(n, mu=mean_C, Sigma=cov_nd)
nd <- array(c(A_nd, B_nd, C_nd), dim=c(100, 3, 3))
dimnames(nd)[[3]] <- c("A", "B", "C")


# power testing function -----------------------------------------------------------

est.power=function(N, n, meanA, meanB, meanC, cov)
{
  sim.pvalues.pillai=rep(0,N)
  sim.pvalues.wilks=rep(0,N)
  sim.pvalues.hl=rep(0,N)
  sim.pvalues.roy=rep(0,N)
  
  for(i in 1:N){
    A = mvrnorm(n, mu=meanA, Sigma=cov)
    B = mvrnorm(n, mu=meanB, Sigma=cov)
    C = mvrnorm(n, mu=meanC, Sigma=cov)
    K <- array(c(A, B, C), dim=c(n, 3, 3))
    #dimnames(K)[[2]] <- c("mu1", "mu2", "mu3")
    dimnames(K)[[3]] <- c("A", "B", "C") # generate MVNormal data

    #rearrange the data such as the response matrix is an n-by-p matrix
    y <- cbind(
      mu1=c(K[,1,1],K[,1,2],K[,1,3]),
      mu2=c(K[,2,1],K[,2,2],K[,2,3]), 
      mu3=c(K[,3,1],K[,3,2],K[,3,3])
      )
    x <- rep(c("A", "B", "C"), each = n)
    # x = factor(c(rep("A", n), rep("B", n), rep("c", n)))
    obj = manova(y ~ x)

    # empirical power of manova ~ABC = p-value of each test statistic
    # manova() has default test as Pillai trace
    sim.pvalues.pillai[i] = summary(obj, test="Pillai")$stats[1,6]
    sim.pvalues.wilks[i] = summary(obj, test="Wilks")$stats[1,6]
    sim.pvalues.hl[i] = summary(obj, test="Hotelling-Lawley")$stats[1,6]
    sim.pvalues.roy[i] = summary(obj, test="Roy")$stats[1,6]
  }

  return(list(Pillai=mean(sim.pvalues.pillai<0.05), 
              Wilks=mean(sim.pvalues.wilks<0.05),
              HotellingLawley=mean(sim.pvalues.hl<0.05),
              Roy=mean(sim.pvalues.roy<0.05)
             )
        )
}

# set paramaters, trial and error -------------------------------------------------------
# number of sample population = 3 ... A, B, C
#p = 3 # number of features 
#N = 1000 # number of simulations 
#n = 100 # number of observations, n = n1 = n2 = n3 

# all mean and covariance values are set arbitrary (by trial-and-error)

cov_iev = matrix(1:9, nrow=3, ncol=3) # indp. and equal variance 
cov_iev[1,] <- c(1, 0, 0)
cov_iev[2,] <- c(0, 1, 0)
cov_iev[3,] <- c(0, 0, 1)
is.positive.definite(cov_iev) # posdef check

cov_iuv = matrix(1:9, nrow=3, ncol=3) # indp. and unequal variance 
cov_iuv[1,] <- c(1, 0, 0)
cov_iuv[2,] <- c(0, 0.5, 0)
cov_iuv[3,] <- c(0, 0, 0.1)
is.positive.definite(cov_iuv) # posdef check

cov_pd = matrix(1:9, nrow=3, ncol=3) # positively dependant  
cov_pd[1,] <- c(1, 0.5, 0.5)
cov_pd[2,] <- c(0.5, 1, 0.5)
cov_pd[3,] <- c(0.5, 0.5, 1)
is.positive.definite(cov_pd) # posdef check

cov_nd = matrix(1:9, nrow=3, ncol=3) # negatively dependant 
cov_nd[1,] <- c(1, -0.2, -0.2)
cov_nd[2,] <- c(-0.2, 1, -0.2)
cov_nd[3,] <- c(-0.2, -0.2, 1)
is.positive.definite(cov_nd) # posdef check

meanA = c(2.0,5.0,8.0) # mean vector for population A
meanB = c(2.1,5.1,8.1) # mean vector for population B
meanC = c(2.2,5.2,8.2) # mean vector for population C


#trial and error results 

# n = 50, 100, 200, 300, 500, 800 
# cov = iev, iuv, pd, nd 
# happens to have same power if the difference of mean is equal(as 0.0, 0.1, 0.2), cov controlled
# ... that means, same result if meanA=c(2.0,3.0,4.0)

est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_iev)
# 1, 2, 3 - 40%, 4 - 63% rejected 
est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_iev)
# 1, 2, 3 - 74%, 4 - 90% 
est.power(N=1000, n=300, meanA, meanB, meanC, cov=cov_iev)
# 1, 2, 3 - 90%, 4 - 97% 
est.power(N=1000, n=500, meanA, meanB, meanC, cov=cov_iev)
# 1, 2, 3 - 99.7%, 4 - 99.9% 

est.power(N=1000, n=50, meanA, meanB, meanC, cov=cov_iuv)
# 1, 2, 3 - 75%, 4 - 89% 
est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_iuv)
# 1, 2, 3 - 99.2%, 4 - 99% 
est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_iuv)
# 1, 2, 3 - 99%, 4 - 100% 

est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_pd)
# 1, 2, 3 - 18%, 4 - 40% 
est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_pd)
# 1, 2, 3 - 40%, 4 - 63% 
est.power(N=1000, n=300, meanA, meanB, meanC, cov=cov_pd)
# 1, 2, 3 - 57%, 4 - 77% 
est.power(N=1000, n=500, meanA, meanB, meanC, cov=cov_pd)
# 1, 2, 3 - 84%, 4 - 94% 
est.power(N=1000, n=800, meanA, meanB, meanC, cov=cov_pd)
# 1, 2, 3 - 97.6%, 4 - 99%
est.power(N=1000, n=1000, meanA, meanB, meanC, cov=cov_pd)
# 1, 2, 3 - 99%, 4 - 100%

est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_nd)
# 1, 2, 3 - 63%, 4 - 81% 
est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_nd)
# 1, 2, 3 - 94%, 4 - 98% 


# plot the result ---------------------------------------------------------
# could have been organized as function 

# make data frames 
iev <- as.data.frame(
  cbind(est.power(N=1000, n=50, meanA, meanB, meanC, cov=cov_iev),
          est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_iev),
          est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_iev),
          est.power(N=1000, n=300, meanA, meanB, meanC, cov=cov_iev),
          est.power(N=1000, n=500, meanA, meanB, meanC, cov=cov_iev),
          est.power(N=1000, n=800, meanA, meanB, meanC, cov=cov_iev),
          est.power(N=1000, n=1000, meanA, meanB, meanC, cov=cov_iev)
  )
)

iuv <- as.data.frame(
  cbind(est.power(N=1000, n=50, meanA, meanB, meanC, cov=cov_iuv),
          est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_iuv),
          est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_iuv),
          est.power(N=1000, n=300, meanA, meanB, meanC, cov=cov_iuv),
          est.power(N=1000, n=500, meanA, meanB, meanC, cov=cov_iuv),
          est.power(N=1000, n=800, meanA, meanB, meanC, cov=cov_iuv),
          est.power(N=1000, n=1000, meanA, meanB, meanC, cov=cov_iuv)
  )
)

nd <- as.data.frame(
  cbind(est.power(N=1000, n=50, meanA, meanB, meanC, cov=cov_nd),
          est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_nd),
          est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_nd),
          est.power(N=1000, n=300, meanA, meanB, meanC, cov=cov_nd),
          est.power(N=1000, n=500, meanA, meanB, meanC, cov=cov_nd),
          est.power(N=1000, n=800, meanA, meanB, meanC, cov=cov_nd),
          est.power(N=1000, n=1000, meanA, meanB, meanC, cov=cov_nd)
  )
)

pd <- as.data.frame(
  cbind(est.power(N=1000, n=50, meanA, meanB, meanC, cov=cov_pd),
          est.power(N=1000, n=100, meanA, meanB, meanC, cov=cov_pd),
          est.power(N=1000, n=200, meanA, meanB, meanC, cov=cov_pd),
          est.power(N=1000, n=300, meanA, meanB, meanC, cov=cov_pd),
          est.power(N=1000, n=500, meanA, meanB, meanC, cov=cov_pd),
          est.power(N=1000, n=800, meanA, meanB, meanC, cov=cov_pd),
          est.power(N=1000, n=1000, meanA, meanB, meanC, cov=cov_pd)
  )
)

colnames(iev) <- c(50, 100, 200, 300, 500, 800, 1000)
colnames(iuv) <- c(50, 100, 200, 300, 500, 800, 1000)
colnames(pd) <- c(50, 100, 200, 300, 500, 800, 1000)
colnames(nd) <- c(50, 100, 200, 300, 500, 800, 1000)

# print result as LaTeX table 
xtable(t(iev)) 
xtable(t(iuv))
xtable(t(pd))
xtable(t(nd))

# plot result 
iev <- rbind(c(50, 100, 200, 300, 500, 800, 1000), iev)
rownames(iev)[1] <- "n"
iev <- as.data.frame(t(iev))
iuv <- rbind(c(50, 100, 200, 300, 500, 800, 1000), iuv)
rownames(iuv)[1] <- "n"
iuv <- as.data.frame(t(iuv))
pd <- rbind(c(50, 100, 200, 300, 500, 800, 1000), pd)
rownames(pd)[1] <- "n"
pd <- as.data.frame(t(pd))
nd <- rbind(c(50, 100, 200, 300, 500, 800, 1000), nd)
rownames(nd)[1] <- "n"
nd <- as.data.frame(t(nd))

par(mfrow=c(2,2))
plot(iev$n, iev$Pillai, type = "l", col = 1, main = "Independent and Equal", xlab = "observations (count)", ylab = "power")
lines(iev$n, iev$Wilks, type = "l", col = 2)
lines(iev$n, iev$HotellingLawley, type = "l", col = 3)
lines(iev$n, iev$Roy, type = "l", col = 4)

plot(iuv$n, iuv$Pillai, type = "l", col = 1, main = "Independent and Unequal", xlab = "observations (count)", ylab = "power")
lines(iuv$n, iuv$Wilks, type = "l", col = 2)
lines(iuv$n, iuv$HotellingLawley, type = "l", col = 3)
lines(iuv$n, iuv$Roy, type = "l", col = 4)

plot(pd$n, pd$Pillai, type = "l", col = 1, main = "Positively dependent", xlab = "observations (count)", ylab = "power")
lines(pd$n, pd$Wilks, type = "l", col = 2)
lines(pd$n, pd$HotellingLawley, type = "l", col = 3)
lines(pd$n, pd$Roy, type = "l", col = 4)

plot(nd$n, nd$Pillai, type = "l", col = 1, main = "Negatively dependent", xlab = "observations (count)", ylab = "power")
lines(nd$n, nd$Wilks, type = "l", col = 2)
lines(nd$n, nd$HotellingLawley, type = "l", col = 3)
lines(nd$n, nd$Roy, type = "l", col = 4)


# plot single trial data  -------------------------------------------------
# generate MV normal data - an example (means different for proposing, sigma=iev) 

n = 100 # as the sample size increases, the mean does not comes out properly. Why??
meanA = c(2.0,5.0,8.0) # mean vector for population A
meanB = c(2.1,5.1,8.1) # mean vector for population B
meanC = c(2.2,5.2,8.2) # mean vector for population C

A_iev = mvrnorm(n, mu=meanA, Sigma=cov_iev)
B_iev = mvrnorm(n, mu=meanB, Sigma=cov_iev)
C_iev = mvrnorm(n, mu=meanC, Sigma=cov_iev)
iev <- array(c(A_iev, B_iev, C_iev), dim=c(100, 3, 3))
dimnames(iev)[[3]] <- c("A", "B", "C")

# boxplot 
par(mfrow=c(3,1))
boxplot(iev[,1,], ylim=c(0,10), main="mu1")
boxplot(iev[,2,], ylim=c(0,10), main="mu2")
boxplot(iev[,3,], ylim=c(0,10), main="mu3")

# density plot 
iev.1.A <- density(iev[,1,1])
iev.2.A <- density(iev[,2,1])
iev.3.A <- density(iev[,3,1])

iev.1.B <- density(iev[,1,2])
iev.2.B <- density(iev[,2,2])
iev.3.B <- density(iev[,3,2])

iev.1.C <- density(iev[,1,3])
iev.2.C <- density(iev[,2,3])
iev.3.C <- density(iev[,3,3])

par(mfrow=c(3,1))
plot(iev.1.C, xlim=c(-2,12), ylim=c(0,0.6), col="red", main="Density plot of mu1", xlab=NA)
lines(iev.1.B, col="magenta")
lines(iev.1.A, col="blue") # how do I include legend here!?!?

plot(iev.2.C, xlim=c(-2,12), ylim=c(0,0.6), col="red", main="Density plot of mu2", xlab=NA)
lines(iev.2.B, col="magenta")
lines(iev.2.A, col="blue")

plot(iev.3.C, xlim=c(-2,12), ylim=c(0,0.6), col="red", main="Density plot of mu3", xlab=NA)
lines(iev.3.B, col="magenta")
lines(iev.3.A, col="blue")
