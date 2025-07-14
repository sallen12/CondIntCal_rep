library(isodistrreg)
library(scoringRules)
set.seed(123)
n <- 1000
alpha <- 0.2
q<-qnorm(0.9)

#define the observations: Y from N(0,10)/N(0,100) w.e.p.
x <- rbinom(n,1,0.5)
y <- rnorm(n,mean= 0,sd=10)*x+(1-x)*rnorm(n, mean=0,sd=100)

# Ideal forecaster: F = N(mu, 1)
#there are only 2 different types of interval forecasts that are nested
#all intervals are comparable w.r.t. "icx" order
L <- qnorm(alpha/2,mean= 0,sd=10)*x+(1-x)*qnorm(alpha/2, mean=0,sd=100)
U <- qnorm(1-alpha/2,mean= 0,sd=10)*x+(1-x)*qnorm(1-alpha/2, mean=0,sd=100)
interval <- cbind(L, U)
print(cbind(L,U)[1:20])

#fit IDR with icx order and calculate recalibrated intervals
orders <- c("comp" = 1, "icx" = 1)
interval <- as.data.frame(interval)
fit <- predict(idr(
  y = y ,
  X = interval,
  groups = ,
  orders = orders[2],
  stoch="sd"),interval)
int_new<-cbind(qpred(fit,alpha/2),qpred(fit,1-alpha/2))

#plot two intervals, one belongs to group with 100 sd, one with 10 sd
plot(fit,index=1,ylim=c(0,1))
plot(fit,index=2,ylim=c(0,1),add=TRUE)
#calculate interval score and components
IVS <- mean(ints_quantiles(y, L, U, target_coverage = 1 - alpha))               
IVS_mg <- mean(ints_quantiles(y, quantile(y,alpha/2),quantile(y,1-alpha/2), target_coverage = 1-alpha))
IVS_ISO <- mean(ints_quantiles(y, int_new[,1], int_new[,2], target_coverage = 1-alpha))
MCB_ISO <- IVS-IVS_ISO
DSC_ISO <- IVS_mg-IVS_ISO
sum <- as.data.frame(cbind(IVS,IVS_mg,IVS_ISO,MCB_ISO,DSC_ISO))
print(sum)
#adapted forecaster for midpoint order
#here we need to have different means for the midpoints of two sets to be different
L2 <- qnorm(alpha/2,mean= 0.1,sd=10)*x+(1-x)*qnorm(alpha/2, mean=0,sd=100)
U2 <- qnorm(1-alpha/2,mean= 0.1,sd=10)*x+(1-x)*qnorm(1-alpha/2, mean=0,sd=100)
midpoint <- (L2+U2)/2