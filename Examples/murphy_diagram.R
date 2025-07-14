
library(forecast)
library(purrr)
library(readr)
library(murphydiagram)
library(MetricsWeighted)

#define the observations
n=1000
mu <- rnorm(n)
y <- rnorm(n,mean= mu,sd=1)


# Ideal forecaster: F = N(mu, 1)
L_id <- qnorm(alpha/2, mu)
U_id <- qnorm(1 - alpha/2, mu)


# Unconditional/Climatological forecaster: F = N(0, 2)
L_cl <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
U_cl <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)

# Unfocused forecaster: F = 0.5N(mu, 1) + 0.5N(mu +/- 1, 1)
tau <- sample(c(-1, 1), n, replace = T)
F_un <- function(q, m, t) 0.5*pnorm(q, m) + 0.5*pnorm(q, m + t)
#used to find quantiles, find zero in a certain interval, of this function take the root
Q_un <- function(p, m, t) uniroot(function(x) F_un(x, m, t) - p, interval = c(-10, 10))$root
L_un <- sapply(1:n, function(i) Q_un(alpha/2, mu[i], tau[i]))
U_un <- sapply(1:n, function(i) Q_un(1 - alpha/2, mu[i], tau[i]))

#define a function for each forecaster that is the weighted sum of the elementary scores
f_id<-function(x){ elementary_score_quantile(
  actual=y,
  predicted=L_id,
  w = NULL,
  alpha = 0.05,
  theta = x)
  +elementary_score_quantile(
    actual=y,
    predicted=U_id,
    w = NULL,
    alpha = 0.95,
    theta = x)
}
f_cl<-function(x){ elementary_score_quantile(
  actual=y,
  predicted=L_cl,
  w = NULL,
  alpha = 0.05,
  theta = x)
  +elementary_score_quantile(
    actual=y,
    predicted=U_cl,
    w = NULL,
    alpha = 0.95,
    theta = x)
}
f_un<-function(x){ elementary_score_quantile(
  actual=y,
  predicted=L_un,
  w = NULL,
  alpha = 0.05,
  theta = x)
  +elementary_score_quantile(
    actual=y,
    predicted=U_un,
    w = NULL,
    alpha = 0.95,
    theta = x)
}
#plot the three functions
x<-seq(length.out=1000,from=-5,to=5 )
plot(x,map(x,f_id),col="blue",type="l",xlab=expression(theta),ylab=expression(paste("S","_",theta,",",alpha)),ylim=c(0,0.05))

lines(x,map(x,f_cl),col="red",type="l")
lines(x,map(x,f_un),col="green",type="l")
legend("topleft", legend=c("Ideal", "Climatological","Unfocused"),
       col=c("blue","red","green"), lty=1, cex=0.8)



