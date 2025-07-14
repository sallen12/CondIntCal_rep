library(EnvStats)
library(isodistrreg)
library(scoringRules)
library(calibrate)

#####TO ADAPT IF N OR ALPHA CHANGE##############################################
#define number of observations and determine alpha
set.seed(123)
n <- 1000
alpha <- 0.1


################################################################################
#define the observations
mu <- rnorm(n)
y <- rnorm(n,mean= mu,sd=1)

# Ideal forecaster: F = N(mu, 1)
L_id <- qnorm(alpha/2, mu)
U_id <- qnorm(1 - alpha/2, mu)
int_id <- cbind(L_id, U_id)

# Unconditional/Climatological forecaster: F = N(0, 2)
L_cl <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
U_cl <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)
int_cl <- cbind(L_cl, U_cl)

# Unfocused forecaster: F = 0.5 N(mu, 1) + 0.5 N(mu +tau, 1), tau=-1/+1 w.e.p
tau <- sample(c(-1, 1), n, replace = T)
#define CDF and inverse
F_un <- function(q, m, t) 0.5*pnorm(q, m) + 0.5*pnorm(q, m + t)
Q_un <- function(p, m, t) uniroot(function(x) F_un(x, m, t) - p, interval = c(-10, 10))$root
L_un <- sapply(1:n, function(i) Q_un(alpha/2, mu[i], tau[i]))
U_un <- sapply(1:n, function(i) Q_un(1 - alpha/2, mu[i], tau[i]))
int_un <- cbind(L_un, U_un)


# mean_based forecaster: F = N(mu+tau, 1),tau=-1/+1 w.e.p
tau <- sample(c(-1, 1), n, replace = T)
#define CDF and inverse
F_me <- function(q, m, t) pnorm(q, m + t)
Q_me <- function(p, m, t) uniroot(function(x) F_me(x, m, t) - p, interval = c(-10, 10))$root
L_me <- sapply(1:n, function(i) Q_me(alpha/2, mu[i], tau[i]))
U_me <- sapply(1:n, function(i) Q_me(1 - alpha/2, mu[i], tau[i]))
int_me <- cbind(L_me, U_me)

# sign-based forecaster: F = N(-mu, 1)
L_si <- qnorm(alpha/2,mean= -mu,sd=1)
U_si <- qnorm(1 - alpha/2,mean= -mu,sd=1)
int_si <- cbind(L_si, U_si)


# mixed (climatological/signed) forecaster: F =  N(0, 2)/N(-mu,1) with prob 1/2
L_1 <- qnorm(alpha/2, 0, sqrt(2)) |> rep(n)
U_1 <- qnorm(1 - alpha/2, 0, sqrt(2)) |> rep(n)
L_2 <- qnorm(alpha/2,mean= -mu,sd=1)
U_2 <- qnorm(1 - alpha/2,mean= -mu,sd=1)
r <- rbern(n, 1/2)
L_mi <- r*L_1+(1-r)*L_2
U_mi <- r*U_1+(1-r)*U_2
int_mi <- cbind(L_mi,U_mi)

#determine mean length and coverage of original intervals
cov_id <- mean(y >= L_id & y <= U_id)
cov_cl <- mean(y >= L_cl & y <= U_cl)
cov_un <- mean(y >= L_un & y <= U_un)
cov_me <- mean(y >= L_me & y <= U_me)
cov_si <- mean(y >= L_si & y <= U_si)
cov_mi <- mean(y >= L_mi & y <= U_mi)
length_id <- mean(U_id-L_id)
length_cl <- mean(U_cl-L_cl)
length_un <- mean(U_un-L_un)
length_me <- mean(U_me-L_me)
length_si <- mean(U_si-L_si)
length_mi <- mean(U_mi-L_mi)

coverage <- cbind(cov_id,cov_cl ,cov_un ,cov_me ,cov_si,cov_mi)
length <- cbind(length_id,length_cl ,length_un ,length_me ,length_si,length_mi)
names <- c("int_id","int_cl","int_un", "int_me","int_si","int_mi")
fit_names_comp <- c("y_1_fitted_comp","y_2_fitted_comp","y_3_fitted_comp","y_4_fitted_comp","y_5_fitted_comp","y_6_fitted_comp")
#On R^2, only the componentwise order and the increasing convex order are different
orders <- c("comp" = 1, "icx" = 1)

#fit IDR with respect to all 6 interval forecasts with regard to "componentwise" order
#plot the predicted distribution of the first observation

#Recall X needs to be a data.frame, y a numeric vector
for (k in 1:length(names))
{ interval<-get(names[k])
  interval<-as.data.frame(interval)
  assign(fit_names_comp[k], predict(idr(
    y = y ,
    X = interval,
    groups = ,
    orders = orders[1],
    stoch="sd"),interval))

  #plots first F_hat as index is set to 1
  plot(get(fit_names_comp[k]),index=1,ylim=c(0,1))}

#generate the new prediction intervals
interval_names_comp <- c("y_1_interval_comp","y_2_interval_comp","y_3_interval_comp","y_4_interval_comp","y_5_interval_comp","y_6_interval_comp")

for (k in 1:length(names))
{ prediction <- get(fit_names_comp[k])
assign(interval_names_comp[k],cbind(qpred(prediction,alpha/2),qpred(prediction,1-alpha/2)))}

#get new coverage and lengths
coverage_new<-cbind(cov_comp_id,cov_comp_cl,cov_comp_un,cov_comp_me,cov_comp_si,cov_comp_mi)
length_new<-cbind(length_comp_id,length_comp_cl,length_comp_un,length_comp_me,length_comp_si,length_comp_mi)
for(i in 1:6){
  assign(coverage_new[i], mean(y >= get(interval_names[i])[,1] & y <= get(interval_names[i])[,2]))
  assign(length_new[i], mean(get(interval_names[i])[,2]-get(interval_names[i])[,1]))}

len_cov <- rbind(coverage,coverage_new,length,length_new)
len_cov <- as.data.frame(len_cov,row.names=c("coverage","coverage_new","length","length_new"))
names(len_cov)[]<-c("ideal", "clima.", "unfocused", "mean-biased","signed","mixed")

#generate a summary output of all the in-sample IVS, MCB and DSC-components
summary<-matrix(0,6,5)
for (k in 1:length(names))
{ interval_comp <- get(interval_names_comp[k])
  old_interval <- get(names[k])

  IVS <- mean(ints_quantiles(y, old_interval[,1], old_interval[,2], target_coverage = 1-alpha))
  IVS_mg <- mean(ints_quantiles(y, quantile(y,probs=alpha/2,type=1),quantile(y,probs=1-alpha/2,type=1), target_coverage = 1-alpha))
  IVS_ISO_comp <- mean(ints_quantiles(y, interval_comp[,1], interval_comp[,2], target_coverage = 1-alpha))

  MCB_ISO_comp <- IVS-IVS_ISO_comp
  DSC_ISO_comp <- IVS_mg-IVS_ISO_comp

  summary[k,] <- cbind(IVS,IVS_ISO_comp,MCB_ISO_comp,DSC_ISO_comp,IVS_mg)}

summary <- as.data.frame(summary,row.names=c("ideal", "clima.", "unfocused", "mean-biased","signed","mixed"))
names(summary)[]<-c("IVS","IVS_ISO_comp","MCB_ISO_comp","DSC_ISO_comp","UNC_0")

######TO ADAPT IF N OR ALPHA CHANGED############################################
#sink variables "summary" and "len_cov" in separate file, save other variables if needed
sink(file="C:/Users/ju/Documents/overview_simulation_1000.txt")
print(summary)
print(len_cov)
sink(file=NULL)
#delete the numbers stored from the last for loop
remove(list=c("IVS","IVS_mg","IVS_ISO_comp","MCB_ISO_comp","DSC_ISO_comp"))
attach(summary)

IVS_1000<-IVS
IVS_mg_1000<-IVS_mg
IVS_ISO_comp_1000<-IVS_ISO_comp
MCB_ISO_comp_1000<-MCB_ISO_comp
DSC_ISO_comp_1000<-DSC_ISO_comp
#save the components under unique names
save(IVS_1000,IVS_ISO_comp_1000,IVS_mg_1000,MCB_ISO_comp_1000,DSC_ISO_comp_1000, file ="C:/Users/ju/Documents/n_1000.txt")
################################################################################
