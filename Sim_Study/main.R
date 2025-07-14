################################################################################
############################### Simulation study ###############################
################################################################################

################################################################################
## Set up

set.seed(123)

library(EnvStats)
library(isodistrreg)
library(scoringRules)
library(calibrate)
library(ggplot2)
source("Sim_Study/utility_funcs.R")

n <- 1000
alpha <- 0.1

tau <- sample(c(-1, 1), n, replace = T)
mu <- rnorm(n)

y <- rnorm(n, mean = mu, sd = 1) # simulate observations


get_ints(mu, tau, alpha) # get interval forecasts


################################################################################
## Evaluate


##### marginal coverage

cov_id <- coverage(y, int_id)
cov_cl <- coverage(y, int_cl)
cov_un <- coverage(y, int_un)
cov_me <- coverage(y, int_me)
cov_si <- coverage(y, int_si)
cov_mi <- coverage(y, int_mi)

cov_all <- c("Id." = cov_id, "Unc." = cov_cl, "Unf." = cov_un, "MB" = cov_me, "SB" = cov_si, "Mix" = cov_mi)


##### average length

len_id <- ilength(int_id)
len_cl <- ilength(int_cl)
len_un <- ilength(int_un)
len_me <- ilength(int_me)
len_si <- ilength(int_si)
len_mi <- ilength(int_mi)

len_all <- c("Id." = len_id, "Unc." = len_cl, "Unf." = len_un, "MB" = len_me, "SB" = len_si, "Mix" = len_mi)


##### interval score decompositions

dcmp_id <- is_decomp(y, int_id, alpha)
dcmp_cl <- is_decomp(y, int_cl, alpha)
dcmp_un <- is_decomp(y, int_un, alpha)
dcmp_me <- is_decomp(y, int_me, alpha)
dcmp_si <- is_decomp(y, int_si, alpha)
dcmp_mi <- is_decomp(y, int_mi, alpha)

is_all <- c("Id." = dcmp_id$decomp[1], "Unc." = dcmp_cl$decomp[1], "Unf." = dcmp_un$decomp[1], "MB" = dcmp_me$decomp[1], "SB" = dcmp_si$decomp[1], "Mix" = dcmp_mi$decomp[1])
is_all <- is_all |> unname()


##### marginal coverage of recalibrated intervals

cov_rc_id <- coverage(y, dcmp_id$int_rc)
cov_rc_cl <- coverage(y, dcmp_cl$int_rc)
cov_rc_un <- coverage(y, dcmp_un$int_rc)
cov_rc_me <- coverage(y, dcmp_me$int_rc)
cov_rc_si <- coverage(y, dcmp_si$int_rc)
cov_rc_mi <- coverage(y, dcmp_mi$int_rc)

cov_rc_all <- c("Id." = cov_rc_id, "Unc." = cov_rc_cl, "Unf." = cov_rc_un, "MB" = cov_rc_me, "SB" = cov_rc_si, "Mix" = cov_rc_mi)


##### average length of recalibrated intervals

len_rc_id <- ilength(dcmp_id$int_rc)
len_rc_cl <- ilength(dcmp_cl$int_rc)
len_rc_un <- ilength(dcmp_un$int_rc)
len_rc_me <- ilength(dcmp_me$int_rc)
len_rc_si <- ilength(dcmp_si$int_rc)
len_rc_mi <- ilength(dcmp_mi$int_rc)

len_rc_all <- c("Id." = len_rc_id, "Unc." = len_rc_cl, "Unf." = len_rc_un, "MB" = len_rc_me, "SB" = len_rc_si, "Mix" = len_rc_mi)


##### display all values

rbind(is_all, cov_all, len_all, cov_rc_all, len_rc_all) |> t()


##### plot decomposition terms

dcmp_all <- rbind(dcmp_id$decomp, dcmp_cl$decomp, dcmp_un$decomp, dcmp_me$decomp, dcmp_si$decomp, dcmp_mi$decomp) |> as.data.frame()
dcmp_all$forecast <- c("Ideal", "Uncon.", "Unfoc.", "MB", "SB", "Mixed")

mcbdsc(dcmp_all, MCB_lim = c(0, 10), DSC_lim = c(0, 2.5))



