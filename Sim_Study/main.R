################################################################################
############################### Simulation study ###############################
################################################################################

################################################################################
## Set up

set.seed(1234)

source("utility_funcs.R") # source functions and packages needed to evaluate interval forecasts

n <- 1000 # sample size
alpha <- 0.1 # nominal coverage level

tau <- sample(c(-1, 1), n, replace = T)
mu <- rnorm(n)
y <- rnorm(n, mean = mu, sd = 1) # simulate observations

data_list <- get_ints(mu, tau, alpha) # get interval forecasts


################################################################################
## IDR

idr_all <- lapply(data_list, is_decomp, y = y, alpha = alpha) # recalibrate interval forecasts using IDR

data_list_rc <- lapply(idr_all, function(x) x[["int_rc"]]) # store recalibrated intervals in list


################################################################################
## Evaluate

# assess original interval forecasts
cov_all <- sapply(data_list, coverage, y = y) # unconditional coverage
ocov_all <- sapply(data_list, coverage, y = y, closed = F) # unconditional coverage of open interval forecasts
len_all <- sapply(data_list, ilength) # average length

# assess recalibrated interval forecasts
is_all <- sapply(idr_all, function(x) x[['decomp']][1] |> unname()) # interval score
cov_rc_all <- sapply(data_list_rc, coverage, y = y) # unconditional coverage
ocov_rc_all <- sapply(data_list_rc, coverage, y = y, closed = F) # unconditional coverage of open interval forecasts
len_rc_all <- sapply(data_list, ilength) # average length


##### display all values

rbind(is_all, ocov_all, cov_all, len_all, ocov_rc_all, cov_rc_all, len_rc_all) |> t()


##### plot decomposition terms

dcmp_all <- sapply(idr_all, function(x) x[['decomp']]) |> t() |> as.data.frame()
dcmp_all$forecast <- rownames(dcmp_all)
dcmp_all

mcbdsc(dcmp_all, MCB_lim = c(-0.01, 9), DSC_lim = c(-0.01, 2.5)) + ggtitle(bquote(alpha == .(alpha)))
ggsave(paste0("Figures/fig1_", 100*(1 - alpha), ".png"), height = 3.2, width = 4.5)

