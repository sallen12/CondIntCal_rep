################################################################################
########################## Case study (STAR dataset) ###########################
################################################################################

################################################################################
## Set up

set.seed(789634)

source("utility_funcs.R") # source functions and packages needed to evaluate interval forecasts

# read interval forecasts
path <- "Case_Study/data/star/"
file_list <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
names(data_list) <- c("CQRNet", "Net (L)", "Net", "QNet", "RF (L)", "RF", "Ridge (L)", "Ridge")
data_list <- data_list[c("Ridge", "Ridge (L)", "RF", "RF (L)", "Net", "Net (L)", "CQRNet", "QNet")]

# check observations are the same for all methods
y <- data_list[[1]]$Obs
sapply(data_list, function(df) all.equal(y, df$Obs)) |> all()

alpha <- 0.1 # nominal coverage level


################################################################################
## IDR

# recalibrate interval forecasts using IDR
idr_all <- lapply(data_list, is_decomp, y = y, alpha = alpha)

# store recalibrated intervals in list
data_list_rc <- lapply(idr_all, function(x) {
  df <- x[["int_rc"]] |> as.data.frame()
  colnames(df) <- c("Lower", "Upper")
  return(df)
})


################################################################################
## Evaluate

# assess original interval forecasts
comp_all <- sapply(data_list, count_comparables) # percentage of distinct intervals that are comparable
cov_all <- sapply(data_list, coverage, y = y) # unconditional coverage
ocov_all <- sapply(data_list, coverage, y = y) # unconditional coverage of open interval forecasts
len_all <- sapply(data_list, ilength) # average length

# assess recalibrated interval forecasts
is_all <- sapply(idr_all, function(x) x[['decomp']][1] |> unname()) # interval score
cov_rc_all <- sapply(data_list_rc, coverage, y = y) # unconditional coverage
ocov_rc_all <- sapply(data_list_rc, coverage, y = y, closed = F) # unconditional coverage of open interval forecasts
len_rc_all <- sapply(data_list, ilength) # average length


##### display all values

rbind(is_all, ocov_all, cov_all, len_all, ocov_rc_all, cov_rc_all, len_rc_all) |> t()


##### plot decomposition terms

dcmp_all <- sapply(idr_all, function(x) x$decomp) |> t() |> as.data.frame()
dcmp_all$forecast <- names(data_list)
dcmp_all

mcbdsc(dcmp_all, MCB_lim = c(-0.001, 0.05), DSC_lim = c(-0.001, 0.02)) + ggtitle("STAR")
ggsave("Figures/fig2_star.png", height = 3.2, width = 4.5)

