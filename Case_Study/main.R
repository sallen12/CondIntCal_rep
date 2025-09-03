################################################################################
################################# Case study ###################################
################################################################################

################################################################################
## Set up

set.seed(789634)

#devtools::install_github("sallen12/CondIntCal")
library(CondIntCal)

# variable (one of "star", "bike", "facebook_1")
var <- "star"

# read interval forecasts
path <- paste0("Case_Study/data/", var)
file_list <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(file_list, read.csv)
names(data_list) <- c("CQRNet", "Net", "Net (L)", "QNet", "RF", "RF (L)", "Ridge", "Ridge (L)")
data_list <- data_list[c("Ridge", "Ridge (L)", "RF", "RF (L)", "Net", "Net (L)", "CQRNet", "QNet")]

# check observations are the same for all methods
y <- data_list[[1]]$Obs
sapply(data_list, function(df) all.equal(y, df$Obs)) |> all()
data_list <- lapply(data_list, function(x) {x$Obs <- NULL; x})

alpha <- 0.1 # nominal coverage level


################################################################################
## IDR

# recalibrate interval forecasts using IDR
idr_all <- lapply(data_list, is_decomp, y = y, level = 1 - alpha, return_fit = T)

# store recalibrated intervals in list
data_list_rc <- lapply(idr_all, function(x) x[["int_rc"]])


################################################################################
## Evaluate

# assess original interval forecasts
comp_all <- sapply(data_list, count_comparables) # percentage of distinct intervals that are comparable
cov_all <- sapply(data_list, coverage, y = y) # unconditional coverage
ocov_all <- sapply(data_list, coverage, y = y, closed = F) # unconditional coverage of open interval forecasts
len_all <- sapply(data_list, ilength) # average length

# assess recalibrated interval forecasts
is_all <- sapply(idr_all, function(x) x[['decomp']][1] |> unname()) # interval score
cov_rc_all <- sapply(data_list_rc, coverage, y = y) # unconditional coverage
ocov_rc_all <- sapply(data_list_rc, coverage, y = y, closed = F) # unconditional coverage of open interval forecasts
len_rc_all <- sapply(data_list_rc, ilength) # average length


##### display all values

rbind(comp_all*100, is_all, cov_all, len_all, ocov_rc_all, cov_rc_all, len_rc_all) |> t() |> round(2)


##### plot decomposition terms

dcmp_all <- sapply(idr_all, function(x) x$decomp) |> t() |> as.data.frame()
dcmp_all$forecast <- rownames(dcmp_all)
dcmp_all

if (var == "star") {
  MCB_lim <- c(-0.001, 0.08)
  DSC_lim <- c(-0.001, 0.06)
  tit <- "STAR"
} else if (var == "bike") {
  MCB_lim <- c(-0.001, 1)
  DSC_lim <- c(-0.001, 3)
  tit <- "Bike"
} else if (var == "facebook_1") {
  MCB_lim <- c(-0.001, 8)
  DSC_lim <- c(-0.001, 12)
  tit <- "Facebook"
}
plot_mcbdsc(dcmp_all, MCB_lim = MCB_lim, DSC_lim = DSC_lim) + ggplot2::ggtitle(tit) + ggplot2::theme(aspect.ratio = 1)
ggplot2::ggsave(paste0("Figures/fig2_", var,".png"), height = 4, width = 4)

