library(purrr)
library(readr)
library(EnvStats)
library(isodistrreg)
library(scoringRules)

star_CQRNet_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_CQRNet_int0.csv",show_col_types = FALSE)
star_Net_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Net_int0.csv",show_col_types = FALSE)
star_Net_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Net_L_int0.csv",show_col_types = FALSE)
star_QNet_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_QNet_int0.csv",show_col_types = FALSE)
star_RF_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_RF_int0.csv",show_col_types = FALSE)
star_RF_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_RF_L_int0.csv",show_col_types = FALSE)
star_Ridge_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Ridge_int0.csv",show_col_types = FALSE)
star_Ridge_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Ridge_L_int0.csv",show_col_types = FALSE)

alpha=0.1


df_names<-c("star_CQRNet_int0","star_Net_int0","star_Net_L_int0","star_QNet_int0","star_RF_int0","star_RF_L_int0", "star_Ridge_int0","star_Ridge_L_int0")
m<-length(df_names)

#mean interval lengths
#lapply gibt liste aus, use unlist() to get a vector
length<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  mean(df[[2]]-df[[1]])
}))
print(paste("Mean lengths", length))


#mean difference li-yi, mean difference ui-yi
mean_dist_li_yi<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  mean(abs(df[[3]]-df[[1]]))
  
}))
mean_dist_ui_yi<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  mean(abs(df[[3]]-df[[2]]))
  
}))
print(paste("Mean distance of observations from lowerbound",mean_dist_li_yi))
print(paste("Mean distance of observations from upperbound",mean_dist_ui_yi))



#count how many y_i do not fall into [l_i, u_i]
coverage<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  1-sum(df[[3]]>df[[2]]|df[[3]]<df[[1]])/433
  
}))
print(paste("instances where y_i not in [l_i,u_i]",count))

#calculate IDR regression fit
fitted_names<-c("star_CQRNet_fitted","star_Net_fitted","star_Net_L_fitted","star_QNet_fitted","star_RF_fitted","star_RF_L_fitted", "star_Ridge_fitted","star_Ridge_L_fitted")
predict_names<-c("star_CQRNet_predict","star_Net_predict","star_Net_L_predict","star_QNet_predict","star_RF_predict","star_RF_L_predict", "star_Ridge_predict","star_Ridge_L_predict")
interval_names<-c("star_CQRNet_interval","star_Net_interval","star_Net_L_interval","star_QNet_interval","star_RF_interval","star_RF_L_interval", "star_Ridge_interval","star_Ridge_L_interval")

#m is the number of schemes
#in-sample predictions
for (i in 1:m){
  df<-get(df_names[i])
  assign(predict_names[i],predict(idr(y = df[[3]] ,
                                      X = data.frame(df[[1]],df[[2]]),
                                      groups = ,
                                      orders = c("icx"=1),
                                      stoch="sd"),data.frame(df[[1]],df[[2]])))
}

for (i in 1:m){
  assign(interval_names[i],cbind(qpred(get(predict_names[i]),alpha/2),qpred(get(predict_names[i]),1-alpha/2)))
}

#count how many y_i do not fall into \hat{l}_i, \hat{u}_i
coverage_new<-unlist(lapply(interval_names,function(name) {
  interval<-get(name)
  1-sum(y>interval[,2]|y<interval[,1])/433
}))
#mean lengths of updated intervals
length_new<-unlist(lapply(interval_names,function(name) {
  interval<-get(name)
  mean(interval[,2]-interval[,1])
}))


#calculate interval score of the old intervals, recalibrated intervals and empirical quantiles of the observations

IVS<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  mean(ints_quantiles(df[[3]], df[[1]], df[[2]], target_coverage = 1-alpha))
}))
print(IVS)

IVS_ISO_icx<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, interval[,1], interval[,2], target_coverage = 1-alpha))}))
print(IVS_ISO_icx)

IVS_mg<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, quantile(y,alpha/2),quantile(y,1-alpha/2), target_coverage = 1-alpha))}))
print(IVS_mg)

MCB_ISO_icx<-IVS-IVS_ISO_icx
DSC_ISO_icx<-IVS_mg-IVS_ISO_icx
UNC<-IVS_mg
#save all the results
sink(file="C:/Users/ju/Documents/icx_star_overview.txt")
decomp<-data.frame(IVS,IVS_ISO_icx,MCB_ISO_icx ,DSC_ISO_icx, UNC, row.names=df_names)
print(decomp)
data<-data.frame(coverage, coverage_new,length, length_new,row.names=df_names)
print(data)
sink(file=NULL)

cov_icx<-coverage_new
len_icx<-length_new

save(IVS_ISO_icx, MCB_ISO_icx,DSC_ISO_icx,cov_icx,len_icx, file ="C:/Users/ju/Documents/icx_star_MCB.txt")

