library(purrr)
library(readr)
library(EnvStats)
library(isodistrreg)
library(scoringRules)
library(ggplot2)
library(ggrepel)

#read in all the data
star_CQRNet_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_CQRNet_int0.csv",show_col_types = FALSE)
star_Net_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Net_int0.csv",show_col_types = FALSE)
star_Net_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Net_L_int0.csv",show_col_types = FALSE)
star_QNet_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_QNet_int0.csv",show_col_types = FALSE)
star_RF_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_RF_int0.csv",show_col_types = FALSE)
star_RF_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_RF_L_int0.csv",show_col_types = FALSE)
star_Ridge_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Ridge_int0.csv",show_col_types = FALSE)
star_Ridge_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Ridge_L_int0.csv",show_col_types = FALSE)


alpha=0.1
label<-c("CQRNet","Net","Net_L","QNet","RF","RF_L", "Ridge","Ridge_L")
df_names<-c("star_CQRNet_int0","star_Net_int0","star_Net_L_int0","star_QNet_int0","star_RF_int0","star_RF_L_int0", "star_Ridge_int0","star_Ridge_L_int0")
m<-length(df_names)


#compute the means of the intervals and save as dataframes
for (i in 1:m){
  df<-get(df_names[i])
  assign(label[i],as.data.frame(cbind((df[[1]]+df[[2]])/2,df[[3]])))
}

#use IDR to generate recalibrated intervals
predict_names<-c("star_CQRNet_predict","star_Net_predict","star_Net_L_predict","star_QNet_predict","star_RF_predict","star_RF_L_predict", "star_Ridge_predict","star_Ridge_L_predict")
interval_names<-c("star_CQRNet_interval","star_Net_interval","star_Net_L_interval","star_QNet_interval","star_RF_interval","star_RF_L_interval", "star_Ridge_interval","star_Ridge_L_interval")
#m is the number of schemes
#in-sample predictions
for (i in 1:m){
  df<-get(label[i])
  assign(predict_names[i],predict(idr(y = df[[2]] ,
                                      X = as.data.frame(df[[1]]),
                                      groups = ,
                                      orders = c("comp"=1),
                                      stoch="sd"),as.data.frame(df[[1]])))
}

for (i in 1:m){
  assign(interval_names[i],cbind(qpred(get(predict_names[i]),alpha/2),qpred(get(predict_names[i]),1-alpha/2)))
}

#count how many y_i do not fall into \hat{l}_i, \hat{u}_i
count_new<-unlist(lapply(interval_names,function(name) {
  interval<-get(name)
  sum(y>interval[,2]|y<interval[,1])
}))
#mean lengths of updated intervals
length_new<-unlist(lapply(interval_names,function(name) {
  interval<-get(name)
  mean(interval[,2]-interval[,1])
}))


#calculate interval score of the old intervals, recalibrated intervals and of empirical quantiles of the observations

IVS<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  mean(ints_quantiles(df[[3]], df[[1]], df[[2]], target_coverage = 1-alpha))
}))
print(IVS)

IVS_ISO_1dim<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, interval[,1], interval[,2], target_coverage = 1-alpha))}))
print(IVS_ISO_1dim)

IVS_mg<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, quantile(y,alpha/2),quantile(y,1-alpha/2), target_coverage = 1-alpha))}))
print(IVS_mg)

MCB_ISO_1dim<-IVS-IVS_ISO_1dim
DSC_ISO_1dim<-IVS_mg-IVS_ISO_1dim
UNC<-IVS_mg

#save the variables in seperate file
sink(file="C:/Users/ju/Documents/1dim_star_overview.txt")

decomp<-data.frame(IVS,IVS_ISO_1dim,MCB_ISO_1dim,DSC_ISO_1dim,UNC, row.names=df_names)
print(decomp)

sink(file=NULL)
cov_1dim<-1-count_new/433
len_1dim<-length_new

save(IVS_ISO_1dim,MCB_ISO_1dim,DSC_ISO_1dim,cov_1dim,len_1dim, file ="C:/Users/ju/Documents/1dim_star_MCB.txt")


