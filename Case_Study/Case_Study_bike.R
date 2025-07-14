library(purrr)
library(readr)
library(EnvStats)
library(isodistrreg)
library(scoringRules)
library(ggplot2)
library(ggrepel)

#load the data from external files
bike_CQRNet_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_CQRNet_int0.csv",show_col_types = FALSE)
bike_Net_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_Net_int0.csv",show_col_types = FALSE)
bike_Net_L_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_Net-L_int0.csv",show_col_types = FALSE)
bike_QNet_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_QNet_int0.csv",show_col_types = FALSE)
bike_RF_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_RF_int0.csv",show_col_types = FALSE)
bike_RF_L_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_RF-L_int0.csv",show_col_types = FALSE)
bike_Ridge_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_Ridge_int0.csv",show_col_types = FALSE)
bike_Ridge_L_int0 <- read_csv("C:/Users/ju/Downloads/bike0/bike_Ridge-L_int0.csv",show_col_types = FALSE)



df_names<-c("bike_CQRNet_int0","bike_Net_int0","bike_Net_L_int0","bike_QNet_int0","bike_RF_int0","bike_RF_L_int0", "bike_Ridge_int0","bike_Ridge_L_int0")
m<-length(df_names) #numbers of methods
label<-c("CQRNet","Net","Net_L","QNet","RF","RF_L", "Ridge","Ridge_L")
df_length<-length(bike_CQRNet_int0[[3]]) #number of observations

alpha=0.1

#check if the observations are the same for all methods
y<-bike_CQRNet_int0[[3]]
for (i in 1:m){
  df<-get(df_names[i])
  if(identical(y,df[[3]])){print("true ")} 
}


#mean interval lengths
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



#count how many y_i do not fall into l_i, u_i
count<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  sum(df[[3]]>df[[2]]|df[[3]]<df[[1]])
}))
print(paste("instances where y_i not in [l_i,u_i]",count))

#calculate IDR regression fit using the "isodistrreg"-package and from that functions "idr", "predict" and "qpred"
fitted_names<-c("facebook_CQRNet_fitted","facebook_Net_fitted","facebook_Net_L_fitted","facebook_QNet_fitted","facebook_RF_fitted","facebook_RF_L_fitted", "facebook_Ridge_fitted","facebook_Ridge_L_fitted")
predict_names<-c("facebook_CQRNet_predict","facebook_Net_predict","facebook_Net_L_predict","facebook_QNet_predict","facebook_RF_predict","facebook_RF_L_predict", "facebook_Ridge_predict","facebook_Ridge_L_predict")
interval_names<-c("facebook_CQRNet_interval","facebook_Net_interval","facebook_Net_L_interval","facebook_QNet_interval","facebook_RF_interval","facebook_RF_L_interval", "facebook_Ridge_interval","facebook_Ridge_L_interval")

#m is the number of schemes
#in-sample predictions
for (i in 1:m){
  df<-get(df_names[i])
  assign(predict_names[i],predict(idr(y = df[[3]] ,
                                      X = data.frame(df[[1]],df[[2]]),
                                      groups = ,
                                      orders = c("comp"=1),
                                      stoch="sd"),data.frame(df[[1]],df[[2]])))
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


#calculate interval score of the old intervals (df[[3]] is the y_i), the new intervals and of the empirical quantiles of the observations

IVS<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  mean(ints_quantiles(df[[3]], df[[1]], df[[2]], target_coverage = 1-alpha))
}))
print(IVS)

IVS_ISO<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, interval[,1], interval[,2], target_coverage = 1-alpha))}))
print(IVS_ISO)

IVS_mg<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, quantile(y,alpha/2),quantile(y,1-alpha/2), target_coverage = 1-alpha))}))
print(IVS_mg)


MCB_ISO<-IVS-IVS_ISO
DSC_ISO<-IVS_mg-IVS_ISO
UNC<-IVS_mg

#load the results from the "icx" and midpoint order
load("C:/Users/ju/Documents/icx_bike_MCB.txt")
load("C:/Users/ju/Documents/1dim_bike_MCB.txt")

coverage<-1-count/2178
coverage_new<-1-count_new/2178

sink(file="C:/Users/ju/Documents/comp_bike_overview.txt")
decomp<-data.frame(IVS,IVS_ISO,MCB_ISO,DSC_ISO,UNC, row.names=df_names)
print(decomp)
data<-data.frame(coverage, coverage_new,coverage_icx,coverage_1dim, row.names=df_names)
data_2<-data.frame(length, length_new,length_icx,length_1dim, row.names=df_names)
names(data)<-c("coverage", "coverage_comp","coverage_icx","coverage_1dim")
names(data_2)<-c("length", "length_comp","length_icx","length_1dim")
print(data)
print(data_2)
sink(file=NULL)

#plot the results in an MCB/DSC-plot
labelshort<-c("CQRNet","","Net_L","QNet","","RF_L","" ,"Ridge_L")

sample_data1<-data.frame(MCB_ISO,DSC_ISO,label)
sample_data2<-data.frame(MCB_ISO_icx,DSC_ISO_icx,labelshort)
sample_data3<-data.frame(MCB_ISO_1dim,DSC_ISO_1dim,labelshort)

#plot all three
Gplot<-ggplot() +
  geom_point(data=sample_data1, aes(x=MCB_ISO, y=DSC_ISO),col="black") +
  geom_text_repel(data=sample_data1,
                  aes(x=MCB_ISO, y=DSC_ISO,label=label),
                  nudge_x=0.2, nudge_y=0.3 ) +
  geom_point(data=sample_data2, aes(x=MCB_ISO_icx,y=DSC_ISO_icx),col="green")+
  geom_text_repel(data=sample_data2,
                  aes(x=MCB_ISO_icx, y=DSC_ISO_icx,label=labelshort),
                  nudge_x=-0.1, nudge_y=0.6,col="green")+
  geom_point(data=sample_data3, aes(x=MCB_ISO_1dim,y=DSC_ISO_1dim),col="blue")+
  geom_text_repel(data=sample_data3,
                  aes(x=MCB_ISO_1dim, y=DSC_ISO_1dim,label=labelshort),
                  nudge_x=0, nudge_y=-0.3,col="blue")+
  labs(x = "MCB_ISO", y = "DSC_ISO")
#add isolines
for(i in -30:30)
{Gplot<-Gplot+geom_abline(intercept =-0.3*i, slope = 1,alpha=0.4)}

show(Gplot)

#all three without Ridge und Ridge_L
Gplot2<-Gplot+
  ggplot2::scale_y_continuous(limits = c(2, 4), expand = c(0, 0))+
  ggplot2::scale_x_continuous(limits = c(-0.1, 0.7), expand = c(0, 0))
show(Gplot2)

#only componentwise order
Gplot3<-ggplot() +
  geom_point(data=sample_data1, aes(x=MCB_ISO, y=DSC_ISO),col="black") +
  geom_text_repel(data=sample_data1,
                  aes(x=MCB_ISO, y=DSC_ISO,label=label))
#add isolines
for(i in -30:30)
{Gplot3<-Gplot3+geom_abline(intercept =-0.3*i, slope = 1,alpha=0.4)}
show(Gplot3)
