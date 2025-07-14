library(purrr)
library(readr)
library(EnvStats)
library(isodistrreg)
library(scoringRules)
library(ggplot2)
library(ggrepel)

#read in all the generated interval forecasts
star_CQRNet_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_CQRNet_int0.csv",show_col_types = FALSE)
star_Net_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Net_int0.csv",show_col_types = FALSE)
star_Net_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Net_L_int0.csv",show_col_types = FALSE)
star_QNet_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_QNet_int0.csv",show_col_types = FALSE)
star_RF_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_RF_int0.csv",show_col_types = FALSE)
star_RF_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_RF_L_int0.csv",show_col_types = FALSE)
star_Ridge_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Ridge_int0.csv",show_col_types = FALSE)
star_Ridge_L_int0 <- read_csv("C:/Users/ju/Downloads/star0/star_Ridge_L_int0.csv",show_col_types = FALSE)



label<-c("CQRNet","Net","Net_L","QNet","RF","RF_L", "Ridge","Ridge_L")
df_names<-c("star_CQRNet_int0","star_Net_int0","star_Net_L_int0","star_QNet_int0","star_RF_int0","star_RF_L_int0", "star_Ridge_int0","star_Ridge_L_int0")

#m is number of methods,num_obs is number of observations (i.e. n)
m<-length(df_names)
num_obs<-433
alpha=0.1

#save the observations and check if all are the same.
y<-star_CQRNet_int0[[3]]
for (i in 1:m){
  df<-get(df_names[i])
 if(identical(y,df[[3]])){print("true ")} 
}
#Count comparable pairs of intervals in each data set, yields a max of n^2 comparable pairs unless some datapoints are non distinct
pairs<-matrix(0, m, 3)
for( l in 1:m){
  df<-get(df_names[l])
comp=0
icx=0
rep=0
for (i in 1:num_obs)
{for(j in 1:num_obs)
{if((df[[1]][i] <= df[[1]][j]) && (df[[2]][i] <= df[[2]][j])){comp<-comp+1}
  if((df[[1]][i] <= df[[1]][j]) && ((df[[1]][i]+df[[2]][i])<=(df[[1]][j]+df[[2]][j]))){icx<-icx+1}
  if((df[[1]][i] == df[[1]][j]) && (df[[2]][i]==df[[2]][j])){rep<-rep+1}
}}
pairs[l,]<-rbind(comp,icx,rep)
}

percentages<-(pairs[,1:2]*2-433)/(num_obs^2)
pairs<-as.data.frame(cbind(2*pairs-433,percentages),row.names=label)
colnames(pairs)<-c("componentwise","icx","repetition","percentage comp","percentage icx")

#create overview over comparable pairs
sink(file="C:/Users/ju/Documents/comparable_pairs.txt")
print(paste("There are ",num_obs^2,"possible ordered pairs. If x_i>=x_j and x_i<=x_j when i=j, 2 are added,x_i=x_i 1 is added"))
print(pairs)
sink(file=NULL)

#calculate mean interval lengths
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
coverage<-unlist(lapply(df_names,function(name) {
  df <- get(name)
  sum(df[[3]]>df[[2]]|df[[3]]<df[[1]])
       
}))
print(paste("instances where y_i not in [l_i,u_i]",coverage))

#calculate IDR regression fit using the "isodistrreg" package and from that the functions "idr", "predict" and "qpred"
predict_names<-c("star_CQRNet_pred","star_Net_pred","star_Net_L_pred","star_QNet_pred","star_RF_pred","star_RF_L_pred", "star_Ridge_pred","star_Ridge_L_pred")
interval_names<-c("star_CQRNet_fitted","star_Net_fitted","star_Net_L_fitted","star_QNet_fitted","star_RF_fitted","star_RF_L_fitted", "star_Ridge_fitted","star_Ridge_L_fitted")
fitted_names<-c("star_CQRNet_fitted","star_Net_fitted","star_Net_L_fitted","star_QNet_fitted","star_RF_fitted","star_RF_L_fitted", "star_Ridge_fitted","star_Ridge_L_fitted")

for (i in 1:m){
    df<-get(df_names[i])
    assign(predict_names[i],predict(idr(y = df[[3]] ,
                               X = data.frame(df[[1]],df[[2]]),
                               groups = ,
                               orders = c("comp"=1),
                               stoch="sd"),data.frame(df[[1]],df[[2]])))

   assign(interval_names[i],cbind(qpred(get(predict_names[i]),alpha/2),qpred(get(predict_names[i]),1-alpha/2)))
}

#count how many y_i do not fall into [\hat{l}_i, \hat{u}_i]
coverage_new<-unlist(lapply(interval_names,function(name) {
  interval<-get(name)
  sum(y>interval[,2]|y<interval[,1])
}))
#mean lengths of updated intervals
length_new<-unlist(lapply(interval_names,function(name) {
  interval<-get(name)
  mean(interval[,2]-interval[,1])
}))


#calculate interval score of the original intervals (df[[3]] is the y_i), of updated intervals and of the empirical quantiles of the observations
IVS<-unlist(lapply(df_names,function(name) {
     df <- get(name)
     mean(ints_quantiles(df[[3]], df[[1]], df[[2]], target_coverage = 1-alpha))
}))
print(IVS)

#calculate recalibrated interval scores
IVS_ISO<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, interval[,1], interval[,2], target_coverage = 1-alpha))}))
print(IVS_ISO)

#calculate interval score of the empirical quantiles
IVS_mg<-unlist(lapply(interval_names,function(name) {
  interval <- get(name)
  mean(ints_quantiles(y, quantile(y,alpha/2),quantile(y,1-alpha/2), target_coverage = 1-alpha))}))
print(IVS_mg)


MCB_ISO<-IVS-IVS_ISO
DSC_ISO<-IVS_mg-IVS_ISO
UNC<-IVS_mg

#load the icx version and the 1dim version
load("C:/Users/ju/Documents/icx_star_MCB.txt")
load("C:/Users/ju/Documents/1dim_star_MCB.txt")

#generate file with interval scores, MCB- and DSC-components
coverage<-1-coverage/433
coverage_new<-1-coverage_new/433
sink(file="C:/Users/ju/Documents/comp_star_overview.txt")
decomp<-data.frame(IVS,IVS_ISO,MCB_ISO,DSC_ISO,UNC, row.names=df_names)
print(decomp)
data<-data.frame(coverage, coverage_new,cov_icx,cov_1dim, row.names=df_names)
data_2<-data.frame(length, length_new,len_icx,len_1dim, row.names=df_names)
names(data)<-c("coverage", "coverage_comp","coverage_icx","coverage_1dim")
names(data_2)<-c("length", "length_comp","length_icx","length_1dim")
print(data)
print(data_2)
sink(file=NULL)

#how big are the differences between "comp" and "icx"
print(abs(MCB_ISO-MCB_ISO_icx))
print(abs(DSC_ISO-DSC_ISO_icx))
#how big are the differences between "comp" and "1dim"
print(abs(MCB_ISO-MCB_ISO_1dim))
print(abs(DSC_ISO-DSC_ISO_1dim))

#plot MCB-DSC of all 8 methods regarding the IDR-methods "comp","icx" and "1dim" in black, green and blue

#for easier legibility only plot labels in green/blue if there is a visible difference
labelshort<-c("CQRNet","","Net_L","QNet","","" ,"","Ridge_L")

sample_data1<-data.frame(MCB_ISO,DSC_ISO,label)
sample_data2<-data.frame(MCB_ISO_icx,DSC_ISO_icx,labelshort)
sample_data3<-data.frame(MCB_ISO_1dim,DSC_ISO_1dim,labelshort)

#all three
Gplot<-ggplot() +
  geom_point(data=sample_data1, aes(x=MCB_ISO, y=DSC_ISO),col="black") +
  geom_text_repel(data=sample_data1,
                  aes(x=MCB_ISO, y=DSC_ISO,label=label),
    nudge_x=c(0.01,0.01,0.01,0.01,-0.01,-0.01,-0.01,0.01), nudge_y=c(0,0,0,0,-0.01,-0.01,-0.005,0) ) +
  geom_point(data=sample_data2, aes(x=MCB_ISO_icx,y=DSC_ISO_icx),col="green")+
  geom_text_repel(data=sample_data2,
                  aes(x=MCB_ISO_icx, y=DSC_ISO_icx,label=labelshort),
  nudge_x=-0.01, nudge_y=0,col="green")+
  geom_point(data=sample_data3, aes(x=MCB_ISO_1dim,y=DSC_ISO_1dim),col="blue")+
  geom_text_repel(data=sample_data3,
                  aes(x=MCB_ISO_1dim, y=DSC_ISO_1dim,label=labelshort),
    nudge_x=0, nudge_y=-0.01,col="blue")+
labs(x = "MCB_ISO", y = "DSC_ISO")

#add isolines
  for(i in -10:10)
  {Gplot<-Gplot+geom_abline(intercept =-0.01*i, slope = 1,alpha=0.4)}
show(Gplot)

#plot only with regard to the componentwise order
Gplot4<-ggplot() +
  geom_point(data=sample_data1, aes(x=MCB_ISO, y=DSC_ISO),col="black") +
  geom_text_repel(data=sample_data1,
                  aes(x=MCB_ISO, y=DSC_ISO,label=label))
#add isolines
for(i in -10:10)
{Gplot4<-Gplot4+geom_abline(intercept =-0.01*i, slope = 1,alpha=0.4)}
show(Gplot4)
