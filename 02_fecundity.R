########## Fecundity analyzes ########## 
########################################

setwd("C:/Users/sibazin/Desktop/Simon/phd/chapitre_2/article/data")
#--------------------------------------------------

# This script analyzes fish fecundity

#--------------------------------------------------

# Packages
library(ggplot2)
library(plyr)
library(lmerTest)
library(car)
library(multcomp)
library(stringr)
library(gridExtra)

# Data
fec<-read.csv(file = "df_fecundity.csv",header=T,sep=";",dec=".")

# Temp as factor
fec$temp<-as.factor(fec$temp)
tank<-unique(fec$tank)


#--------------------------------------------------
# Sexual maturity

# First sexual maturity (age and length) for each tank
mat<-data.frame(age=as.numeric(),
                tank=as.numeric(),
                food=as.numeric(),
                temp=as.numeric(),
                condition=as.numeric(),
                length_mm=as.numeric(),
                fecundity=as.numeric())
for(i in 1:length(tank)){
  res<-fec[fec$tank==tank[i] & fec$age==min(fec[fec$tank==tank[i],]$age),]
  mat<-rbind(mat,res)
}
mat

# Mean sexual maturity and length at age per tank
df1 <- ddply(mat, .(condition,tank,temp,food), summarise,
             age       = mean(age),
             length_mean    = mean(length_mm))
df1
# Mean sexual maturity and length at age per treatment
df2 <- ddply(df1, .(condition,temp,food), summarise,
             mean_age       = mean(age),
             sd_age         = sd(age),
             mean_length    = mean(length_mean),
             sd_length      = sd(length_mean))
df2


#--------------------------------------------------
# Clutch size


# Mean clutch size per female and day
clutches <- ddply(fec, .(age,condition,tank,temp,food), summarise,
                  n_eggs  =             sum(fecundity),
                  n_fish  =             length(fecundity),
                  mean_clutch_size =  sum(fecundity)/length(fecundity))
clutches

# Linear mixed effect model
M1<-lmer(log(mean_clutch_size)~temp*food+(1|tank),
         na.action=na.omit, data=clutches)
Anova(M1)
summary(M1)
# Model asumptions
res_lme=residuals(M1)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)



#--------------------------------------------------
# Clutch size plot (figure 3 of article)

# Form ggplot
clutches$temp<-str_replace(clutches$temp,"20","20 \u00B0C")
clutches$temp<-as.factor(str_replace(clutches$temp,"30","30 \u00B0C"))
clutches$food<-str_replace(clutches$food,"conti","Continuous")
clutches$food<-str_replace(clutches$food,"inter","Intermittent")

# Plot fecundity
p1<-ggplot(clutches,aes(x=age,y=mean_clutch_size,color=temp,shape=food))+
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(values = c("black","red"))+
  scale_shape_manual(values = c(19,1))+
  labs(y=" Mean clutch size (n eggs)",x="Age (days)",color="Temperature",shape="Food")+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=12,color="black"),
        legend.title = element_text(size=13),
        legend.position="none",
        legend.text = element_text(size=13))
# Boxplot fecundity
p2<-ggplot(clutches,aes(x=condition,y=mean_clutch_size,color=temp,fill=condition))+
  geom_boxplot(aes(alpha=0.7))+
  theme_bw()+
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values=c("#333333","brown1","white","white"))+
  labs(y=" Mean clutch size (n eggs)",x="Treatment",color="Temperature",shape="Food")+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=12,color="black"),
        legend.title = element_text(size=13),
        legend.position="none",
        axis.title.y = element_blank(),
        legend.text = element_text(size=13))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
lgd<-get_legend(p1)

grid.arrange(p1,p2,lgd,ncol=3,nrow=1,widths=c(1.2, 0.8, 0.4))




# Von Bertallanfy posteriors

setwd("C:/Users/sibazin/Desktop/Simon/phd/chapitre_2/article/data/ipm")
cfs_bayes<-read.csv(file = "cfs_bayes.csv",header=T,sep=";",dec=".")





# Mean clutch size per female and day

clutches2 <- ddply(clutches, .(condition,temp,food), summarise,
                   n_eggs  =             sum(n_eggs,na.rm=T),
                   clutch_size= mean(mean_clutch_size,na.rm=T))
clutches2


clutches2$K<-NA
clutches2$Linf<-NA

for(i in 1:nrow(clutches2)){
  clutches2[i,]$K<-cfs_bayes[cfs_bayes$cond==clutches2[i,]$condition & cfs_bayes$coef=="K",]$median
  clutches2[i,]$Linf<-cfs_bayes[cfs_bayes$cond==clutches2[i,]$condition & cfs_bayes$coef=="Linf",]$median
}

clutches2$N_fish<-c(64,84,52,80)
clutches2$n_eggs.F<-clutches2$n_eggs/clutches2$N_fish
