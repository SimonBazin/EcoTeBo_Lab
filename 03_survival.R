############# Survival curve analyses #############
###################################################


#--------------------------------------------------

# This script analyzes fish survival data by fitting Cox Proportional Hazard Model.

#--------------------------------------------------

# Package
library(survival)
library(survminer)
library(car)
library(stringr)
library(coxme)

# Data
data<-read.csv(file = "df_survival.csv",header=T,sep=";",dec=".")

# Life Data Classification
# http://reliawiki.org/index.php/Life_Data_Classification
# Right censored implies that the event of interest (i.e., the death) is to the right of our data point. 
# If the units were to keep on operating, the failure would occur at some time after our data point (or to the right on the time scale).

# Our data : Right Censored Data (event : 1 = dead, 0 = alive)


#--------------------------------------------------
# Cox Proportional Hazard Model with random effect


# Temp, food and tank as factors
data$temp<-as.factor(data$temp)
data$food<-as.factor(data$food)
data$tank<-as.factor(data$tank)


# Interactive Cox Proportional Hazard Model with random effect
cox.interact<-coxme(Surv(age,event)~temp*food+(1|tank),data=data)
summary(cox.interact)
Anova(cox.interact, test.statistic = "Wald")

# Hazard proportionality
check_PH <- cox.zph(cox.interact)
check_PH 
plot(check_PH)
abline(h = coef(cox.interact)[3], col = "red", lwd = 2)
# Hazard proportionality is ok (p > 0.05 for each var and global)


# Additive Cox Proportional Hazard Model with random effect
cox.additive<-coxme(Surv(age,event)~temp+food+(1|tank),data=data)
summary(cox.additive)
Anova(cox.additive,  test.statistic = "Wald")
# Hazard proportionality (SM 3 of article)
check_PH <- cox.zph(cox.additive)
check_PH 
par(mfrow = c(1, 2))
plot(check_PH,var=1)
abline(h = coef(cox.additive)[1], col = "red", lwd = 2)
plot(check_PH,var=2)
abline(h = coef(cox.additive)[2], col = "red", lwd = 2)
# Hazard proportionality is ok (p > 0.05 for each var and global)

# Best model is additive model
anova(cox.interact,cox.additive)


#--------------------------------------------------
# Fitted survival curves (figure 4 of article)

# Rename variables for plot
data$temp<-str_replace(data$temp,"20","20 \u00B0C")
data$temp<-str_replace(data$temp,"30","30 \u00B0C")
data$food<-str_replace(data$food,"conti","Continuous")
data$food<-str_replace(data$food,"inter","Intermittent")

# Kaplan Meier survival curves
km.fit<-survfit(Surv(age,event)~food+temp,data=data)

# Plot Kaplan Meier survival curves
ggsurv<-ggsurvplot(km.fit,conf.int = TRUE,data=data,
                   col="temp",palette=c("black","red"),
                   ggtheme=theme_bw(),legend="right",legend.title="Temperature",
                   conf.int.alpha =0.2,
                   ylim = c(0.20,1),
                   xlim = c(60,350))
ggsurv$plot +
  facet_wrap(food~.)+
  #theme (legend.position = "none")+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  scale_x_continuous(breaks = c(60,100,200,300))+
  labs(palette="Temperature",x="Age (days)")
