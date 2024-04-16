########## Von Bertalanffy growth curves ########## 
###################################################


#--------------------------------------------------

# This script estimates the Von Bertallanfy growth curves coefficients (K and Linf) 
# for each treatment (inter_20, conti_20, inter_30, conti_30) by a Bayesian model.

#--------------------------------------------------


# Packages
library(R2jags)
library(ggmcmc)
library(ggplot2)
library(plyr)
library(FSA)
library(dplyr)   # for filter(), mutate()
library(nlme)
library(rjags)
library(corrplot)
library(patchwork)

# Data
dta <- read.csv(file = "df_growth.csv",header=T,sep=";",dec=".")

#--------------------------------------------------
# Number of fish measured at different ages (Fig S1 of article)

# N fish measured by age
df1<-ddply(dta, .(age,condition,tank,temp,food), summarise,
           N           = length(na.omit(age)))
df2<-ddply(df1, .(age,condition,temp,food), summarise,
           N           = mean(N))

# Plot number of fish measured by age
ggplot(df2, aes(x=as.factor(age), y=N,color=as.character(temp),shape=food))+
  geom_jitter(position = position_dodge(0.4),size=2.5)+
  theme_bw()+
  scale_color_manual(values = c("black","red"))+
  scale_shape_manual(values = c(1,16))+
  scale_y_continuous(limits = c(0,17))+
  labs(x="Age (days)", y="Mean number of fish measured")+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=13),
        legend.position = "none")


#--------------------------------------------------
# Bayesian model

# Von Bertallanfy model
vb1 <- vbFuns() 

# Multiplicative non linear model
# General model to estimate initial values
auxi <- nls(length_mm ~ vb1(age, Linf, K, t0), data = dta,
            start = list(Linf = 35.8, K = 0.0085, t0 = 2.6))

# Per condition model to estimate gnls initial values
auxi2 <- coef(nlsList(length_mm ~ vb1(age, Linf, K, t0)|condition, data = dta,
                      start = coef(auxi)))

# Per tank estimation
auxi3 <- coef(nlsList(length_mm ~ vb1(age, Linf, K, t0)|tank, data = dta,
                      start = coef(auxi)))

# Experimental design
conds <- unique(dta[, c("tank", "condition")])
conds$gr <- match(conds$condition, c("inter_20", "conti_20", "inter_30", "conti_30"))
conds$C <- 11:1

# Standard deviation approximation
f1 <- function(x) x - mean(x)
conds <- cbind(conds, auxi3[conds$tank,])
conds %>% 
  group_by(condition) %>% 
  mutate_at(5:7,f1) %>% 
  ungroup() %>% 
  summarise_at(5:7, sd) 

# Jags data

# Length at age
dta2 <- dta[, c("length_mm", "age")]
colnames(dta2)<-c("l","age")
# N
N <- nrow(dta2)
# As list
dta2 <- as.list(dta2)
dta2 <- c(dta2, N = N)
# Replicate number
dta2$C <- conds$C[match(dta$tank, conds$tank)]
# Condition
dta2$gr <- conds$gr
dta2$Ci <- conds$C
# Condition
# 1. inter_20
# 2. conti_20
# 3. inter_30
# 4. conti_30

# Names
jags.data <- names(dta2)

# Inits
jags.inits4 <- function() {
  list(K = rnorm(4, c(0.0062279, 0.00742, 0.01, 0.01),
                 0.001),
       t0 = rnorm(4, c(2, 4)),
       Linf = rnorm(4, c(39, 39, 33, 34), 0.9),
       s = 1/3,
       sK = 0.0005,
       sL = 0.6,
       st = 3.1,
       rho1 = -0.9,
       rho2 = 0.85,
       rho3 = -0.9
  )
}

# Param
jags.param <- c("K", "t0", "Linf", 
                "sL", "st", "sK", "tau",
                "rho1", "rho2", "rho3")

load.module("glm")
# Run bayesian model
resto <- jags.parallel(data = dta2, inits = jags.inits4,
                       jags.param, n.iter = 4e5+5e4,
                       n.chains = 5, n.burnin = 5e4, n.thin = 200,
                       model.file = "bayesian_model.txt",
                       DIC = TRUE)

# Gelman diagnostics
cfs <- as.data.frame(resto$BUGSoutput$sims.matrix[, !grepl("deviance", colnames(resto$BUGSoutput$sims.matrix))])
l2 <- (apply(resto$BUGSoutput$sims.array,2, I, simplify = F))
l2 <- lapply(l2, mcmc, 5e4, 4e5, 200)
gelman.diag(mcmc.list(l2))

# Summary outputs
resto$BUGSoutput$summary

# Quantiles
t(apply(cfs, 2, quantile, c(0.025, 0.975)))

# Coefs names
cfsnames<-c("K[1]","K[2]","K[3]","K[4]","Linf[1]","Linf[2]","Linf[3]","Linf[4]","t0[1]","t0[2]","t0[3]","t0[4]")

# Mean and quantiles
df<-data.frame(coef=rep(c("K","Linf","t0"),each=4),
               temp=rep(c("20 \u00B0C","20 \u00B0C","30 \u00B0C","30 \u00B0C"),3),
               food=rep(c("Intermittent","Continuous", "Intermittent","Continuous"),3),
               condition=rep(c("inter_20","conti_20","inter_30","conti_30"),3))
df$mean<-apply(cfs%>%select(c(cfsnames)), 2, mean)
df$median<-apply(cfs%>%select(c(cfsnames)), 2, median)
df$q2.5<-apply(cfs%>%select(c(cfsnames)), 2, quantile, c(0.025))
df$q97.5<-apply(cfs%>%select(c(cfsnames)), 2, quantile, c(0.975))

# Bayesian model diagnostic
diag<-as.mcmc(resto)
s<-ggs(diag)

# Save pdf bayesian model diags
#ggmcmc(s,file="bayesian_model_diags.pdf",
#       plot = c("histogram","density","traceplot","running","compare_partial","autocorrelation","crosscorrelation"))



#--------------------------------------------------
# Fitted Von Bertallanfy growth curve (Fig 2 of article)

# Start a list
res <- vector("list", 4L)

# Only VB parameters
cfs<-cfs%>%select(c(cfsnames))

# Estimated length at age (mean and quantiles) for each condition
# Replace 350 (days) by 700 (days) to obtain Fig S3 of article
for(i in 1:4) {
  cf <- cfs[, i+c(0,4,8)]
  re <- matrix(numeric(), 350, 3)
  for (j in 1:350) {
    le <- cf[,2]*(1-exp((-1*cf[,1])*(j-cf[,3])))
    re[j,] <- c(mean(le), quantile(le, c(0.025, 0.975)))
  }
  re <- cbind(i, 1:350, re)
  colnames(re) <- c("condition","age", "mean", "q2.5", "q97.5")
  res[[i]] <- re
}

res <- as.data.frame(do.call("rbind", res))
res$condition <- as.factor(res$condition)
res$temp <- factor(rep(c("20 \u00B0C", "30 \u00B0C"), each = 2)[res$condition],levels = c("30 \u00B0C", "20 \u00B0C"))
res$food <- factor(rep(c("Intermittent", "Continuous"), 2)[res$condition])

# Check ages when overlap CIs
which(res[res$temp=="20 \u00B0C" & res$food=="Intermittent",]$q97.5>res[res$temp=="20 \u00B0C" & res$food=="Continuous",]$q2.5)

# Renames variables for plot
dta[dta$temp=="30",]$temp<-"30 \u00B0C"
dta[dta$temp=="20",]$temp<-"20 \u00B0C"
dta[dta$food=="conti",]$food<-"Continuous"
dta[dta$food=="inter",]$food<-"Intermittent"

# Plot predicted Von Bertallanfy growth curve
res %>% ggplot(aes(age, mean, color = temp,
                   fill = temp,lty=food,shape=food)) +
  geom_jitter(data=dta,aes(x=age,y=length_mm),size=1.5,position=position_dodge(15))+
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), alpha=0.1) +
  geom_line()+
  scale_x_continuous(breaks=c(30,60,100,150,200,300,350,700))+ #add 700 days for Fig S3
  scale_color_manual(values=c("black","red"))+
  scale_fill_manual(values=c("black","red"))+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_shape_manual(values = c(19,1))+
  theme_bw()+
  labs(x="Age (days)",y="Length (mm)", color="Temperature",fill="Temperature",lty="Food",shape="Food ")+
  theme(axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        #legend.position="none",
        legend.text = element_text(size=12))+
  guides(shape=guide_legend(override.aes = list(size = 2.5)),
         lty=guide_legend(override.aes = list(size = 0.7)))+
  geom_segment(aes(x = 67.3,xend=67.3 ,y = 0, yend = 16.8),size=1,lty=1)+
  geom_segment(aes(x = 60,xend=60 ,y = 0, yend = 17.2),size=1,lty=2)+
  geom_segment(aes(x = 169.7,xend=169.7 ,y = 0, yend = 27.8),size=1,lty=1,col="black")+
  geom_segment(aes(x = 186.5,xend=186.5 ,y = 0, yend = 26.7),size=1,lty=2,col="black")


  

#--------------------------------------------------
# Estimated Von Bertallanfy parameters (Fig S2 of article)

# Plot Von Bertallanfy parameters (median as the standard deviations and correlations distributions are asymmetrics)
ggplot(df,aes(x=condition,y=median,color=temp,shape=food))+
  geom_point()+
  geom_pointrange(aes(ymin=q2.5,ymax=q97.5))+
  theme_bw()+
  facet_wrap(~coef,scale="free")+
  labs(color="Temperature",shape="Food",x="Treatment",y="Median value of parameter")+
  scale_color_manual(values = c("black","red"))+
  scale_shape_manual(values = c(16,1))+
  theme(axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11))


#--------------------------------------------------
# Autocorrelation and parameters correlation

# Density plot of coefs
pdf("to_2.pdf")
for (i in 1:ncol(resto$BUGSoutput$sims.matrix)) {
  bob <- resto$BUGSoutput$sims.array[,,i]
  bob <- as.data.frame(bob)
  bob <- bob %>% gather()
  bob <- bob %>% mutate(key = factor(key))
  # bob <- bob %>% group_by(key) %>% mutate(x = 1:n()) %>% ungroup()
  p1 <- ggplot(bob[,], aes(value, group = key))  + geom_density(alpha = 0.7) +  
    ggtitle(label = colnames(resto$BUGSoutput$sims.matrix)[i])
  plot(p1)
}
dev.off()

# Autocorrelation of parameters
pdf("ACF_to7.pdf")
for (i in 1:ncol(resto$BUGSoutput$sims.matrix)) {
  bob <- resto$BUGSoutput$sims.array[,,i]
  par(mfrow = c(2, 3), mar = c(2,1,3,0.5))
  for (j in 1:5) {
    acf(bob[,j], lag.max = 300, main = colnames(resto$BUGSoutput$sims.matrix)[i])
  }
  mtext(colnames(resto$BUGSoutput$sims.matrix)[i], outer = T, padj = 0.5, line = -1.5)
}
dev.off()


# Data for parameters correlation
tab<-data.frame(condi=rep(c(1:4),each=nrow(cfs)),
                K=NA,Linf=NA,t0=NA)
for(i in 1:4){
  tab[tab$condi==i,]$K<-cfs[,i]
  tab[tab$condi==i,]$Linf<-cfs[,i+4]
  tab[tab$condi==i,]$t0<-cfs[,i+8]
}
# Parameters correlation
l3 <- vector("list", 4L)
for (i in 1:4) {
  cor1 <- cor(tab[tab$condi==i,2:4])
  l3[[i]] <- cor1
  ted <- data.frame(v1 = c("K", "K", "Linf"),
                    v2 = c("t0", "Linf", "t0"),
                    y =  cor1[lower.tri(cor1)])
  bob <- ted %>% mutate(v2 = factor(v2, levels = unique(v2))) %>%
    ggplot(aes(v1, v2, fill = y)) + geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Pearson\nCorrelation") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.title = element_blank())+
    coord_fixed() +
    ggtitle(c("inter_20", "conti_20", "inter_30", "conti_30")[i])
  assign(sprintf("p%d", i), bob)
  rm(bob)
}

# Plot correlation
combined <- (p1 + p2)/(p3+p4) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")

