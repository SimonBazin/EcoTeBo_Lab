########## IPM model ########## 
###############################

rm(list=ls()) 

#--------------------------------------------------

# This script is used to perform a sensitivity analysis of IPM model for each population

#--------------------------------------------------

# Packages
library(ipmr)
library(Rage)
library(ggplot2)
library(dplyr)
library(scales)
library(stringr)

# Data
# Fecundity
fec<-read.csv(file = "df_fecundity.csv",header=T,sep=";",dec=".")
# Reproduction
reproduction<-read.csv(file = "df_reproduction.csv",header=T,sep=";",dec=".")
# Survival
survival<-read.csv(file = "df_survival_probabilities.csv",header=T,sep=";",dec=".")
# Von Bertallanfy posteriors
cfs_bayes<-read.csv(file = "cfs_bayes.csv",header=T,sep=";",dec=".")


#################################### continuous 30 #########################################
############################################################################################


# Select only one condition: conti_30
fec_i<-fec[which(fec$cond=="conti_30"),]  
fec_i$log_size<-log(fec_i$length_mm)
repro_i<-reproduction[which(reproduction$cond=="conti_30"),]
surv_i<-survival[which(survival$cond=="conti_30"),] 


################# plot biological rate as a function of log body size #############################

# Size dependent vital rates

par(mfrow = c(2,2))

# Survival analyses
s_model_1  <- glm(surv ~ log_size, data = surv_i, family = binomial("logit"), weights = total)
summary(s_model_1)
(s_pars <- coef(s_model_1))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(s_model_1, Predicted_data, type="response")

# Plot survival
plot((surv) ~ log_size,
     data = surv_i,
     ylim= c(0,1),
     xlim=c(2.5,4),
     xlab = "Log Size (t)",
     ylab = "Survival (t + 1)")
lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")

# Growth

# Linf & K from the Von Bertalanffy model
# posterior distribution
# 1. inter_20
# 2. conti_20
# 3. inter_30
# 4. conti_30

# Load data
cfs_bayes<-as.data.frame(cfs_bayes)

# VB parameters for conti_30 condition
cfs <- cfs_bayes[, grep("\\[4\\]", colnames(cfs_bayes))]

Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)

# Predicted size at time t
vb1 <-function(age,Linf, K, t0){Linf*(1-exp(-K*(age-t0)))}
lt <- vb1(time, Linf, K, t0)
plot(lt~time)

# Remove fish with lower size that observed in experiment
lt <- lt[which (lt>8.359746)]
dta <- data.frame(y = log(lt[-1]), lt = log(lt[-length(lt)]))
g_model <- lm(y ~ lt, data = dta)
(g_pars <- c(coef(g_model), sd(resid(g_model))))

summary(g_model)
log((Linf*(1-exp(-K)))) #intercept
log(exp(-K)) ##slope

(g_pars <- c(coef(g_model), sd(resid(g_model))))
g_m<-function(bl){
  g_pars[1] + g_pars[2] * bl
}
body <-seq(0,50)

plot(y ~ lt,
     data = dta,
     xlab = "Size (t)",
     ylab = " Size (t + 1)")
lines(g_m(body)~body,lwd=2, col="blue")


# predict size at time t+1 from size at t
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
str(vbst)
# size at time t
vbst <- vbst[rowMeans(vbst)>8.359746, ]
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd size at t+1 
sdsl1 <- apply(log(vbst1), 1, sd)
mean(sdsl1)
# standard deviation as a function of Lt
sdg <- approxfun(rowMeans(log(vbst1)), sdsl1)
# sdg : function which allows to calculate the standard deviation of the predicted value of lt+1 for a given value of Lt
pred <- sdg(rowMeans(log(vbst1)))
plot(pred, sdsl1)
plot(rowMeans(vbst1),sdsl1)
plot(rowMeans(vbst1),pred)


# Reproduction

f_p_model <- glm(repro ~ log_size , data = repro_i, family = binomial("logit"))
summary(f_p_model)
(f_p_pars <- coef(f_p_model))
 
Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(f_p_model, Predicted_data, type="response")
 
plot(jitter(repro, amount = 0.05) ~ log_size,
     data = repro_i,
     ylim= c(0,1),
     xlim=c(-1,5),
     xlab = "Log Size (t)",
     ylab = "Reproductive (t)")
 
lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")



# Fecundity

f_n_model <- glm(fecundity ~ log_size , data = fec_i, family = poisson())
summary(f_n_model)
f_n_pars <- coef(f_n_model)

Predicted_data <- data.frame(log_size=seq(2.5, 3.8,len=500))

# Fill predicted values using regression model
Predicted_data$fec = predict(f_n_model, Predicted_data, type="response")


plot(fecundity ~ log_size, 
     data = fec_i,
     ylim= c(0,30),
     xlim=c(2.6,3.8),
     xlab = "Size (t)",
     ylab = "# of eggs (t)")
lines(fec ~ log_size, Predicted_data, lwd=2, col="blue")




# Size independent vital rates

# Mean size at t0 in log scale and sd
(f_d_pars  <- c(1.4926754,0.1129072))

# Hatching rate * survival t0
f_g_pars<-0.468*0.553

# List pars
all_pars <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "alpha_G" = g_pars[1], #intercept growth
  "beta_G" = g_pars[2], #slope growth
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars              #survival until from egg to day 30
)

# Init ipm
carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

# Define kernel
carpobrotus_ipm <- define_kernel( #subkernel for existing fish
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G,                                   
  s         = boot::inv.logit(alpha_s + beta_s * z_1),
  g_mu = log((Linf*(1-exp(-K)) + exp(-K)*exp(z_1))),
  G         = dnorm(z_2, g_mu, sigma_G), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( #subkernel for new generation
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,        
  f_d       = dnorm(z_2,mu_f_d,sigma_f_d),
  f_p       = plogis(alpha_f_p + beta_f_p_1* z_1),
  f_n       = exp(alpha_f_n + beta_f_n_1* z_1),
  data_list = all_pars,
  states    = list(c("z")),
  uses_par_sets = FALSE,
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "f_d")
)

carpobrotus_ipm <- define_impl(
  proto_ipm = carpobrotus_ipm,
  kernel_impl_list = list(
    
    P = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    ),
    F = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    )
  )
)



# Integration domain
L<-repro_i$log_size[repro_i$age == 30]
U<-max(fec_i$log_size,na.rm=TRUE) # larger size in the fecundity dataset
m <- 400

carpobrotus_ipm <- define_domains(
  proto_ipm = carpobrotus_ipm,
  z = c(L, U, m)
)

carpobrotus_ipm <- define_pop_state(
  proto_ipm = carpobrotus_ipm,
  n_z       = rep(1/m, m)
)

carpobrotus_ipm_unif <- make_ipm(
  proto_ipm       = carpobrotus_ipm,
  return_all_envs = TRUE,
  iterations      = 1500
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 10:1500)

# Generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)
lambda_obs

# Sensitivity analyses

pars_in_base <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,             #survival until from egg to day 30
  para_sensi = "baseline"
  )

s1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]+s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,          #survival until from egg to day 30
  para_sensi = "survival_+1"
)

s1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]-s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,              #survival until from egg to day 30
  para_sensi = "survival_-1"
)

g1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf+Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "Linf_+1"
)

g1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf-Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                ##survival until from egg to day 30
  para_sensi = "Linf_-1"
)

gi1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K+K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_+1"
)

gi1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K-K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_-1"
)

fp1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]+f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_+1"
)

fp1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]-f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_-1"
)

fn1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro,
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]+f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_+1"
)

fn1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]-f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_-1"
)

mu1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1]+f_d_pars[1]*0.01,        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_+1"
)

mu1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1]-f_d_pars[1]*0.01,   #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_-1"
)

fg1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
   f_g = f_g_pars+f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_+1"
)

fg1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars-f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_-1"
)

pars_in<-rbind(pars_in_base,s1p,s1m,g1p,g1m,gi1p,gi1m,fp1p,fp1m,fn1p,fn1m,mu1p,mu1m,fg1p,fg1m)

use_proto <- carpobrotus_ipm_unif$proto_ipm

Result_CONTI_30<-data.frame (condition= numeric(0),lambda=numeric(0),gen_time=numeric(0),sensi=numeric(0))


for(i in 1:nrow(pars_in)) {
  
  all_pars <- list(
    "alpha_s" = pars_in$alpha_s[i], #intercept survival
    "beta_s" = pars_in$beta_s[i], #slope survival
    "sigma_G" = pars_in$sigma_G[i], #sd resid growth
    "Linf" = pars_in$Linf[i], #slope survival
    "K" = pars_in$K[i], #sd resid growth
    "alpha_f_p" = pars_in$alpha_f_p[i], #intercept proba repro
    "beta_f_p_1" = pars_in$beta_f_p_1[i], #slope proba repro
    "alpha_f_n" = pars_in$alpha_f_n[i], #intercept fecundity
    "beta_f_n_1" = pars_in$beta_f_n_1[i], #slope fecundity
    "mu_f_d" = pars_in$mu_f_d[i],        #mean size
    "sigma_f_d" = pars_in$sigma_f_d[i],    #size sd
    "f_g" = pars_in$f_g[i]              #survival until from egg to day 30
  )
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  
  parameters(use_proto) <- all_pars
  
  boot_ipm <- use_proto %>%
    make_ipm(iterate = TRUE,
             iterations = 2000)
  
  #all_lambdas[i] <- lambda(boot_ipm)
  Tempfile <- data.frame(condition= "CONTI_30",lambda = lambda(boot_ipm),
                         gen_time = gen_time(boot_ipm$sub_kernels$P,boot_ipm$sub_kernels$F),
                         sensi=pars_in$para_sensi[i])
  Result_CONTI_30 <- rbind(Result_CONTI_30,Tempfile)
  
  
}

Result_CONTI_30$lambda_ratio<-Result_CONTI_30$lambda/Result_CONTI_30$lambda[1]
Result_CONTI_30$gen_time_ratio<-Result_CONTI_30$gen_time/Result_CONTI_30$gen_time[1]

write.table(Result_CONTI_30, "conti_30_sensitivity.csv", row.names=FALSE,sep=";")





#################################### conti 20 #########################################
############################################################################################



rm(list=setdiff(ls(), c("cfs_bayes","fec","reproduction","survival"))) 

# Select only one condition: conti_20
fec_i<-fec[which(fec$cond=="conti_20"),]  
fec_i$log_size<-log(fec_i$length_mm)
repro_i<-reproduction[which(reproduction$cond=="conti_20"),]
surv_i<-survival[which(survival$cond=="conti_20"),] 


################# plot biological rate as a function of log body size #############################

# Size dependent vital rates

par(mfrow = c(2,2))

# Survival analyses
s_model_1  <- glm(surv ~ log_size, data = surv_i, family = binomial("logit"), weights = total)
summary(s_model_1)
(s_pars <- coef(s_model_1))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(s_model_1, Predicted_data, type="response")

# Plot survival
plot((surv) ~ log_size,
     data = surv_i,
     ylim= c(0,1),
     xlim=c(2.5,4),
     xlab = "Log Size (t)",
     ylab = "Survival (t + 1)")
lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")

# Growth

# Load data
cfs_bayes<-as.data.frame(cfs_bayes)

# VB parameters for conti_20 condition
cfs <- cfs_bayes[, grep("\\[2\\]", colnames(cfs_bayes))]

Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)

# Predicted size at time t
vb1 <-function(age,Linf, K, t0){Linf*(1-exp(-K*(age-t0)))}
lt <- vb1(time, Linf, K, t0)
plot(lt~time)

# Remove fish with lower size that observed in experiment
lt <- lt[which (lt>7.435182)]
dta <- data.frame(y = log(lt[-1]), lt = log(lt[-length(lt)]))
g_model <- lm(y ~ lt, data = dta)
(g_pars <- c(coef(g_model), sd(resid(g_model))))

summary(g_model)
log((Linf*(1-exp(-K)))) #intercept
log(exp(-K)) ##slope

(g_pars <- c(coef(g_model), sd(resid(g_model))))
g_m<-function(bl){
  g_pars[1] + g_pars[2] * bl
}
body <-seq(0,50)

plot(y ~ lt,
     data = dta,
     xlab = "Size (t)",
     ylab = " Size (t + 1)")
lines(g_m(body)~body,lwd=2, col="blue")


# predict size at time t+1 from size at t
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
str(vbst)
# size at time t
vbst <- vbst[rowMeans(vbst)>7.435182, ]
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd size at t+1 
sdsl1 <- apply(log(vbst1), 1, sd)
mean(sdsl1)
# standard deviation as a function of Lt
sdg <- approxfun(rowMeans(log(vbst1)), sdsl1)
# sdg : function which allows to calculate the standard deviation of the predicted value of lt+1 for a given value of Lt
pred <- sdg(rowMeans(log(vbst1)))
plot(pred, sdsl1)
plot(rowMeans(vbst1),sdsl1)
plot(rowMeans(vbst1),pred)


# Reproduction

head(repro_i)
f_p_model <- glm(repro ~ log_size , data = repro_i, family = binomial("logit"))
summary(f_p_model)
(f_p_pars <- coef(f_p_model))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(f_p_model, Predicted_data, type="response")


plot(jitter(repro, amount = 0.05) ~ log_size,
     data = repro_i,
     ylim= c(0,1),
     xlim=c(-1,5),
     xlab = "Log Size (t)",
     ylab = "Reproductive (t)")

lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")


# Fecundity

f_n_model <- glm(fecundity ~ log_size , data = fec_i, family = poisson())
summary(f_n_model)
f_n_pars <- coef(f_n_model)

Predicted_data <- data.frame(log_size=seq(2.5, 3.8,len=500))

# Fill predicted values using regression model
Predicted_data$fec = predict(f_n_model, Predicted_data, type="response")


plot(fecundity ~ log_size, 
     data = fec_i,
     ylim= c(0,30),
     xlim=c(2.6,3.8),
     xlab = "Size (t)",
     ylab = "# of eggs (t)")
lines(fec ~ log_size, Predicted_data, lwd=2, col="blue")


# Size independent vital rates

# Mean size at t0 in log scale and sd
(f_d_pars  <- c(1.58770120,0.05455395))

# Hatching rate * survival t0
f_g_pars<-0.507*0.489



# List pars
all_pars <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "alpha_G" = g_pars[1], #intercept growth
  "beta_G" = g_pars[2], #slope growth
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars              #survival until from egg to day 30
)

# Init ipm
carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

# Define kernel
carpobrotus_ipm <- define_kernel( #subkernel for existing fish
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G,                                   
  s         = boot::inv.logit(alpha_s + beta_s * z_1),
  g_mu = log((Linf*(1-exp(-K)) + exp(-K)*exp(z_1))),
  G         = dnorm(z_2, g_mu, sigma_G), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( #subkernel for new generation
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,        
  f_d       = dnorm(z_2,mu_f_d,sigma_f_d),
  f_p       = plogis(alpha_f_p + beta_f_p_1* z_1),
  f_n       = exp(alpha_f_n + beta_f_n_1* z_1),
  data_list = all_pars,
  states    = list(c("z")),
  uses_par_sets = FALSE,
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "f_d")
)

carpobrotus_ipm <- define_impl(
  proto_ipm = carpobrotus_ipm,
  kernel_impl_list = list(
    
    P = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    ),
    F = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    )
  )
)



# Integration domain
L<-repro_i$log_size[repro_i$age == 30]
U<-max(fec_i$log_size,na.rm=TRUE) # larger size in the fecundity dataset
m <- 400

carpobrotus_ipm <- define_domains(
  proto_ipm = carpobrotus_ipm,
  z = c(L, U, m)
)

carpobrotus_ipm <- define_pop_state(
  proto_ipm = carpobrotus_ipm,
  n_z       = rep(1/m, m)
)

carpobrotus_ipm_unif <- make_ipm(
  proto_ipm       = carpobrotus_ipm,
  return_all_envs = TRUE,
  iterations      = 10000
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 10:10000)

# Generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)
lambda_obs

# Sensitivity analyses

pars_in_base <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,             #survival until from egg to day 30
  para_sensi = "baseline"
)

s1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]+s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,          #survival until from egg to day 30
  para_sensi = "survival_+1"
)

s1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]-s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,              #survival until from egg to day 30
  para_sensi = "survival_-1"
)

g1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf+Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "Linf_+1"
)

g1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf-Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                ##survival until from egg to day 30
  para_sensi = "Linf_-1"
)

gi1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K+K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_+1"
)

gi1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K-K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_-1"
)

fp1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]+f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_+1"
)

fp1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]-f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_-1"
)

fn1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro,
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]+f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_+1"
)

fn1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]-f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_-1"
)

mu1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1]+f_d_pars[1]*0.01,        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_+1"
)

mu1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1]-f_d_pars[1]*0.01,   #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_-1"
)

fg1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars+f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_+1"
)

fg1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars-f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_-1"
)

pars_in<-rbind(pars_in_base,s1p,s1m,g1p,g1m,gi1p,gi1m,fp1p,fp1m,fn1p,fn1m,mu1p,mu1m,fg1p,fg1m)

use_proto <- carpobrotus_ipm_unif$proto_ipm

Result_CONTI_20<-data.frame (condition= numeric(0),lambda=numeric(0),gen_time=numeric(0),sensi=numeric(0))


for(i in 1:nrow(pars_in)) {
  
  all_pars <- list(
    "alpha_s" = pars_in$alpha_s[i], #intercept survival
    "beta_s" = pars_in$beta_s[i], #slope survival
    "sigma_G" = pars_in$sigma_G[i], #sd resid growth
    "Linf" = pars_in$Linf[i], #slope survival
    "K" = pars_in$K[i], #sd resid growth
    "alpha_f_p" = pars_in$alpha_f_p[i], #intercept proba repro
    "beta_f_p_1" = pars_in$beta_f_p_1[i], #slope proba repro
    "alpha_f_n" = pars_in$alpha_f_n[i], #intercept fecundity
    "beta_f_n_1" = pars_in$beta_f_n_1[i], #slope fecundity
    "mu_f_d" = pars_in$mu_f_d[i],        #mean size
    "sigma_f_d" = pars_in$sigma_f_d[i],    #size sd
    "f_g" = pars_in$f_g[i]              #survival until from egg to day 30
  )
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  
  parameters(use_proto) <- all_pars
  
  boot_ipm <- use_proto %>%
    make_ipm(iterate = TRUE,
             iterations = 11000)
  
  #all_lambdas[i] <- lambda(boot_ipm)
  Tempfile <- data.frame(condition= "CONTI_20",lambda = lambda(boot_ipm),
                         gen_time = gen_time(boot_ipm$sub_kernels$P,boot_ipm$sub_kernels$F),
                         sensi=pars_in$para_sensi[i])
  Result_CONTI_20 <- rbind(Result_CONTI_20,Tempfile)
  
  
}

Result_CONTI_20$lambda_ratio<-Result_CONTI_20$lambda/Result_CONTI_20$lambda[1]
Result_CONTI_20$gen_time_ratio<-Result_CONTI_20$gen_time/Result_CONTI_20$gen_time[1]

write.table(Result_CONTI_20, "conti_20_sensitivity.csv", row.names=FALSE,sep=";")



#################################### intermittent 20 #########################################
#############################################################################################



rm(list=setdiff(ls(), c("cfs_bayes","fec","reproduction","survival"))) 

# Select only one condition: res_20
fec_i<-fec[which(fec$cond=="inter_20"),]  
fec_i$log_size<-log(fec_i$length_mm)
repro_i<-reproduction[which(reproduction$cond=="inter_20"),]
surv_i<-survival[which(survival$cond=="inter_20"),] 


################# plot biological rate as a function of log body size #############################

# Size dependent vital rates

par(mfrow = c(2,2))

# Survival analyses
s_model_1  <- glm(surv ~ log_size, data = surv_i, family = binomial("logit"), weights = total)
summary(s_model_1)
(s_pars <- coef(s_model_1))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(s_model_1, Predicted_data, type="response")

# Plot survival
plot((surv) ~ log_size,
     data = surv_i,
     ylim= c(0,1),
     xlim=c(2.5,4),
     xlab = "Log Size (t)",
     ylab = "Survival (t + 1)")
lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")

# Growth

# Load data
cfs_bayes<-as.data.frame(cfs_bayes)

# VB parameters for res_20 condition
cfs <- cfs_bayes[, grep("\\[1\\]", colnames(cfs_bayes))]

Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)

# Predicted size at time t
vb1 <-function(age,Linf, K, t0){Linf*(1-exp(-K*(age-t0)))}
lt <- vb1(time, Linf, K, t0)
plot(lt~time)

# Remove fish with lower size that observed in experiment
lt <- lt[which (lt>5.823352)]
dta <- data.frame(y = log(lt[-1]), lt = log(lt[-length(lt)]))
g_model <- lm(y ~ lt, data = dta)
(g_pars <- c(coef(g_model), sd(resid(g_model))))

summary(g_model)
log((Linf*(1-exp(-K)))) #intercept
log(exp(-K)) ##slope

(g_pars <- c(coef(g_model), sd(resid(g_model))))
g_m<-function(bl){
  g_pars[1] + g_pars[2] * bl
}
body <-seq(0,50)

plot(y ~ lt,
     data = dta,
     xlab = "Size (t)",
     ylab = " Size (t + 1)")
lines(g_m(body)~body,lwd=2, col="blue")


# predict size at time t+1 from size at t
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
str(vbst)
# size at time t
vbst <- vbst[rowMeans(vbst)>5.823352, ]
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd size at t+1 
sdsl1 <- apply(log(vbst1), 1, sd)
mean(sdsl1)
# standard deviation as a function of Lt
sdg <- approxfun(rowMeans(log(vbst1)), sdsl1)
# sdg : function which allows to calculate the standard deviation of the predicted value of lt+1 for a given value of Lt
pred <- sdg(rowMeans(log(vbst1)))
plot(pred, sdsl1)
plot(rowMeans(vbst1),sdsl1)
plot(rowMeans(vbst1),pred)


# Reproduction

head(repro_i)
f_p_model <- glm(repro ~ log_size , data = repro_i, family = binomial("logit"))
summary(f_p_model)
(f_p_pars <- coef(f_p_model))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(f_p_model, Predicted_data, type="response")


plot(jitter(repro, amount = 0.05) ~ log_size,
     data = repro_i,
     ylim= c(0,1),
     xlim=c(-1,5),
     xlab = "Log Size (t)",
     ylab = "Reproductive (t)")

lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")


# Fecundity

f_n_model <- glm(fecundity ~ log_size , data = fec_i, family = poisson())
summary(f_n_model)
f_n_pars <- coef(f_n_model)

Predicted_data <- data.frame(log_size=seq(2.5, 3.8,len=500))

# Fill predicted values using regression model
Predicted_data$fec = predict(f_n_model, Predicted_data, type="response")


plot(fecundity ~ log_size, 
     data = fec_i,
     ylim= c(0,30),
     xlim=c(2.6,3.8),
     xlab = "Size (t)",
     ylab = "# of eggs (t)")
lines(fec ~ log_size, Predicted_data, lwd=2, col="blue")


# Size independent vital rates

# Mean size at t0 in log scale and sd
(f_d_pars  <- c(1.58770120,0.05455395))

# Hatching rate * survival t0
f_g_pars<-0.486*0.489



# List pars
all_pars <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "alpha_G" = g_pars[1], #intercept growth
  "beta_G" = g_pars[2], #slope growth
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars              #survival until from egg to day 30
)

# Init ipm
carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

# Define kernel
carpobrotus_ipm <- define_kernel( #subkernel for existing fish
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G,                                   
  s         = boot::inv.logit(alpha_s + beta_s * z_1),
  g_mu = log((Linf*(1-exp(-K)) + exp(-K)*exp(z_1))),
  G         = dnorm(z_2, g_mu, sigma_G), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( #subkernel for new generation
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,        
  f_d       = dnorm(z_2,mu_f_d,sigma_f_d),
  f_p       = plogis(alpha_f_p + beta_f_p_1* z_1),
  f_n       = exp(alpha_f_n + beta_f_n_1* z_1),
  data_list = all_pars,
  states    = list(c("z")),
  uses_par_sets = FALSE,
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "f_d")
)

carpobrotus_ipm <- define_impl(
  proto_ipm = carpobrotus_ipm,
  kernel_impl_list = list(
    
    P = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    ),
    F = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    )
  )
)



# Integration domain
L<-repro_i$log_size[repro_i$age == 30]
U<-max(fec_i$log_size,na.rm=TRUE) # larger size in the fecundity dataset
m <- 400

carpobrotus_ipm <- define_domains(
  proto_ipm = carpobrotus_ipm,
  z = c(L, U, m)
)

carpobrotus_ipm <- define_pop_state(
  proto_ipm = carpobrotus_ipm,
  n_z       = rep(1/m, m)
)

carpobrotus_ipm_unif <- make_ipm(
  proto_ipm       = carpobrotus_ipm,
  return_all_envs = TRUE,
  iterations      = 4000
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 10:4000)

# Generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)
lambda_obs

# Sensitivity analyses

pars_in_base <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,             #survival until from egg to day 30
  para_sensi = "baseline"
)

s1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]+s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,          #survival until from egg to day 30
  para_sensi = "survival_+1"
)

s1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]-s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,              #survival until from egg to day 30
  para_sensi = "survival_-1"
)

g1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf+Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "Linf_+1"
)

g1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf-Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                ##survival until from egg to day 30
  para_sensi = "Linf_-1"
)

gi1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K+K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_+1"
)

gi1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K-K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_-1"
)

fp1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]+f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_+1"
)

fp1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]-f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_-1"
)

fn1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro,
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]+f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_+1"
)

fn1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]-f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_-1"
)

mu1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1]+f_d_pars[1]*0.01,        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_+1"
)

mu1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1]-f_d_pars[1]*0.01,   #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_-1"
)

fg1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars+f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_+1"
)

fg1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars-f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_-1"
)

pars_in<-rbind(pars_in_base,s1p,s1m,g1p,g1m,gi1p,gi1m,fp1p,fp1m,fn1p,fn1m,mu1p,mu1m,fg1p,fg1m)

use_proto <- carpobrotus_ipm_unif$proto_ipm

Result_RES_20<-data.frame (condition= numeric(0),lambda=numeric(0),gen_time=numeric(0),sensi=numeric(0))


for(i in 1:nrow(pars_in)) {
  
  all_pars <- list(
    "alpha_s" = pars_in$alpha_s[i], #intercept survival
    "beta_s" = pars_in$beta_s[i], #slope survival
    "sigma_G" = pars_in$sigma_G[i], #sd resid growth
    "Linf" = pars_in$Linf[i], #slope survival
    "K" = pars_in$K[i], #sd resid growth
    "alpha_f_p" = pars_in$alpha_f_p[i], #intercept proba repro
    "beta_f_p_1" = pars_in$beta_f_p_1[i], #slope proba repro
    "alpha_f_n" = pars_in$alpha_f_n[i], #intercept fecundity
    "beta_f_n_1" = pars_in$beta_f_n_1[i], #slope fecundity
    "mu_f_d" = pars_in$mu_f_d[i],        #mean size
    "sigma_f_d" = pars_in$sigma_f_d[i],    #size sd
    "f_g" = pars_in$f_g[i]              #survival until from egg to day 30
  )
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  
  parameters(use_proto) <- all_pars
  
  boot_ipm <- use_proto %>%
    make_ipm(iterate = TRUE,
             iterations = 14000)
  
  #all_lambdas[i] <- lambda(boot_ipm)
  Tempfile <- data.frame(condition= "INTER_20",lambda = lambda(boot_ipm),
                         gen_time = gen_time(boot_ipm$sub_kernels$P,boot_ipm$sub_kernels$F),
                         sensi=pars_in$para_sensi[i])
  Result_INTER_20 <- rbind(Result_INTER_20,Tempfile)
  
  
}

Result_INTER_20$lambda_ratio<-Result_INTER_20$lambda/Result_INTER_20$lambda[1]
Result_INTER_20$gen_time_ratio<-Result_INTER_20$gen_time/Result_INTER_20$gen_time[1]

write.table(Result_INTER_20, "inter_20_sensitivity.csv", row.names=FALSE,sep=";")




#################################### intermittent 30 #########################################
#############################################################################################



rm(list=setdiff(ls(), c("cfs_bayes","fec","reproduction","survival"))) 

# Select only one condition: inter_30
fec_i<-fec[which(fec$cond=="inter_30"),]  
fec_i$log_size<-log(fec_i$length_mm)
repro_i<-reproduction[which(reproduction$cond=="inter_30"),]
surv_i<-survival[which(survival$cond=="inter_30"),] 


################# plot biological rate as a function of log body size #############################

# Size dependent vital rates

par(mfrow = c(2,2))

# Survival analyses
s_model_1  <- glm(surv ~ log_size, data = surv_i, family = binomial("logit"), weights = total)
summary(s_model_1)
(s_pars <- coef(s_model_1))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(s_model_1, Predicted_data, type="response")

# Plot survival
plot((surv) ~ log_size,
     data = surv_i,
     ylim= c(0,1),
     xlim=c(2.5,4),
     xlab = "Log Size (t)",
     ylab = "Survival (t + 1)")
lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")

# Growth

# Load data
cfs_bayes<-as.data.frame(cfs_bayes)

# VB parameters for res_30 condition
cfs <- cfs_bayes[, grep("\\[3\\]", colnames(cfs_bayes))]

Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)

# Predicted size at time t
vb1 <-function(age,Linf, K, t0){Linf*(1-exp(-K*(age-t0)))}
lt <- vb1(time, Linf, K, t0)
plot(lt~time)

# Remove fish with lower size that observed in experiment
lt <- lt[which (lt>7.597057)]
dta <- data.frame(y = log(lt[-1]), lt = log(lt[-length(lt)]))
g_model <- lm(y ~ lt, data = dta)
(g_pars <- c(coef(g_model), sd(resid(g_model))))

summary(g_model)
log((Linf*(1-exp(-K)))) #intercept
log(exp(-K)) #slope

(g_pars <- c(coef(g_model), sd(resid(g_model))))
g_m<-function(bl){
  g_pars[1] + g_pars[2] * bl
}
body <-seq(0,50)

plot(y ~ lt,
     data = dta,
     xlab = "Size (t)",
     ylab = " Size (t + 1)")
lines(g_m(body)~body,lwd=2, col="blue")


# predict size at time t+1 from size at t
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
str(vbst)
# size at time t
vbst <- vbst[rowMeans(vbst)>7.597057, ]
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd size at t+1 
sdsl1 <- apply(log(vbst1), 1, sd)
mean(sdsl1)
# standard deviation as a function of Lt
sdg <- approxfun(rowMeans(log(vbst1)), sdsl1)
# sdg : function which allows to calculate the standard deviation of the predicted value of lt+1 for a given value of Lt
pred <- sdg(rowMeans(log(vbst1)))
plot(pred, sdsl1)
plot(rowMeans(vbst1),sdsl1)
plot(rowMeans(vbst1),pred)


# Reproduction

head(repro_i)
f_p_model <- glm(repro ~ log_size , data = repro_i, family = binomial("logit"))
summary(f_p_model)
(f_p_pars <- coef(f_p_model))

Predicted_data <- data.frame(log_size=seq(min(surv_i$log_size-3), max(surv_i$log_size+2),len=500))

# Fill predicted values using regression model
Predicted_data$var1 = predict(f_p_model, Predicted_data, type="response")


plot(jitter(repro, amount = 0.05) ~ log_size,
     data = repro_i,
     ylim= c(0,1),
     xlim=c(-1,5),
     xlab = "Log Size (t)",
     ylab = "Reproductive (t)")

lines(var1 ~ log_size, Predicted_data, lwd=2, col="blue")


# Fecundity

f_n_model <- glm(fecundity ~ log_size , data = fec_i, family = poisson())
summary(f_n_model)
f_n_pars <- coef(f_n_model)

Predicted_data <- data.frame(log_size=seq(2.5, 3.8,len=500))

# Fill predicted values using regression model
Predicted_data$fec = predict(f_n_model, Predicted_data, type="response")


plot(fecundity ~ log_size, 
     data = fec_i,
     ylim= c(0,30),
     xlim=c(2.6,3.8),
     xlab = "Size (t)",
     ylab = "# of eggs (t)")
lines(fec ~ log_size, Predicted_data, lwd=2, col="blue")


# Size independent vital rates

# Mean size at t0 in log scale and sd
(f_d_pars  <- c(1.4926754,0.1129072))

# Hatching rate * survival t0
f_g_pars<-0.364*0.553



# List pars
all_pars <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "alpha_G" = g_pars[1], #intercept growth
  "beta_G" = g_pars[2], #slope growth
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars              #survival until from egg to day 30
)

# Init ipm
carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

# Define kernel
carpobrotus_ipm <- define_kernel( #subkernel for existing fish
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G,                                   
  s         = boot::inv.logit(alpha_s + beta_s * z_1),
  g_mu = log((Linf*(1-exp(-K)) + exp(-K)*exp(z_1))),
  G         = dnorm(z_2, g_mu, sigma_G), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( #subkernel for new generation
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,        
  f_d       = dnorm(z_2,mu_f_d,sigma_f_d),
  f_p       = plogis(alpha_f_p + beta_f_p_1* z_1),
  f_n       = exp(alpha_f_n + beta_f_n_1* z_1),
  data_list = all_pars,
  states    = list(c("z")),
  uses_par_sets = FALSE,
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "f_d")
)

carpobrotus_ipm <- define_impl(
  proto_ipm = carpobrotus_ipm,
  kernel_impl_list = list(
    
    P = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    ),
    F = list(
      int_rule    = "midpoint",
      state_start = "z",
      state_end   = "z"
    )
  )
)



# Integration domain
L<-repro_i$log_size[repro_i$age == 30]
U<-max(fec_i$log_size,na.rm=TRUE) # larger size in the fecundity dataset
m <- 400

carpobrotus_ipm <- define_domains(
  proto_ipm = carpobrotus_ipm,
  z = c(L, U, m)
)

carpobrotus_ipm <- define_pop_state(
  proto_ipm = carpobrotus_ipm,
  n_z       = rep(1/m, m)
)

carpobrotus_ipm_unif <- make_ipm(
  proto_ipm       = carpobrotus_ipm,
  return_all_envs = TRUE,
  iterations      = 1500
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 10:1500)

# Generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)
lambda_obs


# Sensitivity analyses

pars_in_base <- list(
  "alpha_s" = s_pars[1], #intercept survival
  "beta_s" = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,             #survival until from egg to day 30
  para_sensi = "baseline"
)

s1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]+s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,          #survival until from egg to day 30
  para_sensi = "survival_+1"
)

s1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2]-s_pars[2]*0.01, #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,              #survival until from egg to day 30
  para_sensi = "survival_-1"
)

g1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf+Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "Linf_+1"
)

g1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf-Linf*0.01,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                ##survival until from egg to day 30
  para_sensi = "Linf_-1"
)

gi1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K+K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_+1"
)

gi1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], ##slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K-K*0.01,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "k_-1"
)

fp1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]+f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_+1"
)

fp1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2]-f_p_pars[2]*0.01, #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "repro_proba_-1"
)

fn1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro,
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]+f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_+1"
)

fn1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2]-f_n_pars[2]*0.01,  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "fecundity_-1"
)

mu1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1]+f_d_pars[1]*0.01,        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_+1"
)

mu1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2],#slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2], #slope fecundity
  "mu_f_d" = f_d_pars[1]-f_d_pars[1]*0.01,   #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  "f_g" = f_g_pars ,                #survival until from egg to day 30
  para_sensi = "size_t0_-1"
)

fg1p <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars+f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_+1"
)

fg1m <- data.frame(
  alpha_s = s_pars[1], #intercept survival
  beta_s = s_pars[2], #slope survival
  "sigma_G" =  mean(sdsl1), #sd resid growth
  "Linf" = Linf,
  "K"=K,
  "alpha_f_p" = f_p_pars[1], #intercept proba repro
  "beta_f_p_1" = f_p_pars[2], #slope proba repro
  "alpha_f_n" = f_n_pars[1], #intercept fecundity
  "beta_f_n_1" = f_n_pars[2],  #slope fecundity
  "mu_f_d" = f_d_pars[1],        #mean size
  "sigma_f_d" = f_d_pars[2],    #size sd
  f_g = f_g_pars-f_g_pars*0.01,              #survival until from egg to day 30
  para_sensi = "survival_t0_-1"
)

pars_in<-rbind(pars_in_base,s1p,s1m,g1p,g1m,gi1p,gi1m,fp1p,fp1m,fn1p,fn1m,mu1p,mu1m,fg1p,fg1m)

use_proto <- carpobrotus_ipm_unif$proto_ipm

Result_INTER_30<-data.frame (condition= numeric(0),lambda=numeric(0),gen_time=numeric(0),sensi=numeric(0))


for(i in 1:nrow(pars_in)) {
  
  all_pars <- list(
    "alpha_s" = pars_in$alpha_s[i], #intercept survival
    "beta_s" = pars_in$beta_s[i], #slope survival
    "sigma_G" = pars_in$sigma_G[i], #sd resid growth
    "Linf" = pars_in$Linf[i], #slope survival
    "K" = pars_in$K[i], #sd resid growth
    "alpha_f_p" = pars_in$alpha_f_p[i], #intercept proba repro
    "beta_f_p_1" = pars_in$beta_f_p_1[i], #slope proba repro
    "alpha_f_n" = pars_in$alpha_f_n[i], #intercept fecundity
    "beta_f_n_1" = pars_in$beta_f_n_1[i], #slope fecundity
    "mu_f_d" = pars_in$mu_f_d[i],        #mean size
    "sigma_f_d" = pars_in$sigma_f_d[i],    #size sd
    "f_g" = pars_in$f_g[i]              #survival until from egg to day 30
  )
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  
  parameters(use_proto) <- all_pars
  
  boot_ipm <- use_proto %>%
    make_ipm(iterate = TRUE,
             iterations = 2000)
  
  #all_lambdas[i] <- lambda(boot_ipm)
  Tempfile <- data.frame(condition= "INTER_30",lambda = lambda(boot_ipm),
                         gen_time = gen_time(boot_ipm$sub_kernels$P,boot_ipm$sub_kernels$F),
                         sensi=pars_in$para_sensi[i])
  Result_INTER_30 <- rbind(Result_INTER_30,Tempfile)
  
  
}

Result_INTER_30$lambda_ratio<-Result_INTER_30$lambda/Result_INTER_30$lambda[1]
Result_INTER_30$gen_time_ratio<-Result_INTER_30$gen_time/Result_INTER_30$gen_time[1]


write.table(Result_INTER_30, "inter_30_sensitivity.csv", row.names=FALSE,sep=";")



################################## end of sensibility analyses ########################################

################################## Plot the results ###################################################

rm(list=ls())

# Data
sensi<-read.table("IPM_model_sensitivity.csv", header=T,sep=";")

sensi$food<-str_split(sensi$condition,"_",simplify=T)[,1]
sensi$temp<-str_split(sensi$condition,"_",simplify=T)[,2]

# Form ggplot
sensi$temp<-str_replace(sensi$temp,"20","20 \u00B0C")
sensi$temp<-as.factor(str_replace(sensi$temp,"30","30 \u00B0C"))
sensi$food<-str_replace(sensi$food,"conti","Continuous")
sensi$food<-str_replace(sensi$food,"inter","Intermittent")

sensi<-sensi[-which(sensi$sensi=="baseline"),]

# Plot FigS5 of article
ggplot(sensi, aes(x=sensi, y=log(lambda_ratio), color=temp,shape=food)) + 
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(values = c("black","red"))+
  scale_shape_manual(values = c(19,1))+
  labs(y=expression(log~(italic(lambda) ~ ratio)),x="",color="Temperature",shape="Food")+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10,color="black"),
        #axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,),
        legend.title = element_text(size=13),
        legend.position="none",
        legend.text = element_text(size=13),
        axis.text.x = element_blank())

ggplot(sensi, aes(x=sensi, y=log(gen_time_ratio), color=temp,shape=food)) + 
  geom_point(size=2)+
  theme_bw()+
  scale_color_manual(values = c("black","red"))+
  scale_shape_manual(values = c(19,1))+
  labs(y=expression(log~(italic(T) ~ ratio)),x="",color="Temperature",shape="Food")+
  theme(axis.title = element_text(size=13),
        axis.text = element_text(size=10,color="black"),
        axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,),
        #axis.text.x = element_blank(),
        legend.title = element_text(size=13),
        legend.position="none",
        legend.text = element_text(size=13))


########################################################### end of the code ###################################################
