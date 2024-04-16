########## IPM model ########## 
###############################

rm(list=ls()) 

#--------------------------------------------------

# This script is used to perform IPM model for each population

#--------------------------------------------------

# Packages
library(ipmr)
library(Rage)
library(dplyr)
library(FSA)
library(doSNOW)
library(foreach)
library(itertools)
library(tidyverse)
library(plyr)
library(ggplot2)
library(gridExtra)

# Data
# Fecundity
fec<-read.csv(file = "df_fecundity.csv",header=T,sep=";",dec=".")
# Reproduction
reproduction<-read.csv(file = "df_reproduction.csv",header=T,sep=";",dec=".")
# Survival
survival<-read.csv(file = "df_survival_probabilities.csv",header=T,sep=";",dec=".")
# Size newborn
size_newborn<-read.csv(file = "df_size_newborn.csv",header=T,sep=";",dec=".")
# Von Bertallanfy posteriors
load("cfs_bayes.Rdata")

# Functions
vb1 <- vbFuns()
invlogit <- binomial()$linkinv

#--------------------------------------------------



#################################### continuous 20 #########################################
############################################################################################
rm(list = setdiff(ls(),c("survival","reproduction","fec","cfs_bayes","vb1", "invlogit","size_newborn")))

# Select one condition
# Fecundity
fec_i <- fec[which(fec$cond=="conti_20"),]
names(fec_i)[names(fec_i) == 'length_mm'] <- 'lt'
# Size and reproduction 
repro_i <- reproduction[which(reproduction$cond=="conti_20"),] 
names(repro_i)[names(repro_i) == 'size'] <- 'lt'

# Survival
surv_i <- survival[which(survival$cond=="conti_20"),] 
names(surv_i)[names(surv_i) == 'size'] <- 'lt'

# Survival analyses
s_model_1 <- glm(surv ~ log(lt), data = surv_i, 
                 family = binomial("logit"), 
                 weights = total)
s_pars <- coef(s_model_1)


# Growth
# Linf & K from the von bertalanffy model
# posterior distribution
# 1. inter_20
# 2. conti_20
# 3. inter_30
# 4. conti_30

cfs <- cfs_bayes[, grep("\\[2\\]", colnames(cfs_bayes))]
Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)
lt <- vb1(time, Linf, K, t0)


# Von Bertallanfy
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
vbst <- vbst[rowMeans(vbst)>7.435182, ]
# Lt+1 from Lt
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd Lt+1
sigma_g <- apply(log(vbst1), 1, sd)
sigma_g <- mean(sigma_g)


# Reproduction
f_p_model <- glm(repro ~ log(lt), data = repro_i, 
                 family = binomial("logit"))
f_p_pars <- coef(f_p_model)

# Fecundity
f_n_model <- glm(fecundity ~ log(lt) , data = fec_i, 
                 family = poisson())
f_n_pars <- coef(f_n_model)

# size independent vital rates
size_m <- log(mean((size_newborn$size_newborn_mm[which(size_newborn$generation == "T20_F6")]))) 
size_sd<-sd(log((size_newborn$size_newborn_mm[which(size_newborn$generation == "T20_F6")])))
f_d_pars  <- c(size_m,size_sd)

# Hatching rate * survival t0
f_g_pars <- 0.507*0.489 

all_pars <- list(
  "alpha_s"   = s_pars[1],
  "beta_s"    = s_pars[2],
  "Linf"      = Linf,
  "K"         = K,
  "sigma_g"   = sigma_g,
  "alpha_fp"  = f_p_pars[1],
  "beta_fp"   = f_p_pars[2],
  "alpha_fn"  = f_n_pars[1],
  "beta_fn"   = f_n_pars[2],
  "mu_fd"     = f_d_pars[1],
  "sigma_fd"  = f_d_pars[2],
  "f_g"       = f_g_pars,
  "invlogit"  = invlogit
)


carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

carpobrotus_ipm <- define_kernel(  #subkernel for existing fish: they survive (s) and grow (G)
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G, # Subkernel formula
  logit_s   = alpha_s + beta_s*z_1, # z_1 = log(lt)
  s         = invlogit(logit_s),
  g_mu      = Linf*(1-exp(-K)) + exp(-K)*exp(z_1), # link between Lt & Lt+1
  gl_mu     = log(g_mu),
  G         = dnorm(z_2, gl_mu, sigma_g), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)


carpobrotus_ipm <- define_kernel( #subkernel for new recruits
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,        
  f_d       = dnorm(z_2, mu_fd, sigma_fd),   
  logit_fp  = alpha_fp + beta_fp * z_1,
  f_p       = invlogit(logit_fp),
  log_n     = alpha_fn + beta_fn * z_1,
  f_n       = exp(log_n),
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
L <-log(repro_i$lt[repro_i$age == 30])
U <-log(max(fec_i$lt,na.rm=TRUE)) 
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
  iterations      = 11000
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 1000:11000)

# Calculate generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)

use_proto <- carpobrotus_ipm_unif$proto_ipm

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl, c("surv_i", "repro_i", "use_proto", "fec_i", 
                    "invlogit", "m", "cfs",
                    "f_g_pars", "sigma_g", "size_newborn"),
              envir = globalenv())

res_conti_20 <- foreach(u = isplitIndices(1000, chunkSize = 6),
                     .combine = rbind,
                     .inorder = FALSE,
                     .packages = c("ipmr", "Rage", "dplyr")) %dopar% {
                       
                       n <- length(u)                 
                       l <- g <- numeric(n)
                       
                       
                       Ng <- nrow(cfs)
                       cL <- grep("Linf", colnames(cfs))
                       cK <- grep("K", colnames(cfs))
                       ct0 <- grep("t0", colnames(cfs))
                       
                       for(i in 1:n) {
                         size_newboot <- sample_n(size_newborn[size_newborn$generation == "T20_F6", ],
                                                  sum(size_newborn$generation == "T20_F6"))
                         size_m <- log(mean(size_newboot$size_newborn_mm))
                         size_sd<-sd(log(size_newboot$size_newborn_mm))  
                         f_d_pars  <- c(size_m,size_sd)
                         
                         # Survival data
                         surv_i_boot<-sample_n(surv_i, size = nrow(surv_i), replace = T)
                         s_model_1  <- glm(surv ~ log(lt),
                                           data = surv_i_boot, 
                                           family = binomial("logit"), 
                                           weights = total)
                         s_pars <- coef(s_model_1)
                         
                         repro_i_boot <- sample_n(repro_i, size = nrow(repro_i), replace = T)
                         
                         # Growth
                         N  <- sample(Ng, 1)
                         Linf <- cfs[N, cL]
                         K <- cfs[N, cK]
                         
                         # Reproduction
                         f_p_model <- glm(repro ~ log(lt) ,
                                          data = repro_i_boot,
                                          family = binomial("logit"))
                         f_p_pars <- coef(f_p_model)
                         
                         # Fecundity
                         fec_i_boot <- sample_n(fec_i, size = nrow(fec_i), replace = T)
                         f_n_model <- glm(fecundity ~ log(lt) ,
                                          data = fec_i_boot,
                                          family = poisson())
                         f_n_pars <- coef(f_n_model)
                         
                         all_pars <- list(
                           "alpha_s"   = s_pars[1],
                           "beta_s"    = s_pars[2],
                           "Linf"      = Linf,
                           "K"         = K,
                           "sigma_g"   = sigma_g,
                           "alpha_fp"  = f_p_pars[1],
                           "beta_fp"   = f_p_pars[2],
                           "alpha_fn"  = f_n_pars[1],
                           "beta_fn"   = f_n_pars[2],
                           "mu_fd"     = f_d_pars[1],
                           "sigma_fd"  = f_d_pars[2],
                           "f_g"       = f_g_pars,
                           "invlogit"  = invlogit
                         )
                         
                         # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
                         parameters(use_proto) <- all_pars
                         boot_ipm <- make_ipm(use_proto, iterate = TRUE,
                                              iterations = 6000)
                         
                         #all_lambdas[i] <- lambda(boot_ipm)
                         l[i] <- lambda(boot_ipm)
                         g[i] <- gen_time(boot_ipm$sub_kernels$P, boot_ipm$sub_kernels$F)
                         
                       }
                       Result_CONTI_20 <- data.frame(condition= "CONTI_20",
                                                  lambda = l,
                                                  gen_time = g)
                       Result_CONTI_20
                     }

write.table(res_conti_20,file="res_conti_20.csv",row.names = F,sep=";")

stopCluster(cl)


#################################### intermittent 20 #########################################
#############################################################################################
rm(list = setdiff(ls(),c("survival","reproduction","fec","cfs_bayes","vb1", "invlogit","size_newborn")))


# Select one condition
# Fecundity
fec_i <- fec[which(fec$cond=="inter_20"),]
names(fec_i)[names(fec_i) == 'length_mm'] <- 'lt'
# Size and reproduction 
repro_i <- reproduction[which(reproduction$cond=="inter_20"),] 
names(repro_i)[names(repro_i) == 'size'] <- 'lt'
# Survival
surv_i <- survival[which(survival$cond=="inter_20"),]
names(surv_i)[names(surv_i) == 'size'] <- 'lt'

# Survival analyses
s_model_1 <- glm(surv ~ log(lt), data = surv_i, 
                 family = binomial("logit"), 
                 weights = total)
s_pars <- coef(s_model_1)


# Growth
# Linf & K from the von bertalanffy model
# posterior distribution
# 1. inter_20
# 2. conti_20
# 3. inter_30
# 4. conti_30

cfs <- cfs_bayes[, grep("\\[1\\]", colnames(cfs_bayes))]
Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)
lt <- vb1(time, Linf, K, t0)

# Von Bertallanfy
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
vbst <- vbst[rowMeans(vbst)>5.82335, ]
# Lt+1 from Lt
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd Lt+1
sigma_g <- apply(log(vbst1), 1, sd)
sigma_g <- mean(sigma_g)


# Reproduction
f_p_model <- glm(repro ~ log(lt), data = repro_i, 
                 family = binomial("logit"))
f_p_pars <- coef(f_p_model)

# Fecundity
f_n_model <- glm(fecundity ~ log(lt) , data = fec_i, 
                 family = poisson())
f_n_pars <- coef(f_n_model)

# Size independent vital rates

size_m <- log(mean((size_newborn$size_newborn_mm[which(size_newborn$generation == "T20_F6")]))) 
size_sd<-sd(log((size_newborn$size_newborn_mm[which(size_newborn$generation == "T20_F6")]))) 
f_d_pars  <- c(size_m,size_sd)

# Hatching rate * survival t0
f_g_pars <- 0.486*0.489

all_pars <- list(
  "alpha_s"   = s_pars[1],
  "beta_s"    = s_pars[2],
  "Linf"      = Linf,
  "K"         = K,
  "sigma_g"   = sigma_g,
  "alpha_fp"  = f_p_pars[1],
  "beta_fp"   = f_p_pars[2],
  "alpha_fn"  = f_n_pars[1],
  "beta_fn"   = f_n_pars[2],
  "mu_fd"     = f_d_pars[1],
  "sigma_fd"  = f_d_pars[2],
  "f_g"       = f_g_pars,
  "invlogit"  = invlogit
)


carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

carpobrotus_ipm <- define_kernel(  #subkernel for existing fish: they survive (s) and grow (G)
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G, # Subkernel formula
  logit_s   = alpha_s + beta_s*z_1, # z_1 = log(lt)
  s         = invlogit(logit_s),
  g_mu      = Linf*(1-exp(-K)) + exp(-K)*exp(z_1), # link between Lt & Lt+1
  gl_mu     = log(g_mu),
  G         = dnorm(z_2, gl_mu, sigma_g), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( #subkernel for new recruits
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,         
  f_d       = dnorm(z_2, mu_fd, sigma_fd),   
  logit_fp  = alpha_fp + beta_fp * z_1,
  f_p       = invlogit(logit_fp), # probability of reproduction
  log_n     = alpha_fn + beta_fn * z_1,
  f_n       = exp(log_n),
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
L <-log(repro_i$lt[repro_i$age == 30])
U <-  log(max(fec_i$lt,na.rm = T)) 
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
  iterations      = 11000
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 2000:4000)

# Calculate generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)

use_proto <- carpobrotus_ipm_unif$proto_ipm

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl, c("surv_i", "repro_i", "use_proto", "fec_i", 
                    "invlogit", "m", "cfs",
                    "f_g_pars", "sigma_g", "size_newborn"),
              envir = globalenv())

res_inter_20 <- foreach(u = isplitIndices(1000, chunkSize = 6),
                     .combine = rbind,
                     .inorder = FALSE,
                     .packages = c("ipmr", "Rage", "dplyr"))%dopar% {
  
  n <- length(u)                 
  l <- g <- numeric(n)
  
  
  Ng <- nrow(cfs)
  cL <- grep("Linf", colnames(cfs))
  cK <- grep("K", colnames(cfs))
  ct0 <- grep("t0", colnames(cfs))
  
  for(i in 1:n) {
    size_newboot <- sample_n(size_newborn[size_newborn$generation == "T20_F6", ],
                             sum(size_newborn$generation == "T20_F6"))
    size_m <- log(mean(size_newboot$size_newborn_mm)) 
    size_sd<-sd(log(size_newboot$size_newborn_mm))
    f_d_pars  <- c(size_m,size_sd)
    
    # Survival data
    surv_i_boot<-sample_n(surv_i, size = nrow(surv_i), replace = T)
    s_model_1  <- glm(surv ~ log(lt),
                      data = surv_i_boot, 
                      family = binomial("logit"), 
                      weights = total)
    s_pars <- coef(s_model_1)
    
    repro_i_boot <- sample_n(repro_i, size = nrow(repro_i), replace = T)
    # Growth
    N  <- sample(Ng, 1)
    Linf <- cfs[N, cL]
    K <- cfs[N, cK]
    
    # Reproduction
    f_p_model <- glm(repro ~ log(lt) ,
                     data = repro_i_boot,
                     family = binomial("logit"))
    f_p_pars <- coef(f_p_model)
    
    # Fecundity
    fec_i_boot <- sample_n(fec_i, size = nrow(fec_i), replace = T)
    f_n_model <- glm(fecundity ~ log(lt) ,
                     data = fec_i_boot,
                     family = poisson())
    f_n_pars <- coef(f_n_model)
    
    all_pars <- list(
      "alpha_s"   = s_pars[1],
      "beta_s"    = s_pars[2],
      "Linf"      = Linf,
      "K"         = K,
      "sigma_g"   = sigma_g,
      "alpha_fp"  = f_p_pars[1],
      "beta_fp"   = f_p_pars[2],
      "alpha_fn"  = f_n_pars[1],
      "beta_fn"   = f_n_pars[2],
      "mu_fd"     = f_d_pars[1],
      "sigma_fd"  = f_d_pars[2],
      "f_g"       = f_g_pars,
      "invlogit"  = invlogit
    )
    
    # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
    parameters(use_proto) <- all_pars
    boot_ipm <- make_ipm(use_proto, iterate = TRUE,
                         iterations = 4000)
    
    #all_lambdas[i] <- lambda(boot_ipm)
    l[i] <- lambda(boot_ipm)
    g[i] <- gen_time(boot_ipm$sub_kernels$P, boot_ipm$sub_kernels$F)
  }
  Result_INTER_20 <- data.frame(condition= "INTER_20",
                              lambda = l,
                              gen_time = g)
  Result_INTER_20
                     }

write.table(res_inter_20,file="res_inter_20.csv",row.names = F,sep=";")
stopCluster(cl)

#################################### continuous 30 ##########################################
#############################################################################################
rm(list = setdiff(ls(),c("survival","reproduction","fec","cfs_bayes","vb1", "invlogit","size_newborn")))


# Select one condition
# Fecundity
fec_i <- fec[which(fec$cond=="conti_30"),]
names(fec_i)[names(fec_i) == 'length_mm'] <- 'lt'
# Size and reproduction 
repro_i <- reproduction[which(reproduction$cond=="conti_30"),] 
names(repro_i)[names(repro_i) == 'size'] <- 'lt'
# Survival
surv_i <- survival[which(survival$cond=="conti_30"),]
names(surv_i)[names(surv_i) == 'size'] <- 'lt'

# Survival analyses
s_model_1 <- glm(surv ~ log(lt), data = surv_i, 
                 family = binomial("logit"), 
                 weights = total)
s_pars <- coef(s_model_1)


# Growth
# Linf & K from the von bertalanffy model
# posterior distribution
# 1. inter_20
# 2. conti_20
# 3. inter_30
# 4. conti_30

cfs <- cfs_bayes[, grep("\\[4\\]", colnames(cfs_bayes))]
Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)
lt <- vb1(time, Linf, K, t0)

# Von Bertallanfy
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
vbst <- vbst[rowMeans(vbst)>8.3597, ]
# Lt+1 from Lt
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd Lt+1
sigma_g <- apply(log(vbst1), 1, sd)
sigma_g <- mean(sigma_g)

# Repro
f_p_model <- glm(repro ~ log(lt), data = repro_i, 
                 family = binomial("logit"))
f_p_pars <- coef(f_p_model)

# Fecundity
f_n_model <- glm(fecundity ~ log(lt) , data = fec_i, 
                 family = poisson()) 
f_n_pars <- coef(f_n_model)

# Size independent vital rates
size_m <- log(mean((size_newborn$size_newborn_mm[which(size_newborn$generation == "T30_F10")]))) 
size_sd<-sd(log((size_newborn$size_newborn_mm[which(size_newborn$generation == "T30_F10")])))
f_d_pars  <- c(size_m,size_sd)

# Hatching rate * survival t0
f_g_pars <- 0.468*0.553

all_pars <- list(
  "alpha_s"   = s_pars[1],
  "beta_s"    = s_pars[2],
  "Linf"      = Linf,
  "K"         = K,
  "sigma_g"   = sigma_g,
  "alpha_fp"  = f_p_pars[1],
  "beta_fp"   = f_p_pars[2],
  "alpha_fn"  = f_n_pars[1],
  "beta_fn"   = f_n_pars[2],
  "mu_fd"     = f_d_pars[1],
  "sigma_fd"  = f_d_pars[2],
  "f_g"       = f_g_pars,
  "invlogit"  = invlogit
)


carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

carpobrotus_ipm <- define_kernel(  #subkernel for existing fish: they survive (s) and grow (G)
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G, # Subkernel formula
  logit_s   = alpha_s + beta_s*z_1, # z_1 = log(lt)
  s         = invlogit(logit_s),
  g_mu      = Linf*(1-exp(-K)) + exp(-K)*exp(z_1), # link between Lt & Lt+1
  gl_mu     = log(g_mu),
  G         = dnorm(z_2, gl_mu, sigma_g), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( #subkernel for new recruits
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,         
  f_d       = dnorm(z_2, mu_fd, sigma_fd),   
  logit_fp  = alpha_fp + beta_fp * z_1,
  f_p       = invlogit(logit_fp), 
  log_n     = alpha_fn + beta_fn * z_1,
  f_n       = exp(log_n),
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


L <-log(repro_i$lt[repro_i$age == 30])
U <-  log(max(fec_i$lt,na.rm = T)) 
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

# conv_plot(carpobrotus_ipm_unif, iterations = 2000:4000)

# Calculate generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)

use_proto <- carpobrotus_ipm_unif$proto_ipm

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl, c("surv_i", "repro_i", "use_proto", "fec_i", 
                    "invlogit", "m", "cfs",
                    "f_g_pars", "sigma_g", "size_newborn"),
              envir = globalenv())

res_conti_30 <- foreach(u = isplitIndices(1000, chunkSize = 6),
                     .combine = rbind,
                     .inorder = FALSE,
                     .packages = c("ipmr", "Rage", "dplyr")) %dopar% {
  
  n <- length(u)                 
  l <- g <- numeric(n)
  
  
  Ng <- nrow(cfs)
  cL <- grep("Linf", colnames(cfs))
  cK <- grep("K", colnames(cfs))
  ct0 <- grep("t0", colnames(cfs))
  for(i in 1:n) {
    size_newboot <- sample_n(size_newborn[size_newborn$generation == "T30_F10", ],
                             sum(size_newborn$generation == "T30_F10"))
    size_m <- log(mean(size_newboot$size_newborn_mm)) 
    size_sd<-sd(log(size_newboot$size_newborn_mm))  
    f_d_pars  <- c(size_m,size_sd)
    
    # Survival data
    surv_i_boot<-sample_n(surv_i, size = nrow(surv_i), replace = T)
    s_model_1  <- glm(surv ~ log(lt),
                      data = surv_i_boot, 
                      family = binomial("logit"), 
                      weights = total)
    s_pars <- coef(s_model_1)
    
    repro_i_boot <- sample_n(repro_i, size = nrow(repro_i), replace = T)
    
    # Growth
    N  <- sample(Ng, 1)
    Linf <- cfs[N, cL]
    K <- cfs[N, cK]
    
    # Reproduction
    f_p_model <- glm(repro ~ log(lt) ,
                     data = repro_i_boot,
                     family = binomial("logit"))
    f_p_pars <- coef(f_p_model)

    # Fecundity
    fec_i_boot <- sample_n(fec_i, size = nrow(fec_i), replace = T)
    f_n_model <- glm(fecundity ~ log(lt) ,
                     data = fec_i_boot,
                     family = poisson())
    f_n_pars <- coef(f_n_model)
    
    all_pars <- list(
      "alpha_s"   = s_pars[1],
      "beta_s"    = s_pars[2],
      "Linf"      = Linf,
      "K"         = K,
      "sigma_g"   = sigma_g,
      "alpha_fp"  = f_p_pars[1],
      "beta_fp"   = f_p_pars[2],
      "alpha_fn"  = f_n_pars[1],
      "beta_fn"   = f_n_pars[2],
      "mu_fd"     = f_d_pars[1],
      "sigma_fd"  = f_d_pars[2],
      "f_g"       = f_g_pars,
      "invlogit"  = invlogit
    )
    
    # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
    parameters(use_proto) <- all_pars
    boot_ipm <- make_ipm(use_proto, iterate = TRUE,
                         iterations = 1500) 
    l[i] <- lambda(boot_ipm)
    g[i] <- gen_time(boot_ipm$sub_kernels$P, boot_ipm$sub_kernels$F)
    }
  Result_CONTI_30 <- data.frame(condition= "CONTI_30",
                              lambda = l,
                              gen_time = g)
  Result_CONTI_30

}

write.table(res_conti_30,file="res_conti_30.csv",row.names = F,sep=";")

stopCluster(cl)

#################################### intermittent 30 ##########################################
##############################################################################################
rm(list = setdiff(ls(),c("survival","reproduction","fec","cfs_bayes","vb1", "invlogit","size_newborn")))


# Select one condition
# Fecundity
fec_i <- fec[which(fec$cond=="inter_30"),]
names(fec_i)[names(fec_i) == 'length_mm'] <- 'lt'
# Size and reproduction 
repro_i <- reproduction[which(reproduction$cond=="inter_30"),] 
names(repro_i)[names(repro_i) == 'size'] <- 'lt'
# Survival
surv_i <- survival[which(survival$cond=="inter_30"),]
names(surv_i)[names(surv_i) == 'size'] <- 'lt'

##survival analyses
s_model_1 <- glm(surv ~ log(lt), data = surv_i, 
                 family = binomial("logit"), 
                 weights = total)
s_pars <- coef(s_model_1)


# Growth
# Linf & K from the von bertalanffy model
# posterior distribution
# 1. inter_20
# 2. conti_20
# 3. inter_30
# 4. conti_30

cfs <- cfs_bayes[, grep("\\[3\\]", colnames(cfs_bayes))]
Linf <- mean(cfs[, grep("Linf", colnames(cfs))])
K <- mean(cfs[, grep("K", colnames(cfs))])
t0 <- mean(cfs[, grep("t0", colnames(cfs))])

time <- seq(0, 350, le = 601)
lt <- vb1(time, Linf, K, t0)


# Von Bertallanfy
vbst <- matrix(numeric(), length(time), nrow(cfs))
for (i in 1:nrow(cfs)) {
  vbst[, i] <- cfs[i, 2]*(1-exp(-cfs[i, 1]*(time-cfs[i, 3])))
}
vbst <- vbst[rowMeans(vbst)>7.5970, ] 
# Lt+1 from Lt
vbst1 <- vbst*exp(-K)+(Linf*(1-exp(-K)))
# sd Lt+1 
sigma_g <- apply(log(vbst1), 1, sd)
sigma_g <- mean(sigma_g)

# Repro
f_p_model <- glm(repro ~ log(lt), data = repro_i, 
                 family = binomial("logit"))
f_p_pars <- coef(f_p_model)

# Fecundity
f_n_model <- glm(fecundity ~ log(lt) , data = fec_i, 
                 family = poisson()) 
f_n_pars <- coef(f_n_model)

# size independent vital rates
size_m <- log(mean((size_newborn$size_newborn_mm[which(size_newborn$generation == "T30_F10")])))
size_sd<-sd(log((size_newborn$size_newborn_mm[which(size_newborn$generation == "T30_F10")])))
f_d_pars  <- c(size_m,size_sd)

# Hatching rate * survival t0
f_g_pars <- 0.364*0.553

all_pars <- list(
  "alpha_s"   = s_pars[1],
  "beta_s"    = s_pars[2],
  "Linf"      = Linf,
  "K"         = K,
  "sigma_g"   = sigma_g,
  "alpha_fp"  = f_p_pars[1],
  "beta_fp"   = f_p_pars[2],
  "alpha_fn"  = f_n_pars[1],
  "beta_fn"   = f_n_pars[2],
  "mu_fd"     = f_d_pars[1],
  "sigma_fd"  = f_d_pars[2],
  "f_g"       = f_g_pars,
  "invlogit"  = invlogit
)


carpobrotus_ipm <- init_ipm(sim_gen = "simple", di_dd = "di", det_stoch = "det")

carpobrotus_ipm <- define_kernel(  #subkernel for existing fish: they survive (s) and grow (G)
  proto_ipm = carpobrotus_ipm,
  name      = "P",
  formula   = s * G, # Subkernel formula
  logit_s   = alpha_s + beta_s*z_1, # z_1 = log(lt)
  s         = invlogit(logit_s),
  g_mu      = Linf*(1-exp(-K)) + exp(-K)*exp(z_1), # link between Lt & Lt+1
  gl_mu     = log(g_mu),
  G         = dnorm(z_2, gl_mu, sigma_g), 
  data_list = all_pars,
  uses_par_sets = FALSE,
  states    = list(c("z")),
  evict_cor = TRUE,
  evict_fun = truncated_distributions(fun = "norm",
                                      target = "G")
)

carpobrotus_ipm <- define_kernel( ##subkernel for new recruits
  proto_ipm = carpobrotus_ipm, 
  name      = "F",
  formula   = f_p*f_n*f_d*f_g,         
  f_d       = dnorm(z_2, mu_fd, sigma_fd),   
  logit_fp  = alpha_fp + beta_fp * z_1,
  f_p       = invlogit(logit_fp), # probability of reproduction
  log_n     = alpha_fn + beta_fn * z_1,
  f_n       = exp(log_n),
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
L <-log(repro_i$lt[repro_i$age == 30])
U <-  log(max(fec_i$lt,na.rm=T)) 
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
  iterations      = 11000
)

is_conv_to_asymptotic(carpobrotus_ipm_unif, tol = 1e-9)

conv_plot(carpobrotus_ipm_unif, iterations = 1:2000)

# Calculate generation time from package rage
gen_time(carpobrotus_ipm_unif$sub_kernels$P,carpobrotus_ipm_unif$sub_kernels$F)

lambda_obs <- lambda(carpobrotus_ipm_unif)

use_proto <- carpobrotus_ipm_unif$proto_ipm

cl <- makeCluster(6, type = "SOCK")
registerDoSNOW(cl)
clusterExport(cl, c("surv_i", "repro_i", "use_proto", "fec_i", 
                    "invlogit", "m", "cfs",
                    "f_g_pars", "sigma_g", "size_newborn"),
              envir = globalenv())

res_inter_30 <- foreach(u = isplitIndices(1000, chunkSize = 6),
                     .combine = rbind,
                     .inorder = FALSE,
                     .packages = c("ipmr", "Rage", "dplyr")) %dopar% {
  
  n <- length(u)                 
  l <- g <- numeric(n)
  
  
  Ng <- nrow(cfs)
  cL <- grep("Linf", colnames(cfs))
  cK <- grep("K", colnames(cfs))
  ct0 <- grep("t0", colnames(cfs))
  
  for(i in 1:n) {
    size_newboot <- sample_n(size_newborn[size_newborn$generation == "T20_F6", ],
                             sum(size_newborn$generation == "T30_F10"))
    size_m <- log(mean(size_newboot$size_newborn_mm)) 
    size_sd<-sd(log(size_newboot$size_newborn_mm))  
    f_d_pars  <- c(size_m,size_sd)
    
    # Survival data
    surv_i_boot<-sample_n(surv_i, size = nrow(surv_i), replace = T)
    s_model_1  <- glm(surv ~ log(lt),
                      data = surv_i_boot, 
                      family = binomial("logit"), 
                      weights = total)
    s_pars <- coef(s_model_1)
    
    repro_i_boot <- sample_n(repro_i, size = nrow(repro_i), replace = T)
    # Growth
    N  <- sample(Ng, 1)
    Linf <- cfs[N, cL]
    K <- cfs[N, cK]
    
    # Reproduction
    f_p_model <- glm(repro ~ log(lt) ,
                     data = repro_i_boot,
                     family = binomial("logit"))
    f_p_pars <- coef(f_p_model)
    
    # Fecundity
    fec_i_boot <- sample_n(fec_i, size = nrow(fec_i), replace = T)
    f_n_model <- glm(fecundity ~ log(lt) ,
                     data = fec_i_boot,
                     family = poisson())
    f_n_pars <- coef(f_n_model)
    
    all_pars <- list(
      "alpha_s"   = s_pars[1],
      "beta_s"    = s_pars[2],
      "Linf"      = Linf,
      "K"         = K,
      "sigma_g"   = sigma_g,
      "alpha_fp"  = f_p_pars[1],
      "beta_fp"   = f_p_pars[2],
      "alpha_fn"  = f_n_pars[1],
      "beta_fn"   = f_n_pars[2],
      "mu_fd"     = f_d_pars[1],
      "sigma_fd"  = f_d_pars[2],
      "f_g"       = f_g_pars,
      "invlogit"  = invlogit
    )
    
    # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
    parameters(use_proto) <- all_pars
    boot_ipm <- make_ipm(use_proto, iterate = TRUE,
                         iterations = 1500)
    
    l[i] <- lambda(boot_ipm)
    g[i] <- gen_time(boot_ipm$sub_kernels$P, boot_ipm$sub_kernels$F)
  }
  Result_INTER_30 <- data.frame(condition= "RES_30",
                             lambda = l,
                             gen_time = g)
  Result_INTER_30
                     }
write.table(res_inter_30,file="res_inter_30.csv",row.names = F,sep=";")
stopCluster(cl)


#--------------------------------------------------
# Plot the demographic parameters (Fig 5 of the article)


rm(list = ls())

# IPM model outputs
df_ipm<-read.csv(file = "IPM_model_outputs.csv",header=T,sep=";",dec=".")

# We added 30 days to the generation time compared to the IPMs outputs 
# since the models started from this age
df_ipm$gen_time<-df_ipm$gen_time+30

# Form to plot
df_ipm<-gather(df_ipm,key="param",value="value",lambda:gen_time)

# df ggplot
df <- ddply(df_ipm, .(condition,temp,food,param), summarise,
            mean       = mean(value),
            median     = median(value),
            q2.5       = quantile(value,0.025),
            q97.5       = quantile(value,0.975))
df
df$temp<-as.character(df$temp)

# Form ggplot
df$temp<-str_replace(df$temp,"20","20 \u00B0C")
df$temp<-as.factor(str_replace(df$temp,"30","30 \u00B0C"))
df$food<-str_replace(df$food,"conti","Continuous")
df$food<-str_replace(df$food,"inter","Intermittent")

# Plot generation time T
p1<-ggplot(df[df$param=="gen_time",],aes(x=condition,y=median,color=temp,shape=food))+
  geom_pointrange(aes(ymin=q2.5,ymax=q97.5),size=0.5)+
  labs(y=bquote(italic(T)),color="Temperature",shape="Food")+
  theme_bw()+
  scale_color_manual(values=c("black","red"))+
  scale_shape_manual(values = c(19,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        axis.title.y = element_text(face="italic"))

# Plot asymptotic per capita population growth rate λ
p2<-ggplot(df[df$param=="lambda",],aes(x=condition,y=median,color=temp,shape=food))+
  geom_pointrange(aes(ymin=q2.5,ymax=q97.5),size=0.5)+
  labs(y=expression(lambda))+
  theme_bw()+
  scale_color_manual(values=c("black","red"))+
  scale_shape_manual(values = c(19,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        axis.title.y = element_text(face="italic"))


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
lgd<-get_legend(p1)

grid.arrange(p1,p2,lgd,ncol=3,nrow=1,widths=c(1.3, 1.3, 0.5))

#--------------------------------------------------
# End of the code

