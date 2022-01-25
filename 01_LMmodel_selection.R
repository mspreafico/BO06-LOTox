##########################################################
# 4.1 Latent Markov model for longitudinal toxicity data #
##########################################################
rm(list=ls())
library(data.table)
library(emdbook)
library(LMest)
library(xtable)
setwd("~/OneDrive - Politecnico di Milano/PhD - Leiden/BO06-LOTox")
load('data/BO06_tox_4_2cat.Rdata')

# Generic toxicities: naus, inf, oral (j = 1, 2, 3)
#   categories: Cj = {0 : none, 1 : mild, 2 : moderate, 3 : severe}

# Drug-specific toxicities: car, oto, neur(j = 4, 5, 6): 
#   categories: Cj = {0 : no, 1 : yes}

#######################
# (I) MODEL SELECTION #
#######################

# M1: Unrestricted LM model without covariates #
#----------------------------------------------#
out.tdep <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                        latentFormula = NULL,
                        version = 'categorical',
                        index = c("patid","cycno"),
                        data = df,
                        k = 1:10,
                        nrep = 20,
                        tol1 = 10^-8,
                        seed = 1234)
out.tdep
g <- NULL
for(k in 1:10){
  g <- c(g, out.tdep$out.single[[k]]$np)
}
data.out.tdep <- cbind.data.frame(k = out.tdep$k,
                                  g = g,
                                  lk = out.tdep$lkv,
                                  AIC = out.tdep$Aic,
                                  BIC = out.tdep$Bic)
print(xtable(data.out.tdep, digits=c(0,0,0,2,2,2)),  include.rownames=FALSE)
# k = 4; BIC = 16728.90 <- BEST k
# k=4, g=111, l=-8035.21, AIC=16292.42, BIC=16728.90 


# M2: Multinomial logit LM model without covariates #
#---------------------------------------------------#
out0 <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                    latentFormula = ~ 1 | 1,
                    version = 'categorical',
                    index = c("patid","cycno"),
                    data = df,
                    k = 4,
                    nrep = 20,
                    tol1 = 10^-8,
                    seed = 1234)
out0
# BIC = 16512.16
# k=4, g=63, l=-8069.21, AIC=16264.43, BIC=16512.16


# M3: M2 + regimen effect on initial prob. #
#------------------------------------------#
out.d.trt <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ I(0+(trt=='Reg-DI')) | NULL,
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.d.trt
# BIC = 16522.5
# k=4, g=66, l=-8065.49, AIC=16262.97, BIC=16522.50


# M4: M2 + gender effect on initial prob. #
#-----------------------------------------#
out.d.sex <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ I(0+(sex=='male')) | NULL,
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.d.sex
# BIC =  16514.98
# k=4, g=66, l=-8061.73, AIC=16255.45, BIC=16514.98


# M5: M2 + age effect on initial prob. #	
#--------------------------------------#
out.d.age <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ I(age_in-15) | NULL,
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.d.age
# BIC = 16502.20 <---- BEST MODEL
# k=4, g=66, l=-8055.35, AIC=16242.69, BIC=16502.20


# M6: M2 + regimen effect on transition prob. #
#---------------------------------------------#
out.t.trt <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ NULL | I(0+(trt=='Reg-DI')),
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.t.trt
# BIC = 16571.66
# k=4, g=75, l=-8063.37, AIC=16276.74, BIC=16571.66


# M7: M2 + gender effect on transition prob. #
#--------------------------------------------#
out.t.sex <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ NULL | I(0+(sex=='male')),
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.t.sex
# BIC = 16565.58
# k=4, g=75, l=-8060.33, AIC=16270.66, BIC=16565.58


# M8: M2 + age effect on transition prob. #
#-----------------------------------------#
out.t.age <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ NULL | I(age_in-15),
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.t.age
# BIC = 16567.06
# k=4, g=75, l=-8061.07, AIC=16272.14, BIC=16567.06


# M9: M2 + time-varying chemotherapy dose on both prob. #
#-------------------------------------------------------#
out.delta <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ I(delta-100),
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = df,
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         seed = 1234)
out.delta
# BIC = 16553.82
# k=4, g=78, l=-8045.551, AIC=16247.1, BIC=16553.82


# M10: M2 + time-varying WBC count on both prob. #
#------------------------------------------------#
out.wbc <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                       latentFormula = ~  I(log(WBC)-1.8),
                       version = 'categorical',
                       index = c("patid","cycno"),
                       data = df,
                       k = 4,
                       nrep = 20,
                       tol1 = 10^-8,
                       seed = 1234)
out.wbc
# BIC = 16587.78
# k=4, g=78, l=-8062.534, AIC=16281.07, BIC=16587.78


# M11: M2 + time-varying PLT count on both prob. #
#------------------------------------------------#
out.plt <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                       latentFormula = ~  I(log(PLT)-5.5),
                       version = 'categorical',
                       index = c("patid","cycno"),
                       data = df,
                       k = 4,
                       nrep = 20,
                       tol1 = 10^-8,
                       seed = 1234)
out.plt
# BIC = 16557.02
# k=4, g=78, l=-8047.152, AIC=16250.3, BIC=16557.02


# M11: M2 + time-varying NEUT count on both prob. #
#-------------------------------------------------#
out.neut <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                        latentFormula = ~  I(log(NEUT)-1.7),
                        version = 'categorical',
                        index = c("patid","cycno"),
                        data = df,
                        k = 4,
                        nrep = 20,
                        tol1 = 10^-8,
                        seed = 1234)
out.neut
# BIC = 16588.05
# k=4, g=78, l=-8062.67, AIC=16281.34, BIC=16588.05


############################
# FINAL MODEL: LM model M5 #
############################
# Re-compute model M5 with SE
out.final <- lmestSearch(responsesFormula = naus + inf + oral + car + oto + neur  ~ NULL,
                         latentFormula = ~ I(age_in-15) | NULL,
                         version = 'categorical',
                         index = c("patid","cycno"),
                         data = data.frame(df),
                         k = 4,
                         nrep = 20,
                         tol1 = 10^-8,
                         out_se = T, # with SE = TRUE
                         seed = 1)

final <- out.final$out.single[[1]]
#rm(list=setdiff(ls(), c("final")))
#save.image('results/LM_final_M5.Rdata')

# Figure 3: Estimated conditional response probabilities psi_jy|u for the final LM model M5
final$Psi
# Based on these results, the following LOTox states labelling were derived:
#   - State 1: quite good conditions (non-toxic) -> no LOTox
#   - State 2: non-severe nausea with possible drug-specific AEs -> moderate LOTox
#   - State 3: moderate/severe nausea/vomiting only -> low LOTox (limited to nausea)
#   - State 4: multiple severe/moderate generic toxicities -> high LOTox


# Table 3: Estimated regression parameters affecting the distribution of the initial 
#          probabilities in Equation (9)
round(final$Be,4)

# Table 3: Estimated transition probabilities in Equation (10)
tau <- final$PI[,,1,2]
round(tau,4)

# Figure 4 [left panel]: Estimated initial probabilities of latent states for patients 
#                        aged 10, 15 and 20 years old and average vector of the initial 
#                        probabilities over all the 377 subjects in the sample.
IP <- matrix(c(as.vector(final$Piv[which(unique(df$patid)==4),]), #10
               as.vector(final$Piv[which(unique(df$patid)==23),]), #15
               as.vector(final$Piv[which(unique(df$patid)==24),]), #20 
               as.vector(colMeans(final$Piv)) # average vector
), nrow=4, ncol=4, byrow=T)
colnames(IP)<-1:4 # Latent States
rownames(IP)=c('10','15','20','mean')
round(IP*100,1)

# Figure 4 [right panel]: Latent states prevalences over cycles averaged over all the subjects.
mPIV <- matrix(NA, nrow=6, ncol=4)
mPIV[1,] <- colMeans(final$Piv)
tau <- final$PI[,,1,2]
for(t in 2:6){
  mPIV[t,] <- mPIV[t-1,]%*%tau
}
colnames(mPIV)=1:4 # Latent states
rownames(mPIV)=1:6 # Cycles
mPIV
round(mPIV*100,1)

