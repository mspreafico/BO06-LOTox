#########################################################
# 4.2 Longitudinal profiles of Latent Overall Toxicity #
#########################################################
rm(list=ls())
library(compositions)
library(data.table)
setwd('~/github/BO06-LOTox/')
load('results/LM_final_M5.Rdata') # NOT PUBLICY AVAILABLE

source('~/github/BO06-LOTox/utils.R')


#####################################################################
# FAKE RESULTS
#-------------------------------------------------------------------
# Un-comment next line if you want to use the final model 
# obtained by 'fake_dataset.Rdata':
# load('fake_data_results/fake_LM_final_M5.Rdata')
#-------------------------------------------------------------------
# WARNING: Results for fake data are different than for original ones
#####################################################################


# Longitudinal Probability profiles of LOTox (P-LOTox) #
#------------------------------------------------------#
# Vector of subjects
pts <- unique(final$data$patid)

# Create dataset with P_LOTox profiles and LOTox sequences
P_LOTox <- NULL
for(i in pts){
  pt_list <- lmestDecoding.Vlogit(final,sequence=which(pts==i))
  P_LOTox<-rbind.data.frame(P_LOTox,
                            cbind.data.frame('patid' = rep(i, final$TT),
                                             'cycno' = c(1:final$TT),
                                             'Ul' = pt_list$Ul, # LOTox sequence by local decoding
                                             'Ug' = pt_list$Ug, # LOTox sequence by global decoding
                                             matrix(unlist(pt_list$V), nrow = final$TT, byrow = TRUE))
  )
}
P_LOTox <- data.table(P_LOTox)
colnames(P_LOTox) <- c('patid','cycno','Ul','Ug','minor1','moderate2','low3','high4')
P_LOTox


# Longitudinal Relative Risk profiles of LOTox (RR-LOTox) #
#---------------------------------------------------------#
# Probability column vectors p_i^(t) as Aitchison compositions
prob_vec = acomp(P_LOTox[, c("minor1","low3", "moderate2", "high4")])

# Compute the additive log ratio transform of the compositions prob_t
# where ivar is the column related to the reference latent state R
rr_vec = alr(prob_vec, ivar=1)

# Create RR_LOTox profiles dataset
RR_LOTox <- cbind(P_LOTox[,.(patid,cycno)], rr_vec)
RR_LOTox

#rm(list=setdiff(ls(), c("P_LOTox","RR_LOTox")))
#save.image('results/LOTox_profiles.Rdata')



#####################################################################
# LOTox profiles for FAKE DATA
#-------------------------------------------------------------------
# Un-comment the next lines if you want to save the LOTox profiles
# obtained using 'fake_dataset.Rdata':
# rm(list=setdiff(ls(), c("P_LOTox","RR_LOTox")))
# save.image('fake_data_results/fake_LOTox_profiles.Rdata')
#####################################################################

