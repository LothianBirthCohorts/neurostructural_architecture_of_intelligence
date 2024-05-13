# R script for "Examining the neurostructural architecture of intelligence: The Lothian Birth Cohort 1936 Study"
# Author: Danielle Page

# Load packages -----------------------------------------------------------
#install.packages("foreign")
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("lavaan")
#install.packages("psych")
#install.packages("dplyr")
#install.packages("semPlot")
#install.packages("tidySEM")
#install.packages("Hmisc")
#install.packages("RColorBrewer")
#install.packages("ggpubr")
library(foreign)
library(tidyverse)
library(ggplot2)
library(lavaan)
library(psych)
library(dplyr)
library(semPlot)
library(tidySEM)
library(Hmisc)
library(RColorBrewer)
library(ggpubr)

# Read NeuralArch data ---------------------------------------------------------------
NeuralArch <- data.frame(read.spss("LBC1936_NeuralArchitectureOfIntelligence_DP_01SEP2022.sav", to.data.frame = TRUE))
names(NeuralArch)

# Demographic NeuralArch cleaning --------------------------------------------------

class(NeuralArch$lbc36no)

#sex
class(NeuralArch$sex)
levels(NeuralArch$sex)

#age
class(NeuralArch$agedays_w2)
summary(NeuralArch$agedays_w2)
plot(NeuralArch$agedays_w2)

#years of education
class(NeuralArch$yrsedu_w1)
summary(NeuralArch$yrsedu_w1)
hist(NeuralArch$yrsedu_w1)

#highest occupation (SES)
class(NeuralArch$hmsonum_w1)
summary(NeuralArch$hmsonum_w1)
NeuralArch$hmsonum_w1[NeuralArch$hmsonum_w1 == 3.5] <- 3
hist(NeuralArch$hmsonum_w1)
NeuralArch$hmsonum_w1 <- ordered(NeuralArch$hmsonum_w1, 
                           levels = c(1,2,3,4,5), 
                           labels = c("professional", "managerial", "skilled", "partly skilled", "unskilled"))
summary(NeuralArch$hmsonum_w1)

#father's occupation (SES)
class(NeuralArch$fathclass_w1)
summary(NeuralArch$fathclass_w1)
hist(NeuralArch$fathclass_w1)
NeuralArch$fathclass_w1 <- ordered(NeuralArch$fathclass_w1, 
                             levels = c(1,2,3,4,5), 
                             labels = c("professional", "managerial", "skilled", "partly skilled", "unskilled"))
summary(NeuralArch$fathclass_w1)

#alcohol units/week
class(NeuralArch$alcunitwk_w2)
summary(NeuralArch$alcunitwk_w2)
levels(NeuralArch$alcunitwk_w2)
NeuralArch$alcunitwk_w2[NeuralArch$alcunitwk_w2 == "very infrequently (ie unable to code)"] <- '0'
NeuralArch$alcunitwk_w2 <- droplevels(NeuralArch$alcunitwk_w2)
NeuralArch$alcunitwk_w2 <- as.numeric(as.character(NeuralArch$alcunitwk_w2))
hist(NeuralArch$alcunitwk_w2)
describe(NeuralArch$alcunitwk_w2)
which(NeuralArch$alcunitwk_w2 > (mean(NeuralArch$alcunitwk_w2, na.rm=T) + (3*sd(NeuralArch$alcunitwk_w2, na.rm=T))))
#outlier of 120 units/week - case 330. More than 5 SDs from mean

#smoker status w1
class(NeuralArch$smokcat_w1)
summary(NeuralArch$smokcat_w1)
plot(NeuralArch$smokcat_w1)
class(NeuralArch$smokever_w1)
summary(NeuralArch$smokever_w1)
plot(NeuralArch$smokever_w1)
class(NeuralArch$smokagestart_w1)
summary(NeuralArch$smokagestart_w1)
hist(NeuralArch$smokagestart_w1)
class(NeuralArch$smokagestop_w1)
summary(NeuralArch$smokagestop_w1)
hist(NeuralArch$smokagestop_w1)
class(NeuralArch$smoknumcigs_w1)
summary(NeuralArch$smoknumcigs_w1)
hist(NeuralArch$smoknumcigs_w1)
which(NeuralArch$smoknumcigs_w1 > 75)

#smoker status w2
class(NeuralArch$smokcurr_w2)
summary(NeuralArch$smokcurr_w2)
plot(NeuralArch$smokcurr_w2)
class(NeuralArch$smokprev_w2)
summary(NeuralArch$smokprev_w2)
plot(NeuralArch$smokprev_w2)
class(NeuralArch$smokagestop_w1)
summary(NeuralArch$smokagestop_w1)
hist(NeuralArch$smokagestop_w1)
class(NeuralArch$smoknumcigs_w2)
summary(NeuralArch$smoknumcigs_w2)
hist(NeuralArch$smoknumcigs_w2)

#high blood pressure
class(NeuralArch$hibp_w2)
summary(NeuralArch$hibp_w2)
plot(NeuralArch$hibp_w2)
summary(NeuralArch$hibpcom_w2)

#diabetes
class(NeuralArch$diab_w2)
summary(NeuralArch$diab_w2)
plot(NeuralArch$diab_w2)
summary(NeuralArch$diabcom_w2)

#high cholesterol
class(NeuralArch$hichol_w2)
summary(NeuralArch$hichol_w2)
plot(NeuralArch$hichol_w2)
summary(NeuralArch$hicholcom_w2)

#cardiovascular disease
class(NeuralArch$cvdhist_w2)
summary(NeuralArch$cvdhist_w2)
plot(NeuralArch$cvdhist_w2)
     
#cardiovascular disease
class(NeuralArch$dement_w2)
summary(NeuralArch$dement_w2) #exclude N=2 with dementia? query
plot(NeuralArch$dement_w2)

# Cognitive data cleaning -------------------------------------------------

#mmse
class(NeuralArch$mmse_w2)
describe(NeuralArch$mmse_w2)
hist(NeuralArch$mmse_w2)

#MEMORY
#logical memory story units 1 & 2
class(NeuralArch$lm1_re_w2)
describe(NeuralArch$lm1_re_w2)
hist(NeuralArch$lm1_re_w2)
class(NeuralArch$lm2_re_w2)
describe(NeuralArch$lm2_re_w2)
hist(NeuralArch$lm2_re_w2)
NeuralArch$lmtot_w2 <- (NeuralArch$lm1_re_w2 + NeuralArch$lm2_re_w2)
class(NeuralArch$lmtot_w2)
describe(NeuralArch$lmtot_w2)
hist(NeuralArch$lmtot_w2)

#verbal pairs 1 & 2
class(NeuralArch$vpa1_w2)
describe(NeuralArch$vpa1_w2)
hist(NeuralArch$vpa1_w2, breaks = 20)
class(NeuralArch$vpa2_w2)
describe(NeuralArch$vpa2_w2)
hist(NeuralArch$vpa2_w2, breaks = 8)
NeuralArch$vpatot_w2 <- (NeuralArch$vpa1_w2 + NeuralArch$vpa2_w2)
class(NeuralArch$vpatot_w2)
describe(NeuralArch$vpatot_w2)
hist(NeuralArch$vpatot_w2, breaks = 8)

#digit span
class(NeuralArch$digback_w2)
describe(NeuralArch$digback_w2)
hist(NeuralArch$digback_w2)

#CRYSTALLISED
#verbal fluency C, F, L, total
class(NeuralArch$vfc_w2)
describe(NeuralArch$vfc_w2)
hist(NeuralArch$vfc_w2)
class(NeuralArch$vff_w2)
describe(NeuralArch$vff_w2)
hist(NeuralArch$vff_w2)
class(NeuralArch$vfl_w2)
describe(NeuralArch$vfl_w2)
hist(NeuralArch$vfl_w2)
class(NeuralArch$vftot_w2)
describe(NeuralArch$vftot_w2)
hist(NeuralArch$vftot_w2)

#nart
class(NeuralArch$nart_w2)
describe(NeuralArch$nart_w2)
hist(NeuralArch$nart_w2, breaks=50)

#wtar
class(NeuralArch$wtar_w2)
describe(NeuralArch$wtar_w2)
hist(NeuralArch$wtar_w2, breaks=50)

#PROCESSING SPEED
#symbol search
class(NeuralArch$symsear_w2)
describe(NeuralArch$symsear_w2)
hist(NeuralArch$symsear_w2, breaks=20)

#digit symbol
class(NeuralArch$digsym_w2)
describe(NeuralArch$digsym_w2)
hist(NeuralArch$digsym_w2)

#inspection time
class(NeuralArch$ittotal_w2)
describe(NeuralArch$ittotal_w2)
hist(NeuralArch$ittotal_w2, breaks=50)

#reaction time
class(NeuralArch$crtmean_w2)
describe(NeuralArch$crtmean_w2)
hist(NeuralArch$crtmean_w2, breaks=50)

#VISUOSPATIAL
#block design
class(NeuralArch$blkdes_w2)
describe(NeuralArch$blkdes_w2)
hist(NeuralArch$blkdes_w2, breaks=20)

#matrix reasoning
class(NeuralArch$matreas_w2)
describe(NeuralArch$matreas_w2)
hist(NeuralArch$matreas_w2)

#spatial span forward, back, total
class(NeuralArch$spanf_w2)
describe(NeuralArch$spanf_w2)
hist(NeuralArch$spanf_w2)
class(NeuralArch$spanb_w2)
describe(NeuralArch$spanb_w2)
hist(NeuralArch$spanb_w2)
class(NeuralArch$spantot_w2)
describe(NeuralArch$spantot_w2)
hist(NeuralArch$spantot_w2, breaks = 20)

# MRI variable cleaning ---------------------------------------------------
class(NeuralArch$ageMRI_w2)
describe(NeuralArch$ageMRI_w2)
hist(NeuralArch$ageMRI_w2, breaks=50)
NeuralArch$MRI_comments_w2

#stroke mask
class(NeuralArch$stroke_mask_w2)
summary(NeuralArch$stroke_mask_w2)
plot(NeuralArch$stroke_mask_w2)

#Intra cranial volume
class(NeuralArch$ICV_mm3_wX) #in mm3
describe(NeuralArch$ICV_mm3_wX)
hist(NeuralArch$ICV_mm3_wX, breaks=50)
#Last wave used for ICV - NeuralArch$ICV_wave_wXn

#total brain volume
class(NeuralArch$brain_mm3_w2) #raw mm3 measure
describe(NeuralArch$brain_mm3_w2)
hist(NeuralArch$brain_mm3_w2, breaks=50)
class(NeuralArch$brainIcv_ratio_w2) #as ratio of ICV volume
describe(NeuralArch$brainIcv_ratio_w2)
hist(NeuralArch$brainIcv_ratio_w2, breaks=50)

#grey matter
class(NeuralArch$gm_mm3_w2) #raw mm3 measure
describe(NeuralArch$gm_mm3_w2)
hist(NeuralArch$gm_mm3_w2, breaks=50)
class(NeuralArch$gmIcv_ratio_w2) #as ratio of ICV volume
describe(NeuralArch$gmIcv_ratio_w2)
hist(NeuralArch$gmIcv_ratio_w2, breaks=50)

#normal appearing white matter
class(NeuralArch$nawm_mm3_w2) #raw mm3 measure
describe(NeuralArch$nawm_mm3_w2)
hist(NeuralArch$nawm_mm3_w2, breaks=50)
class(NeuralArch$nawmIcv_ratio_w2) #as ratio of ICV volume
describe(NeuralArch$nawmIcv_ratio_w2)
hist(NeuralArch$nawmIcv_ratio_w2, breaks=50)

#White matter hyperintensities
class(NeuralArch$wmh_mm3_w2) #raw mm3 measure
describe(NeuralArch$wmh_mm3_w2)
hist(NeuralArch$wmh_mm3_w2, breaks=50)
class(NeuralArch$wmhIcv_ratio_w2) #as ratio of ICV volume
describe(NeuralArch$wmhIcv_ratio_w2)
hist(NeuralArch$wmhIcv_ratio_w2, breaks=50)

#WHITE MATTER FA/MD

#Genu
class(NeuralArch$Genu_FA_w2) #Fractional anisotropy
describe(NeuralArch$Genu_FA_w2)
NeuralArch$Genu_FA_w2[NeuralArch$Genu_FA_w2 == -999] <- NA
hist(NeuralArch$Genu_FA_w2, breaks=50)
class(NeuralArch$Genu_MD_w2) #Mean diffusivity
describe(NeuralArch$Genu_MD_w2)
NeuralArch$Genu_MD_w2[NeuralArch$Genu_MD_w2 == -999] <- NA
hist(NeuralArch$Genu_MD_w2, breaks=50)

#Splenium
class(NeuralArch$Splenium_FA_w2) #Fractional anisotropy
describe(NeuralArch$Splenium_FA_w2)
NeuralArch$Splenium_FA_w2[NeuralArch$Splenium_FA_w2 == -999] <- NA
hist(NeuralArch$Splenium_FA_w2, breaks=50)
class(NeuralArch$Splenium_MD_w2) #Mean diffusivity
describe(NeuralArch$Splenium_MD_w2)
NeuralArch$Splenium_MD_w2[NeuralArch$Splenium_MD_w2 == -999] <- NA
hist(NeuralArch$Splenium_MD_w2, breaks=50)

#Arcuate fasciculus (L+R)
class(NeuralArch$LArc_FA_w2) #LEFT Fractional anisotropy
describe(NeuralArch$LArc_FA_w2)
NeuralArch$LArc_FA_w2[NeuralArch$LArc_FA_w2 == -999] <- NA
hist(NeuralArch$LArc_FA_w2, breaks=50)
class(NeuralArch$LArc_MD_w2) #LEFT Mean diffusivity
describe(NeuralArch$LArc_MD_w2)
NeuralArch$LArc_MD_w2[NeuralArch$LArc_MD_w2 == -999] <- NA
hist(NeuralArch$LArc_MD_w2, breaks=50)
class(NeuralArch$RArc_FA_w2) #RIGHT Fractional anisotropy
describe(NeuralArch$RArc_FA_w2)
NeuralArch$RArc_FA_w2[NeuralArch$RArc_FA_w2 == -999] <- NA
hist(NeuralArch$RArc_FA_w2, breaks=50)
class(NeuralArch$RArc_MD_w2) #RIGHT Mean diffusivity
describe(NeuralArch$RArc_MD_w2)
NeuralArch$RArc_MD_w2[NeuralArch$RArc_MD_w2 == -999] <- NA
hist(NeuralArch$RArc_MD_w2, breaks=50)

#Anterior thalamic radiation (L+R)
class(NeuralArch$LATR_FA_w2) #LEFT Fractional anisotropy
describe(NeuralArch$LATR_FA_w2)
NeuralArch$LATR_FA_w2[NeuralArch$LATR_FA_w2 == -999] <- NA
hist(NeuralArch$LATR_FA_w2, breaks=50)
class(NeuralArch$LATR_MD_w2) #LEFT Mean diffusivity
describe(NeuralArch$LATR_MD_w2)
NeuralArch$LATR_MD_w2[NeuralArch$LATR_MD_w2 == -999] <- NA
hist(NeuralArch$LATR_MD_w2, breaks=50)
class(NeuralArch$RATR_FA_w2) #RIGHT Fractional anisotropy
describe(NeuralArch$RATR_FA_w2)
NeuralArch$RATR_FA_w2[NeuralArch$RATR_FA_w2 == -999] <- NA
hist(NeuralArch$RATR_FA_w2, breaks=50)
class(NeuralArch$RATR_MD_w2) #RIGHT Mean diffusivity
describe(NeuralArch$RATR_MD_w2)
NeuralArch$RATR_MD_w2[NeuralArch$RATR_MD_w2 == -999] <- NA
hist(NeuralArch$RATR_MD_w2, breaks=50)

#Cingulum
class(NeuralArch$LCing_FA_w2) #LEFT Fractional anisotropy
describe(NeuralArch$LCing_FA_w2)
NeuralArch$LCing_FA_w2[NeuralArch$LCing_FA_w2 == -999] <- NA
hist(NeuralArch$LCing_FA_w2, breaks=50)
class(NeuralArch$LCing_MD_w2) #LEFT Mean diffusivity
describe(NeuralArch$LCing_MD_w2)
NeuralArch$LCing_MD_w2[NeuralArch$LCing_MD_w2 == -999] <- NA
hist(NeuralArch$LCing_MD_w2, breaks=50)
class(NeuralArch$RCing_FA_w2) #RIGHT Fractional anisotropy
describe(NeuralArch$RCing_FA_w2)
NeuralArch$RCing_FA_w2[NeuralArch$RCing_FA_w2 == -999] <- NA
hist(NeuralArch$RCing_FA_w2, breaks=50)
class(NeuralArch$RCing_MD_w2) #RIGHT Mean diffusivity
describe(NeuralArch$RCing_MD_w2)
NeuralArch$RCing_MD_w2[NeuralArch$RCing_MD_w2 == -999] <- NA
hist(NeuralArch$RCing_MD_w2, breaks=50)

#Uncinate fasciculus
class(NeuralArch$LUnc_FA_w2) #LEFT Fractional anisotropy
describe(NeuralArch$LUnc_FA_w2)
NeuralArch$LUnc_FA_w2[NeuralArch$LUnc_FA_w2 == -999] <- NA
hist(NeuralArch$LUnc_FA_w2, breaks=50)
class(NeuralArch$LUnc_MD_w2) #LEFT Mean diffusivity
describe(NeuralArch$LUnc_MD_w2)
NeuralArch$LUnc_MD_w2[NeuralArch$LUnc_MD_w2 == -999] <- NA
hist(NeuralArch$LUnc_MD_w2, breaks=50)
class(NeuralArch$RUnc_FA_w2) #RIGHT Fractional anisotropy
describe(NeuralArch$RUnc_FA_w2)
NeuralArch$RUnc_FA_w2[NeuralArch$RUnc_FA_w2 == -999] <- NA
hist(NeuralArch$RUnc_FA_w2, breaks=50)
class(NeuralArch$RUnc_MD_w2) #RIGHT Mean diffusivity
describe(NeuralArch$RUnc_MD_w2)
NeuralArch$RUnc_MD_w2[NeuralArch$RUnc_MD_w2 == -999] <- NA
hist(NeuralArch$RUnc_MD_w2, breaks=50)

#Inferior longitudinal fasciculus
class(NeuralArch$LILF_FA_w2) #LEFT Fractional anisotropy
describe(NeuralArch$LILF_FA_w2)
NeuralArch$LILF_FA_w2[NeuralArch$LILF_FA_w2 == -999] <- NA
hist(NeuralArch$LILF_FA_w2, breaks=50)
class(NeuralArch$LILF_MD_w2) #LEFT Mean diffusivity
describe(NeuralArch$LILF_MD_w2)
NeuralArch$LILF_MD_w2[NeuralArch$LILF_MD_w2 == -999] <- NA
hist(NeuralArch$LILF_MD_w2, breaks=50)
class(NeuralArch$RILF_FA_w2) #RIGHT Fractional anisotropy
describe(NeuralArch$RILF_FA_w2)
NeuralArch$RILF_FA_w2[NeuralArch$RILF_FA_w2 == -999] <- NA
hist(NeuralArch$RILF_FA_w2, breaks=50)
class(NeuralArch$RILF_MD_w2) #RIGHT Mean diffusivity
describe(NeuralArch$RILF_MD_w2)
NeuralArch$RILF_MD_w2[NeuralArch$RILF_MD_w2 == -999] <- NA
hist(NeuralArch$RILF_MD_w2, breaks=50)

# Create analytic sample -------------------------------------
summary(NeuralArch$attend73) #225/1091 did not attend W2
data <- NeuralArch[which(NeuralArch$attend73 == 'Yes'),] #exclude those who didn't attend at W2 (remaining N=866)
summary(NeuralArch$ageMRI_w2) #further 134/866 did not have scan
data <- NeuralArch[!is.na(NeuralArch$ageMRI_w2),] #exclude those who didn't have scan (remaining N=731)

#Check MRI comments for missing data due to terminated scans
data$MRI_comments_w2[data$MRI_comments_w2 =="                                                                           "] <- NA
data$MRI_comments_w2[!is.na(data$MRI_comments_w2)]
subset(data, subset=!is.na(data$MRI_comments_w2), select=c(lbc36no, MRI_comments_w2, brain_mm3_w2, gm_mm3_w2, Genu_FA_w2, RArc_FA_w2))
#remove LBC numbers: 5, 65, 129, 156, 200, 206, 302, 332, 335, 481, 543, 552, 682, 694, 776, 793, 811, 837, 870, 930, 954, 1002, 1021, 1041, 1053, 1125, 1166, 1196, 1214, 1240, 1258, 1259, 1260, 1261

#385/1133 - note says MRI aborted but has MRI values
#688/769/787/983/1076/1099/1192/1222/1296 MRI terminated but have values for WM
data <- data[-c(4,39,75,91,122,126,192,205,208, 285,315,321,389, 395, 441, 453, 463,475,495,527,537, 564,574,587, 597, 630, 651,669,679,692,698,699,700,701),]

#complete.cases(data[, c(9:28, 49:83)]) #complete imaging & cognitive  logical vector


# Descriptives ------------------------------------------------------------

describe(data$agedays_w2/365.25) #mean 72.5y, sd 0.71 (range 70.91-74.11)
describe(data$ageMRI_w2/365.25) #(mean at MRI - 72.7y, sd 0.73)
summary(data$sex) #364 female 333 male
describe(data$yrsedu_w1) #mean 10.79 years edu, sd 1.14 (range 8-14)
summary(data$hmsonum_w1) #138 prof, 257 manag, 261 skill, 24 partskill, 4 unskill, 13 NA
summary(data$fathclass_w1) #42 prof, 124 manag, 359 skill, 69 partskill, 43 unskill, 60 NA

#days between cog & MRI
mean(data$ageMRI_w2 - data$agedays_w2)
sd(data$ageMRI_w2 - data$agedays_w2)

#health & lifestyle:
summary(data[29:46])
summary(data[30:37]) #343 hibp, 75 diab, 289 hichol, 189 cvdhist, 0 dementia

#cognitive
describe(data[9:28])
describe(data[84:85])
colnames(data)

#MRI
describe(data[52:55])
describe(data[60:83])
  mean(data$RILF_MD_w2 , na.rm=T)

# CFA: Extract latent domains & g -----------------------------------------------------

#create temp data frame
tmp = data.frame(it = data$ittotal_w2, digsym = data$digsym_w2, sym = data$symsear_w2, crt = (data$crtmean_w2)*100,
                 nart = data$nart_w2, wtar = data$wtar_w2, vf = data$vftot_w2,
                 lm = data$lmtot_w2, vpa = data$vpatot_w2, digsp = data$digback_w2, 
                 blk = data$blkdes_w2, mat = data$matreas_w2, sspn = data$spantot_w2)

### PROCESSING SPEED
#Speed latent variable predicted by ittotal, digsym, symsear, crt
speed <- 'speed =~ it + digsym + sym + crt'
#fit model
fit<-cfa(speed, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
parTable(fit)
varTable(fit)
#vcov(fit)
FactorScores <- predict(fit, tmp)
data <- cbind(data, FactorScores)
hist(data$speed)
((0.485^2) + (0.827^2) + (0.749^2) + (0.673^2)) / 4 #variance

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##CRYSTALLISED INTELLIGENCE
#Crystallised latent variable predicted by nart, wtar, vftotal
cryst <- 'cryst =~ nart + wtar + vf'
#fit model
fit<-cfa(cryst, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
parTable(fit)
varTable(fit)
cor.test(data$nart,data$wtar)
#vcov(fit)
FactorScores <- predict(fit, tmp)
data <- cbind(data, FactorScores)
hist(data$cryst)
((0.951^2) + (0.951^2) + (0.449^2)) / 3 #variance


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##MEMORY 
#Memory latent variable predicted by vpa, lm, digback
memory <- 'memory =~ lm + vpa + digsp'
#fit model
fit<-cfa(memory, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
parTable(fit)
varTable(fit)
#vcov(fit)
FactorScores <- predict(fit, tmp)
data <- cbind(data, FactorScores)
hist(data$memory)
((0.804^2) + (0.672^2) + (0.384^2)) / 3 #variance

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##VISUOSPATIAL  
#Visuospatial latent variable predicted by block, matrix, spatial span
visuo <- 'visuo =~ blk + mat + sspn'
#fit model
fit<-cfa(visuo, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
parTable(fit)
varTable(fit)
#vcov(fit)
FactorScores <- predict(fit, tmp)
data <- cbind(data, FactorScores)
hist(data$visuo)
((0.803^2) + (0.663^2) + (0.539^2)) / 4 #variance

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## G FACTOR  
#g  latent variable predicted by 13 items
g <- 'g =~ it + digsym + sym + crt + nart + wtar + vf + lm + vpa + digsp + blk + mat + sspn

it~~digsym
it~~sym
it~~crt
digsym~~sym
digsym~~crt
sym~~crt

nart~~wtar
vf~~nart
vf~~wtar

blk~~mat
blk~~sspn
mat~~sspn

lm~~vpa
lm~~digsp
vpa~~digsp'

#fit model
fit<-cfa(g, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
parTable(fit)
varTable(fit)
#vcov(fit)
FactorScores <- predict(fit, tmp)
data <- cbind(data, FactorScores)
hist(data$g)
plot(data$g)
#((0.290^2)+(0.523^2)+(0.487^2)+(0.333^2)+(0.923^2)+(0.923^2)+(0.492^2)+(0.526^2)+(0.465^2)+(0.462^2)+(0.485^2)+(0.484^2)+(0.283^2)) / 13 #variance BEFORE CORRECTION
((0.378^2)+(0.620^2)+(0.573^2)+(0.422^2)+(0.705^2)+(0.707^2)+(0.583^2)+(0.579^2)+(0.516^2)+(0.545^2)+(0.550^2)+(0.593^2)+(0.380^2)) / 13 #variance

#Visualising with semPlot
semPaths(fit, what = "std", style = "lisrel",
         residScale = 6, theme = "colorblind", nCharNodes = 0,
         reorder = F, legend.cex = 0.5, nDigits = 2,
         rotation = 2, layout = "tree2", curvePivot = T,
         sizeMan = 5, sizeLat = 10, edge.label.cex = .8,
         fade=F, trans=F)


# CFA: Heirarchical model (for demo) --------------------------------------


heirarchical <- 'g =~ speed + cryst + memory + visuo
speed =~ it + digsym + sym + crt
cryst =~ nart + wtar + vf
memory =~ lm + vpa + digsp
visuo =~ blk + mat + sspn
nart~~wtar'
fit<-cfa(heirarchical, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
parTable(fit)
varTable(fit)

# CFA: Extract latent gFA & gMD ----------------------------------------------------------------

#Calculate average across left & right tracts
data$Arc_FA_w2 <- (data$LArc_FA_w2 + data$RArc_FA_w2)/2
data$Arc_MD_w2 <- (data$LArc_MD_w2 + data$RArc_MD_w2)/2
data$ATR_FA_w2 <- (data$LATR_FA_w2 + data$RATR_FA_w2)/2
data$ATR_MD_w2 <- (data$LATR_MD_w2 + data$RATR_MD_w2)/2
data$Cing_FA_w2 <- (data$LCing_FA_w2 + data$RCing_FA_w2)/2
data$Cing_MD_w2 <- (data$LCing_MD_w2 + data$RCing_MD_w2)/2
data$Unc_FA_w2 <- (data$LUnc_FA_w2 + data$RUnc_FA_w2)/2
data$Unc_MD_w2 <- (data$LUnc_MD_w2 + data$RUnc_MD_w2)/2
data$ILF_FA_w2 <- (data$LILF_FA_w2 + data$RILF_FA_w2)/2
data$ILF_MD_w2 <- (data$LILF_MD_w2 + data$RILF_MD_w2)/2


tmp = data.frame(genufa = data$Genu_FA_w2, genumd = data$Genu_MD_w2,
                 splenfa = data$Splenium_FA_w2, splenmd = data$Splenium_MD_w2,
                 arcfa = data$Arc_FA_w2, arcmd = data$Arc_MD_w2,
                 atrfa = data$ATR_FA_w2 , atrmd = data$ATR_MD_w2,
                 cingfa = data$Cing_FA_w2, cingmd = data$Cing_MD_w2,
                 uncfa = data$Unc_FA_w2, uncmd = data$Unc_MD_w2,
                 ilffa = data$ILF_FA_w2, ilfmd = data$ILF_MD_w2)

#latent variable for FA
fa <- 'fa =~ genufa + splenfa + arcfa + atrfa + cingfa + uncfa + ilffa'
#fit model
fit<-cfa(fa, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
parTable(fit)
varTable(fit)
#vcov(fit)
FA_value <- predict(fit, tmp)
data <- cbind(data, FA_value)
hist(data$fa)

(0.621^2 + 0.414^2 + 0.671^2 + 0.682^2 + 0.615^2 + 0.742^2 + 0.508^2) / 7


#latent variable for MD
md <- 'md =~ genumd + splenmd + arcmd + atrmd + cingmd + uncmd + ilfmd'
#fit model
fit<-cfa(md, missing="fiml.x", data=tmp)
summary(fit, fit.measures=T, standardize=T, rsquare=T, modindices = TRUE)
parTable(fit)
varTable(fit)
vcov(fit)
MD_value <- predict(fit, tmp)
data <- cbind(data, MD_value)
hist(data$md)

(0.610^2 + 0.284^2 + 0.693^2 + 0.762^2 + 0.723^2 + 0.730^2 + 0.437^2) / 7

#Descriptives
describe(data[101:102])

# Regressions: Domains/g =~ TBV,GM,NAWM,WMH,gFA,gMD ------------------------------------

#sanity checks for unusual beta values when using brain mm3/icv mm3 rather than brain/icv ratios
#summary(data$brain_mm3_w2)
#sd(data$brain_mm3_w2, na.rm=T)
#summary(data$ICV_mm3_wX)
#sd(data$ICV_mm3_wX, na.rm=T)
#summary(data$brainIcv_ratio_w2)
#sd(data$brainIcv_ratio_w2, na.rm=T)
#summary(data$brain_mm3_w2/data$ICV_mm3_wX)
#sd(data$brain_mm3_w2/data$ICV_mm3_wX, na.rm=T)
#sanitycheck <- corr.test(data.frame(ratio = data$brainIcv_ratio_w2, manualcalc = (data$brain_mm3_w2/data$ICV_mm3_wX)))
#sanitycheck
#cor.test(data$speed, data$g)
#plot(data$g)
#hist(data$g)

#library(ltm)
#x<-data.frame(data$g, data$speed, data$memory, data$cryst, data$visuo, data$brain_mm3_w2, data$ICV_mm3_wX, data$brainIcv_ratio_w2)
#para=ltm::rcor.test(as.matrix(x), method="pearson", use="pairwise.complete.obs")
#para$cor.mat
##SOLUTION: high correlation between brain/icv was causing inflated variance which pushed beta up. true value is lower.
##calculate equivalent to betas by creating a residual of brain_mm3 on all covariates and run coorrelation with domain
#model = (lm(scale(data$speed) ~ scale(data$brain_mm3_w2) + scale(data$ICV_mm3_wX) + scale(data$agedays_w2) + as.numeric(data$sex)))
#car::vif(model)

##TBV##
data$x <- resid(lm(data$brain_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$speed, data$x)
cor.test(data$cryst, data$x)
cor.test(data$memory, data$x)
cor.test(data$visuo, data$x)
cor.test(data$g, data$x)
##GM## 
data$x <- resid(lm(data$gm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$speed, data$x)
cor.test(data$cryst, data$x)
cor.test(data$memory, data$x)
cor.test(data$visuo, data$x)
cor.test(data$g, data$x)
##NAWM##
data$x <- resid(lm(data$nawm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$speed, data$x)
cor.test(data$cryst, data$x)
cor.test(data$memory, data$x)
cor.test(data$visuo, data$x)
cor.test(data$g, data$x)
##WMH## 
data$x <- resid(lm(data$wmh_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$speed, data$x)
cor.test(data$cryst, data$x)
cor.test(data$memory, data$x)
cor.test(data$visuo, data$x)
cor.test(data$g, data$x)
##FA
data$x <- resid(lm(data$fa ~ data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$speed, data$x)
cor.test(data$cryst, data$x)
cor.test(data$memory, data$x)
cor.test(data$visuo, data$x)
cor.test(data$g, data$x)
##MD 
data$x <- resid(lm(data$md ~ data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$speed, data$x)
cor.test(data$cryst, data$x)
cor.test(data$memory, data$x)
cor.test(data$visuo, data$x)
cor.test(data$g, data$x)

pvals1 <-c(.000, .000, .000, .000, .000, 
           .000, .000, .000, .000, .000, 
           .000, .009, .006, .000, .000, 
           .000, .208, .036, .006, .000, 
           .000, .115, .076, .011, .000, 
           .000, .084, .028, .009, .006)
p.adjust(pvals1, method="fdr")


# CFA: Bi-factor model of domains & g ---------------------------------------------------------

#create temp data frame
tmp = data.frame(it = data$ittotal_w2, digsym = data$digsym_w2, sym = data$symsear_w2, crt = (data$crtmean_w2)*100,
                 nart = data$nart_w2, wtar = data$wtar_w2, vf = data$vftot_w2,
                 lm = data$lmtot_w2, vpa = data$vpatot_w2, digsp = data$digback_w2, 
                 blk = data$blkdes_w2, mat = data$matreas_w2, sspn = data$spantot_w2)

## BI-FACTOR DOMAINS & G  
bifac <- 'bi_g =~ it + digsym + sym + crt + nart + wtar + vf + lm + vpa + digsp + blk + mat + sspn
bi_speed =~ it + digsym + sym + crt
bi_cryst =~ nart + wtar + vf
bi_memory =~ lm + vpa + digsp
bi_visuo =~ blk + mat + sspn
vpa ~~ 0*vpa' #fix variance of vpa to 0

#bifac <- 'bi_g =~ 1*it + digsym + sym + crt + nart + wtar + vf + lm + vpa + digsp + blk + mat + sspn
#bi_speed =~ 1*it + digsym + sym + crt
#bi_cryst =~ 1*nart + wtar + vf
#bi_memory =~ 1*lm + vpa + digsp
#bi_visuo =~ 1*blk + mat + sspn
#vpa ~~ 0*vpa' #fix variance of vpa to 0
#THIS CODE FIXES FIRST LOADING TO 1 AS IT SHOULD HAVE BEEN TO START

#fit model
fit<-cfa(bifac, missing="fiml.x", data=tmp, orthogonal=T, std.lv=T)
vcov(fit)
parTable(fit)
varTable(fit)
summary(fit, fit.measures=T, standardized=T, rsquare=T, modindices=T)
fitmeasures(fit, c('cfi', 'tli', 'rmsea', 'srmr'))
standardizedsolution(fit)
residuals(fit, type="cor")
FactorScores <- predict(fit, tmp)
head(FactorScores)
data <- cbind(data, FactorScores)
cor(data$g, data$bi_g)
standardizedsolution(fit)
((0.387^2)+(0.640^2)+(0.610^2)+(0.436^2)+(0.664^2)+(0.670^2)+(0.507^2)+(0.552^2)+(0.513^2)+(0.565^2)+(0.591^2)+(0.587^2)+(0.447^2))/13 #variance
((0.303^2)+(0.515^2)+(0.439^2)+(0.538^2))/4 #variance
((0.688^2)+(0.668^2)+(0.131^2))/3 #variance
((0.299^2)+(0.858^2)+(0.036^2))/3 #variance
((0.573^2)+(0.323^2)+(0.293^2))/3 #variance

# Correlation matrix ------------------------------------------------------
tmp = data.frame(it = data$ittotal_w2, digsym = data$digsym_w2, sym = data$symsear_w2, crt = (data$crtmean_w2)*100,
                 nart = data$nart_w2, wtar = data$wtar_w2, vf = data$vftot_w2,
                 lm = data$lmtot_w2, vpa = data$vpatot_w2, digsp = data$digback_w2, 
                 blk = data$blkdes_w2, mat = data$matreas_w2, sspn = data$spantot_w2)
corrmat <- corr.test(tmp)
write.csv(corrmat$r, 'corrmat_r_cogitems.csv')
write.csv(corrmat$p, 'corrmat_p_cogitems.csv')

tmp = data.frame(speed = data$speed, cryst = data$cryst, memory = data$memory,
                 visuo = data$visuo, g = data$g, 
                 bi_speed = data$bi_speed, bi_cryst = data$bi_cryst, bi_memory = data$bi_memory,
                 bi_visuo = data$bi_visuo, bi_g = data$bi_g,
                 TBV = data$brainIcv_ratio_w2, GM = data$gmIcv_ratio_w2, NAWM = data$nawmIcv_ratio_w2, 
                 WMH = data$wmhIcv_ratio_w2, FA = data$fa, MD = data$md)
corrmat <- corr.test(tmp)
write.csv(corrmat$r, 'corrmat_r_cogbrain.csv')
write.csv(corrmat$p, 'corrmat_p_cogbrain.csv')




# Regressions: Domains/g =~ TBV,GM,NAWM,WMH,gFA,gMD -------------------------------------

##TBV##
data$x <- resid(lm(data$brain_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$bi_speed, data$x)
cor.test(data$bi_cryst, data$x)
cor.test(data$bi_memory, data$x)
cor.test(data$bi_visuo, data$x)
cor.test(data$bi_g, data$x)
##GM##
data$x <- resid(lm(data$gm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$bi_speed, data$x)
cor.test(data$bi_cryst, data$x)
cor.test(data$bi_memory, data$x)
cor.test(data$bi_visuo, data$x)
cor.test(data$bi_g, data$x)
##NAWM##
data$x <- resid(lm(data$nawm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$bi_speed, data$x)
cor.test(data$bi_cryst, data$x)
cor.test(data$bi_memory, data$x)
cor.test(data$bi_visuo, data$x)
cor.test(data$bi_g, data$x)
##WMH##
data$x <- resid(lm(data$wmh_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$bi_speed, data$x)
cor.test(data$bi_cryst, data$x)
cor.test(data$bi_memory, data$x)
cor.test(data$bi_visuo, data$x)
cor.test(data$bi_g, data$x)
##FA##
data$x <- resid(lm(data$fa ~ data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$bi_speed, data$x)
cor.test(data$bi_cryst, data$x)
cor.test(data$bi_memory, data$x)
cor.test(data$bi_visuo, data$x)
cor.test(data$bi_g, data$x)
##MD##
data$x <- resid(lm(data$md ~ data$agedays_w2 + data$sex, na.action=na.exclude))
cor.test(data$bi_speed, data$x)
cor.test(data$bi_cryst, data$x)
cor.test(data$bi_memory, data$x)
cor.test(data$bi_visuo, data$x)
cor.test(data$bi_g, data$x)

pvals1 <-c(.000, .876, .234, .008, .000, 
           .003, .895, .861, .008, .000, 
           .000, .198, .005, .263, .000, 
           .000, .094, .207, .452, .000, 
           .000, .417, .318, .371, .001, 
           .007, .005, .011, .139, .006)
p.adjust(pvals1, method="fdr")

# Figure 2: Create model coeff variables ---------------

##INIDIVIDUAL LM VARIABLES##
data$x <- resid(lm(data$brain_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
speed_tbv<-cor.test(data$speed, data$x)
cryst_tbv<-cor.test(data$cryst, data$x)
memory_tbv<-cor.test(data$memory, data$x)
visuo_tbv<-cor.test(data$visuo, data$x)
g_tbv<-cor.test(data$g, data$x)
data$x <- resid(lm(data$gm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
speed_gm<-cor.test(data$speed, data$x)
cryst_gm<-cor.test(data$cryst, data$x)
memory_gm<-cor.test(data$memory, data$x)
visuo_gm<-cor.test(data$visuo, data$x)
g_gm<-cor.test(data$g, data$x)
data$x <- resid(lm(data$nawm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
speed_nawm<-cor.test(data$speed, data$x)
cryst_nawm<-cor.test(data$cryst, data$x)
memory_nawm<-cor.test(data$memory, data$x)
visuo_nawm<-cor.test(data$visuo, data$x)
g_nawm<-cor.test(data$g, data$x)
data$x <- resid(lm(data$wmh_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
speed_wmh<-cor.test(data$speed, data$x)
cryst_wmh<-cor.test(data$cryst, data$x)
memory_wmh<-cor.test(data$memory, data$x)
visuo_wmh<-cor.test(data$visuo, data$x)
g_wmh<-cor.test(data$g, data$x)
data$x <- resid(lm(data$fa ~ data$agedays_w2 + data$sex, na.action=na.exclude))
speed_fa<-cor.test(data$speed, data$x)
cryst_fa<-cor.test(data$cryst, data$x)
memory_fa<-cor.test(data$memory, data$x)
visuo_fa<-cor.test(data$visuo, data$x)
g_fa<-cor.test(data$g, data$x)
data$x <- resid(lm(data$md ~ data$agedays_w2 + data$sex, na.action=na.exclude))
speed_md<-cor.test(data$speed, data$x)
cryst_md<-cor.test(data$cryst, data$x)
memory_md<-cor.test(data$memory, data$x)
visuo_md<-cor.test(data$visuo, data$x)
g_md<-cor.test(data$g, data$x)

##BIFACTOR LM VARIABLES##
data$x <- resid(lm(data$brain_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
bi_speed_tbv<-cor.test(data$bi_speed, data$x)
bi_cryst_tbv<-cor.test(data$bi_cryst, data$x)
bi_memory_tbv<-cor.test(data$bi_memory, data$x)
bi_visuo_tbv<-cor.test(data$bi_visuo, data$x)
bi_g_tbv<-cor.test(data$bi_g, data$x)
data$x <- resid(lm(data$gm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
bi_speed_gm<-cor.test(data$bi_speed, data$x)
bi_cryst_gm<-cor.test(data$bi_cryst, data$x)
bi_memory_gm<-cor.test(data$bi_memory, data$x)
bi_visuo_gm<-cor.test(data$bi_visuo, data$x)
bi_g_gm<-cor.test(data$bi_g, data$x)
data$x <- resid(lm(data$nawm_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
bi_speed_nawm<-cor.test(data$bi_speed, data$x)
bi_cryst_nawm<-cor.test(data$bi_cryst, data$x)
bi_memory_nawm<-cor.test(data$bi_memory, data$x)
bi_visuo_nawm<-cor.test(data$bi_visuo, data$x)
bi_g_nawm<-cor.test(data$bi_g, data$x)
data$x <- resid(lm(data$wmh_mm3_w2 ~ data$ICV_mm3_wX + data$agedays_w2 + data$sex, na.action=na.exclude))
bi_speed_wmh<-cor.test(data$bi_speed, data$x)
bi_cryst_wmh<-cor.test(data$bi_cryst, data$x)
bi_memory_wmh<-cor.test(data$bi_memory, data$x)
bi_visuo_wmh<-cor.test(data$bi_visuo, data$x)
bi_g_wmh<-cor.test(data$bi_g, data$x)
data$x <- resid(lm(data$fa ~ data$agedays_w2 + data$sex, na.action=na.exclude))
bi_speed_fa<-cor.test(data$bi_speed, data$x)
bi_cryst_fa<-cor.test(data$bi_cryst, data$x)
bi_memory_fa<-cor.test(data$bi_memory, data$x)
bi_visuo_fa<-cor.test(data$bi_visuo, data$x)
bi_g_fa<-cor.test(data$bi_g, data$x)
data$x <- resid(lm(data$md ~ data$agedays_w2 + data$sex, na.action=na.exclude))
bi_speed_md<-cor.test(data$bi_speed, data$x)
bi_cryst_md<-cor.test(data$bi_cryst, data$x)
bi_memory_md<-cor.test(data$bi_memory, data$x)
bi_visuo_md<-cor.test(data$bi_visuo, data$x)
bi_g_md<-cor.test(data$bi_g, data$x)

# Figure 2A: Model comparison: G effect size --------------------------------

#g bind coeff tables
g_tbv_r<-g_tbv$estimate
g_gm_r<-g_gm$estimate
g_nawm_r<-g_nawm$estimate
g_wmh_r<-g_wmh$estimate
g_fa_r<-g_fa$estimate
g_md_r<-g_md$estimate
bi_g_tbv_r<-bi_g_tbv$estimate
bi_g_gm_r<-bi_g_gm$estimate
bi_g_nawm_r<-bi_g_nawm$estimate
bi_g_wmh_r<-bi_g_wmh$estimate
bi_g_fa_r<-bi_g_fa$estimate
bi_g_md_r<-bi_g_md$estimate
g_tbv_p<-g_tbv$p.value
g_gm_p<-g_gm$p.value
g_nawm_p<-g_nawm$p.value
g_wmh_p<-g_wmh$p.value
g_fa_p<-g_fa$p.value
g_md_p<-g_md$p.value
bi_g_tbv_p<-bi_g_tbv$p.value
bi_g_gm_p<-bi_g_gm$p.value
bi_g_nawm_p<-bi_g_nawm$p.value
bi_g_wmh_p<-bi_g_wmh$p.value
bi_g_fa_p<-bi_g_fa$p.value
bi_g_md_p<-bi_g_md$p.value
g_tbv_conf1<-g_tbv$conf.int[1]
g_gm_conf1<-g_gm$conf.int[1]
g_nawm_conf1<-g_nawm$conf.int[1]
g_wmh_conf1<-g_wmh$conf.int[1]
g_fa_conf1<-g_fa$conf.int[1]
g_md_conf1<-g_md$conf.int[1]
bi_g_tbv_conf1<-bi_g_tbv$conf.int[1]
bi_g_gm_conf1<-bi_g_gm$conf.int[1]
bi_g_nawm_conf1<-bi_g_nawm$conf.int[1]
bi_g_wmh_conf1<-bi_g_wmh$conf.int[1]
bi_g_fa_conf1<-bi_g_fa$conf.int[1]
bi_g_md_conf1<-bi_g_md$conf.int[1]
g_tbv_conf2<-g_tbv$conf.int[2]
g_gm_conf2<-g_gm$conf.int[2]
g_nawm_conf2<-g_nawm$conf.int[2]
g_wmh_conf2<-g_wmh$conf.int[2]
g_fa_conf2<-g_fa$conf.int[2]
g_md_conf2<-g_md$conf.int[2]
bi_g_tbv_conf2<-bi_g_tbv$conf.int[2]
bi_g_gm_conf2<-bi_g_gm$conf.int[2]
bi_g_nawm_conf2<-bi_g_nawm$conf.int[2]
bi_g_wmh_conf2<-bi_g_wmh$conf.int[2]
bi_g_fa_conf2<-bi_g_fa$conf.int[2]
bi_g_md_conf2<-bi_g_md$conf.int[2]

g_plot<-data.frame(rbind(g_tbv_r,bi_g_tbv_r,
                         g_gm_r,bi_g_gm_r,
                         g_nawm_r,bi_g_nawm_r,
                         g_wmh_r,bi_g_wmh_r,
                         g_fa_r,bi_g_fa_r,
                         g_md_r,bi_g_md_r))
g_plot$p.value<- rbind(g_tbv_p,bi_g_tbv_p,
                         g_gm_p,bi_g_gm_p,
                         g_nawm_p,bi_g_nawm_p,
                         g_wmh_p,bi_g_wmh_p,
                         g_fa_p,bi_g_fa_p,
                         g_md_p,bi_g_md_p)
g_plot$measure <- c('tbv','tbv','gm','gm','nawm','nawm',
                    'wmh','wmh','fa','fa','md','md')
g_plot$confidence1 <- rbind(g_tbv_conf1, bi_g_tbv_conf1,
                           g_gm_conf1, bi_g_gm_conf1,
                           g_nawm_conf1, bi_g_nawm_conf1,
                           g_wmh_conf1, bi_g_wmh_conf1,
                           g_fa_conf1, bi_g_fa_conf1,
                           g_md_conf1, bi_g_md_conf1)
g_plot$confidence2 <- rbind(g_tbv_conf2, bi_g_tbv_conf2,
                            g_gm_conf2, bi_g_gm_conf2,
                            g_nawm_conf2, bi_g_nawm_conf2,
                            g_wmh_conf2, bi_g_wmh_conf2,
                            g_fa_conf2, bi_g_fa_conf2,
                            g_md_conf2, bi_g_md_conf2)
g_plot$model <- c('uni','bi')
g_plot <- mutate(g_plot, psig = case_when(p.value <= 0.001 ~ "***",
                                          p.value <= 0.01 ~ "**",
                                          p.value <= 0.05 ~ "*",
                                          p.value > 0.05 ~ ""))
g_plot
#plot of g estimates and SEs
fig2a<-
  ggplot(g_plot, mapping=aes(x=measure,y=cor,
                           fill=factor(model, levels=c('uni','bi'))))+
  geom_bar(stat="identity",position = "dodge",width = 1.5)+
  geom_errorbar(aes(ymin=(confidence1[,1]), ymax=(confidence2[,1])), 
                width=.3, position=position_dodge(1.5) )+
  scale_x_discrete(limits=g_plot$measure, labels=c('TBV','','GM','',
                                                   'NAWM','','WMH','',
                                                   'gFA','','gMD',''),
                   expand = expansion(add=c(1,0))) + 
  scale_y_continuous(limits = c(-0.3,0.35), breaks=c(-0.2,-0.1,0,0.1,0.2,0.3))+
  scale_fill_manual(values = c("#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3"), 
                    limits=g_plot$model,
                    labels=c('Single-order','Bifactor'),
                    guide = guide_legend(reverse = F))+
  labs(y='Effect size', x='', fill = 'Model')+
  theme_light() + 
  theme(plot.margin=margin(t=1,r=0.25,l=0.25,unit='cm'))+
  geom_text(label=g_plot$psig, size=4, 
            position=position_dodge(1.5),
            vjust= ifelse(g_plot$cor >= 0, 
                          yes=g_plot$cor-2.5*sign(g_plot$cor), 
                          no=g_plot$cor-4*sign(g_plot$cor)))
fig2a
# Figure 2B: Model comparison: Processing speed effect size ----------------------------------------

#speed bind coeff tables

speed_tbv_r<-speed_tbv$estimate
speed_gm_r<-speed_gm$estimate
speed_nawm_r<-speed_nawm$estimate
speed_wmh_r<-speed_wmh$estimate
speed_fa_r<-speed_fa$estimate
speed_md_r<-speed_md$estimate
bi_speed_tbv_r<-bi_speed_tbv$estimate
bi_speed_gm_r<-bi_speed_gm$estimate
bi_speed_nawm_r<-bi_speed_nawm$estimate
bi_speed_wmh_r<-bi_speed_wmh$estimate
bi_speed_fa_r<-bi_speed_fa$estimate
bi_speed_md_r<-bi_speed_md$estimate
speed_tbv_p<-speed_tbv$p.value
speed_gm_p<-speed_gm$p.value
speed_nawm_p<-speed_nawm$p.value
speed_wmh_p<-speed_wmh$p.value
speed_fa_p<-speed_fa$p.value
speed_md_p<-speed_md$p.value
bi_speed_tbv_p<-bi_speed_tbv$p.value
bi_speed_gm_p<-bi_speed_gm$p.value
bi_speed_nawm_p<-bi_speed_nawm$p.value
bi_speed_wmh_p<-bi_speed_wmh$p.value
bi_speed_fa_p<-bi_speed_fa$p.value
bi_speed_md_p<-bi_speed_md$p.value
speed_tbv_conf1<-speed_tbv$conf.int[1]
speed_gm_conf1<-speed_gm$conf.int[1]
speed_nawm_conf1<-speed_nawm$conf.int[1]
speed_wmh_conf1<-speed_wmh$conf.int[1]
speed_fa_conf1<-speed_fa$conf.int[1]
speed_md_conf1<-speed_md$conf.int[1]
bi_speed_tbv_conf1<-bi_speed_tbv$conf.int[1]
bi_speed_gm_conf1<-bi_speed_gm$conf.int[1]
bi_speed_nawm_conf1<-bi_speed_nawm$conf.int[1]
bi_speed_wmh_conf1<-bi_speed_wmh$conf.int[1]
bi_speed_fa_conf1<-bi_speed_fa$conf.int[1]
bi_speed_md_conf1<-bi_speed_md$conf.int[1]
speed_tbv_conf2<-speed_tbv$conf.int[2]
speed_gm_conf2<-speed_gm$conf.int[2]
speed_nawm_conf2<-speed_nawm$conf.int[2]
speed_wmh_conf2<-speed_wmh$conf.int[2]
speed_fa_conf2<-speed_fa$conf.int[2]
speed_md_conf2<-speed_md$conf.int[2]
bi_speed_tbv_conf2<-bi_speed_tbv$conf.int[2]
bi_speed_gm_conf2<-bi_speed_gm$conf.int[2]
bi_speed_nawm_conf2<-bi_speed_nawm$conf.int[2]
bi_speed_wmh_conf2<-bi_speed_wmh$conf.int[2]
bi_speed_fa_conf2<-bi_speed_fa$conf.int[2]
bi_speed_md_conf2<-bi_speed_md$conf.int[2]

speed_plot<-data.frame(rbind(speed_tbv_r,bi_speed_tbv_r,
                         speed_gm_r,bi_speed_gm_r,
                         speed_nawm_r,bi_speed_nawm_r,
                         speed_wmh_r,bi_speed_wmh_r,
                         speed_fa_r,bi_speed_fa_r,
                         speed_md_r,bi_speed_md_r))
speed_plot$p.value<- rbind(speed_tbv_p,bi_speed_tbv_p,
                       speed_gm_p,bi_speed_gm_p,
                       speed_nawm_p,bi_speed_nawm_p,
                       speed_wmh_p,bi_speed_wmh_p,
                       speed_fa_p,bi_speed_fa_p,
                       speed_md_p,bi_speed_md_p)
speed_plot$measure <- c('tbv','tbv','gm','gm','nawm','nawm',
                    'wmh','wmh','fa','fa','md','md')
speed_plot$confidence1 <- rbind(speed_tbv_conf1, bi_speed_tbv_conf1,
                            speed_gm_conf1, bi_speed_gm_conf1,
                            speed_nawm_conf1, bi_speed_nawm_conf1,
                            speed_wmh_conf1, bi_speed_wmh_conf1,
                            speed_fa_conf1, bi_speed_fa_conf1,
                            speed_md_conf1, bi_speed_md_conf1)
speed_plot$confidence2 <- rbind(speed_tbv_conf2, bi_speed_tbv_conf2,
                            speed_gm_conf2, bi_speed_gm_conf2,
                            speed_nawm_conf2, bi_speed_nawm_conf2,
                            speed_wmh_conf2, bi_speed_wmh_conf2,
                            speed_fa_conf2, bi_speed_fa_conf2,
                            speed_md_conf2, bi_speed_md_conf2)
speed_plot$model <- c('uni','bi')
speed_plot <- mutate(speed_plot, psig = case_when(p.value <= 0.001 ~ "***",
                                          p.value <= 0.01 ~ "**",
                                          p.value <= 0.05 ~ "*",
                                          p.value > 0.05 ~ ""))
speed_plot
#plot of speed estimates and SEs
fig2b<-
  ggplot(speed_plot, mapping=aes(x=measure,y=cor,
                             fill=factor(model, levels=c('uni','bi'))))+
  geom_bar(stat="identity",position = "dodge",width = 1.5)+
  geom_errorbar(aes(ymin=(confidence1[,1]), ymax=(confidence2[,1])), 
                width=.3, position=position_dodge(1.5) )+
  scale_x_discrete(limits=speed_plot$measure, labels=c('TBV','','GM','',
                                                   'NAWM','','WMH','',
                                                   'gFA','','gMD',''),
                   expand = expansion(add=c(1,0))) + 
  scale_y_continuous(limits = c(-0.3,0.35), breaks=c(-0.2,-0.1,0,0.1,0.2,0.3))+
  scale_fill_manual(values = c("#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3"), 
                    limits=speed_plot$model,
                    labels=c('Single-order','Bifactor'),
                    guide = guide_legend(reverse = F))+
  labs(y='Effect size', x='', fill = 'Model')+
  theme_light() + 
  theme(plot.margin=margin(t=1,r=0.25,l=0.25,unit='cm'))+
  geom_text(label=speed_plot$psig, size=4, 
            position=position_dodge(1.5),
            vjust= ifelse(speed_plot$cor >= 0, 
                          yes=speed_plot$cor-2.5*sign(speed_plot$cor), 
                          no=speed_plot$cor-4*sign(speed_plot$cor)))
fig2b

# Figure 2C: Model comparison: Crystallised ability effect size ------------

#cryst bind coeff tables

cryst_tbv_r<-cryst_tbv$estimate
cryst_gm_r<-cryst_gm$estimate
cryst_nawm_r<-cryst_nawm$estimate
cryst_wmh_r<-cryst_wmh$estimate
cryst_fa_r<-cryst_fa$estimate
cryst_md_r<-cryst_md$estimate
bi_cryst_tbv_r<-bi_cryst_tbv$estimate
bi_cryst_gm_r<-bi_cryst_gm$estimate
bi_cryst_nawm_r<-bi_cryst_nawm$estimate
bi_cryst_wmh_r<-bi_cryst_wmh$estimate
bi_cryst_fa_r<-bi_cryst_fa$estimate
bi_cryst_md_r<-bi_cryst_md$estimate
cryst_tbv_p<-cryst_tbv$p.value
cryst_gm_p<-cryst_gm$p.value
cryst_nawm_p<-cryst_nawm$p.value
cryst_wmh_p<-cryst_wmh$p.value
cryst_fa_p<-cryst_fa$p.value
cryst_md_p<-cryst_md$p.value
bi_cryst_tbv_p<-bi_cryst_tbv$p.value
bi_cryst_gm_p<-bi_cryst_gm$p.value
bi_cryst_nawm_p<-bi_cryst_nawm$p.value
bi_cryst_wmh_p<-bi_cryst_wmh$p.value
bi_cryst_fa_p<-bi_cryst_fa$p.value
bi_cryst_md_p<-bi_cryst_md$p.value
cryst_tbv_conf1<-cryst_tbv$conf.int[1]
cryst_gm_conf1<-cryst_gm$conf.int[1]
cryst_nawm_conf1<-cryst_nawm$conf.int[1]
cryst_wmh_conf1<-cryst_wmh$conf.int[1]
cryst_fa_conf1<-cryst_fa$conf.int[1]
cryst_md_conf1<-cryst_md$conf.int[1]
bi_cryst_tbv_conf1<-bi_cryst_tbv$conf.int[1]
bi_cryst_gm_conf1<-bi_cryst_gm$conf.int[1]
bi_cryst_nawm_conf1<-bi_cryst_nawm$conf.int[1]
bi_cryst_wmh_conf1<-bi_cryst_wmh$conf.int[1]
bi_cryst_fa_conf1<-bi_cryst_fa$conf.int[1]
bi_cryst_md_conf1<-bi_cryst_md$conf.int[1]
cryst_tbv_conf2<-cryst_tbv$conf.int[2]
cryst_gm_conf2<-cryst_gm$conf.int[2]
cryst_nawm_conf2<-cryst_nawm$conf.int[2]
cryst_wmh_conf2<-cryst_wmh$conf.int[2]
cryst_fa_conf2<-cryst_fa$conf.int[2]
cryst_md_conf2<-cryst_md$conf.int[2]
bi_cryst_tbv_conf2<-bi_cryst_tbv$conf.int[2]
bi_cryst_gm_conf2<-bi_cryst_gm$conf.int[2]
bi_cryst_nawm_conf2<-bi_cryst_nawm$conf.int[2]
bi_cryst_wmh_conf2<-bi_cryst_wmh$conf.int[2]
bi_cryst_fa_conf2<-bi_cryst_fa$conf.int[2]
bi_cryst_md_conf2<-bi_cryst_md$conf.int[2]

cryst_plot<-data.frame(rbind(cryst_tbv_r,bi_cryst_tbv_r,
                             cryst_gm_r,bi_cryst_gm_r,
                             cryst_nawm_r,bi_cryst_nawm_r,
                             cryst_wmh_r,bi_cryst_wmh_r,
                             cryst_fa_r,bi_cryst_fa_r,
                             cryst_md_r,bi_cryst_md_r))
cryst_plot$p.value<- rbind(cryst_tbv_p,bi_cryst_tbv_p,
                           cryst_gm_p,bi_cryst_gm_p,
                           cryst_nawm_p,bi_cryst_nawm_p,
                           cryst_wmh_p,bi_cryst_wmh_p,
                           cryst_fa_p,bi_cryst_fa_p,
                           cryst_md_p,bi_cryst_md_p)
cryst_plot$measure <- c('tbv','tbv','gm','gm','nawm','nawm',
                        'wmh','wmh','fa','fa','md','md')
cryst_plot$confidence1 <- rbind(cryst_tbv_conf1, bi_cryst_tbv_conf1,
                                cryst_gm_conf1, bi_cryst_gm_conf1,
                                cryst_nawm_conf1, bi_cryst_nawm_conf1,
                                cryst_wmh_conf1, bi_cryst_wmh_conf1,
                                cryst_fa_conf1, bi_cryst_fa_conf1,
                                cryst_md_conf1, bi_cryst_md_conf1)
cryst_plot$confidence2 <- rbind(cryst_tbv_conf2, bi_cryst_tbv_conf2,
                                cryst_gm_conf2, bi_cryst_gm_conf2,
                                cryst_nawm_conf2, bi_cryst_nawm_conf2,
                                cryst_wmh_conf2, bi_cryst_wmh_conf2,
                                cryst_fa_conf2, bi_cryst_fa_conf2,
                                cryst_md_conf2, bi_cryst_md_conf2)
cryst_plot$model <- c('uni','bi')
cryst_plot <- mutate(cryst_plot, psig = case_when(p.value <= 0.001 ~ "***",
                                                  p.value <= 0.01 ~ "**",
                                                  p.value <= 0.05 ~ "*",
                                                  p.value > 0.05 ~ ""))
cryst_plot
#plot of cryst estimates and SEs
fig2c<-
  ggplot(cryst_plot, mapping=aes(x=measure,y=cor,
                                 fill=factor(model, levels=c('uni','bi'))))+
  geom_bar(stat="identity",position = "dodge",width = 1.5)+
  geom_errorbar(aes(ymin=(confidence1[,1]), ymax=(confidence2[,1])), 
                width=.3, position=position_dodge(1.5) )+
  scale_x_discrete(limits=cryst_plot$measure, labels=c('TBV','','GM','',
                                                       'NAWM','','WMH','',
                                                       'gFA','','gMD',''),
                   expand = expansion(add=c(1,0))) + 
  scale_y_continuous(limits = c(-0.3,0.35), breaks=c(-0.2,-0.1,0,0.1,0.2,0.3))+
  scale_fill_manual(values = c("#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3"), 
                    limits=cryst_plot$model,
                    labels=c('Single-order','Bifactor'),
                    guide = guide_legend(reverse = F))+
  labs(y='Effect size', x='', fill = 'Model')+
  theme_light() + 
  theme(plot.margin=margin(t=1,r=0.25,l=0.25,unit='cm'))+
  geom_text(label=cryst_plot$psig, size=4, 
            position=position_dodge(1.5),
            vjust= ifelse(cryst_plot$cor >= 0, 
                          yes=cryst_plot$cor-2.5*sign(cryst_plot$cor), 
                          no=cryst_plot$cor-4*sign(cryst_plot$cor)))
fig2c

# Figure 2D: Model comparison: Memory effect size -------------------------

#memory bind coeff tables

memory_tbv_r<-memory_tbv$estimate
memory_gm_r<-memory_gm$estimate
memory_nawm_r<-memory_nawm$estimate
memory_wmh_r<-memory_wmh$estimate
memory_fa_r<-memory_fa$estimate
memory_md_r<-memory_md$estimate
bi_memory_tbv_r<-bi_memory_tbv$estimate
bi_memory_gm_r<-bi_memory_gm$estimate
bi_memory_nawm_r<-bi_memory_nawm$estimate
bi_memory_wmh_r<-bi_memory_wmh$estimate
bi_memory_fa_r<-bi_memory_fa$estimate
bi_memory_md_r<-bi_memory_md$estimate
memory_tbv_p<-memory_tbv$p.value
memory_gm_p<-memory_gm$p.value
memory_nawm_p<-memory_nawm$p.value
memory_wmh_p<-memory_wmh$p.value
memory_fa_p<-memory_fa$p.value
memory_md_p<-memory_md$p.value
bi_memory_tbv_p<-bi_memory_tbv$p.value
bi_memory_gm_p<-bi_memory_gm$p.value
bi_memory_nawm_p<-bi_memory_nawm$p.value
bi_memory_wmh_p<-bi_memory_wmh$p.value
bi_memory_fa_p<-bi_memory_fa$p.value
bi_memory_md_p<-bi_memory_md$p.value
memory_tbv_conf1<-memory_tbv$conf.int[1]
memory_gm_conf1<-memory_gm$conf.int[1]
memory_nawm_conf1<-memory_nawm$conf.int[1]
memory_wmh_conf1<-memory_wmh$conf.int[1]
memory_fa_conf1<-memory_fa$conf.int[1]
memory_md_conf1<-memory_md$conf.int[1]
bi_memory_tbv_conf1<-bi_memory_tbv$conf.int[1]
bi_memory_gm_conf1<-bi_memory_gm$conf.int[1]
bi_memory_nawm_conf1<-bi_memory_nawm$conf.int[1]
bi_memory_wmh_conf1<-bi_memory_wmh$conf.int[1]
bi_memory_fa_conf1<-bi_memory_fa$conf.int[1]
bi_memory_md_conf1<-bi_memory_md$conf.int[1]
memory_tbv_conf2<-memory_tbv$conf.int[2]
memory_gm_conf2<-memory_gm$conf.int[2]
memory_nawm_conf2<-memory_nawm$conf.int[2]
memory_wmh_conf2<-memory_wmh$conf.int[2]
memory_fa_conf2<-memory_fa$conf.int[2]
memory_md_conf2<-memory_md$conf.int[2]
bi_memory_tbv_conf2<-bi_memory_tbv$conf.int[2]
bi_memory_gm_conf2<-bi_memory_gm$conf.int[2]
bi_memory_nawm_conf2<-bi_memory_nawm$conf.int[2]
bi_memory_wmh_conf2<-bi_memory_wmh$conf.int[2]
bi_memory_fa_conf2<-bi_memory_fa$conf.int[2]
bi_memory_md_conf2<-bi_memory_md$conf.int[2]

memory_plot<-data.frame(rbind(memory_tbv_r,bi_memory_tbv_r,
                             memory_gm_r,bi_memory_gm_r,
                             memory_nawm_r,bi_memory_nawm_r,
                             memory_wmh_r,bi_memory_wmh_r,
                             memory_fa_r,bi_memory_fa_r,
                             memory_md_r,bi_memory_md_r))
memory_plot$p.value<- rbind(memory_tbv_p,bi_memory_tbv_p,
                           memory_gm_p,bi_memory_gm_p,
                           memory_nawm_p,bi_memory_nawm_p,
                           memory_wmh_p,bi_memory_wmh_p,
                           memory_fa_p,bi_memory_fa_p,
                           memory_md_p,bi_memory_md_p)
memory_plot$measure <- c('tbv','tbv','gm','gm','nawm','nawm',
                        'wmh','wmh','fa','fa','md','md')
memory_plot$confidence1 <- rbind(memory_tbv_conf1, bi_memory_tbv_conf1,
                                memory_gm_conf1, bi_memory_gm_conf1,
                                memory_nawm_conf1, bi_memory_nawm_conf1,
                                memory_wmh_conf1, bi_memory_wmh_conf1,
                                memory_fa_conf1, bi_memory_fa_conf1,
                                memory_md_conf1, bi_memory_md_conf1)
memory_plot$confidence2 <- rbind(memory_tbv_conf2, bi_memory_tbv_conf2,
                                memory_gm_conf2, bi_memory_gm_conf2,
                                memory_nawm_conf2, bi_memory_nawm_conf2,
                                memory_wmh_conf2, bi_memory_wmh_conf2,
                                memory_fa_conf2, bi_memory_fa_conf2,
                                memory_md_conf2, bi_memory_md_conf2)
memory_plot$model <- c('uni','bi')
memory_plot <- mutate(memory_plot, psig = case_when(p.value <= 0.001 ~ "***",
                                                  p.value <= 0.01 ~ "**",
                                                  p.value <= 0.05 ~ "*",
                                                  p.value > 0.05 ~ ""))
memory_plot
#plot of memory estimates and SEs
fig2d<-
  ggplot(memory_plot, mapping=aes(x=measure,y=cor,
                                 fill=factor(model, levels=c('uni','bi'))))+
  geom_bar(stat="identity",position = "dodge",width = 1.5)+
  geom_errorbar(aes(ymin=(confidence1[,1]), ymax=(confidence2[,1])), 
                width=.3, position=position_dodge(1.5) )+
  scale_x_discrete(limits=memory_plot$measure, labels=c('TBV','','GM','',
                                                       'NAWM','','WMH','',
                                                       'gFA','','gMD',''),
                   expand = expansion(add=c(1,0))) + 
  scale_y_continuous(limits = c(-0.3,0.35), breaks=c(-0.2,-0.1,0,0.1,0.2,0.3))+
  scale_fill_manual(values = c("#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3"), 
                    limits=memory_plot$model,
                    labels=c('Single-order','Bifactor'),
                    guide = guide_legend(reverse = F))+
  labs(y='Effect size', x='', fill = 'Model')+
  theme_light() + 
  theme(plot.margin=margin(t=1,r=0.25,l=0.25,unit='cm'))+
  geom_text(label=memory_plot$psig, size=4, 
            position=position_dodge(1.5),
            vjust= ifelse(memory_plot$cor >= 0, 
                          yes=memory_plot$cor-2.5*sign(memory_plot$cor), 
                          no=memory_plot$cor-4*sign(memory_plot$cor)))
fig2d
# Figure 2E: Model comparison: Visuospatial ability effect size -----------

#visuo bind coeff tables

visuo_tbv_r<-visuo_tbv$estimate
visuo_gm_r<-visuo_gm$estimate
visuo_nawm_r<-visuo_nawm$estimate
visuo_wmh_r<-visuo_wmh$estimate
visuo_fa_r<-visuo_fa$estimate
visuo_md_r<-visuo_md$estimate
bi_visuo_tbv_r<-bi_visuo_tbv$estimate
bi_visuo_gm_r<-bi_visuo_gm$estimate
bi_visuo_nawm_r<-bi_visuo_nawm$estimate
bi_visuo_wmh_r<-bi_visuo_wmh$estimate
bi_visuo_fa_r<-bi_visuo_fa$estimate
bi_visuo_md_r<-bi_visuo_md$estimate
visuo_tbv_p<-visuo_tbv$p.value
visuo_gm_p<-visuo_gm$p.value
visuo_nawm_p<-visuo_nawm$p.value
visuo_wmh_p<-visuo_wmh$p.value
visuo_fa_p<-visuo_fa$p.value
visuo_md_p<-visuo_md$p.value
bi_visuo_tbv_p<-bi_visuo_tbv$p.value
bi_visuo_gm_p<-bi_visuo_gm$p.value
bi_visuo_nawm_p<-bi_visuo_nawm$p.value
bi_visuo_wmh_p<-bi_visuo_wmh$p.value
bi_visuo_fa_p<-bi_visuo_fa$p.value
bi_visuo_md_p<-bi_visuo_md$p.value
visuo_tbv_conf1<-visuo_tbv$conf.int[1]
visuo_gm_conf1<-visuo_gm$conf.int[1]
visuo_nawm_conf1<-visuo_nawm$conf.int[1]
visuo_wmh_conf1<-visuo_wmh$conf.int[1]
visuo_fa_conf1<-visuo_fa$conf.int[1]
visuo_md_conf1<-visuo_md$conf.int[1]
bi_visuo_tbv_conf1<-bi_visuo_tbv$conf.int[1]
bi_visuo_gm_conf1<-bi_visuo_gm$conf.int[1]
bi_visuo_nawm_conf1<-bi_visuo_nawm$conf.int[1]
bi_visuo_wmh_conf1<-bi_visuo_wmh$conf.int[1]
bi_visuo_fa_conf1<-bi_visuo_fa$conf.int[1]
bi_visuo_md_conf1<-bi_visuo_md$conf.int[1]
visuo_tbv_conf2<-visuo_tbv$conf.int[2]
visuo_gm_conf2<-visuo_gm$conf.int[2]
visuo_nawm_conf2<-visuo_nawm$conf.int[2]
visuo_wmh_conf2<-visuo_wmh$conf.int[2]
visuo_fa_conf2<-visuo_fa$conf.int[2]
visuo_md_conf2<-visuo_md$conf.int[2]
bi_visuo_tbv_conf2<-bi_visuo_tbv$conf.int[2]
bi_visuo_gm_conf2<-bi_visuo_gm$conf.int[2]
bi_visuo_nawm_conf2<-bi_visuo_nawm$conf.int[2]
bi_visuo_wmh_conf2<-bi_visuo_wmh$conf.int[2]
bi_visuo_fa_conf2<-bi_visuo_fa$conf.int[2]
bi_visuo_md_conf2<-bi_visuo_md$conf.int[2]

visuo_plot<-data.frame(rbind(visuo_tbv_r,bi_visuo_tbv_r,
                             visuo_gm_r,bi_visuo_gm_r,
                             visuo_nawm_r,bi_visuo_nawm_r,
                             visuo_wmh_r,bi_visuo_wmh_r,
                             visuo_fa_r,bi_visuo_fa_r,
                             visuo_md_r,bi_visuo_md_r))
visuo_plot$p.value<- rbind(visuo_tbv_p,bi_visuo_tbv_p,
                           visuo_gm_p,bi_visuo_gm_p,
                           visuo_nawm_p,bi_visuo_nawm_p,
                           visuo_wmh_p,bi_visuo_wmh_p,
                           visuo_fa_p,bi_visuo_fa_p,
                           visuo_md_p,bi_visuo_md_p)
visuo_plot$measure <- c('tbv','tbv','gm','gm','nawm','nawm',
                        'wmh','wmh','fa','fa','md','md')
visuo_plot$confidence1 <- rbind(visuo_tbv_conf1, bi_visuo_tbv_conf1,
                                visuo_gm_conf1, bi_visuo_gm_conf1,
                                visuo_nawm_conf1, bi_visuo_nawm_conf1,
                                visuo_wmh_conf1, bi_visuo_wmh_conf1,
                                visuo_fa_conf1, bi_visuo_fa_conf1,
                                visuo_md_conf1, bi_visuo_md_conf1)
visuo_plot$confidence2 <- rbind(visuo_tbv_conf2, bi_visuo_tbv_conf2,
                                visuo_gm_conf2, bi_visuo_gm_conf2,
                                visuo_nawm_conf2, bi_visuo_nawm_conf2,
                                visuo_wmh_conf2, bi_visuo_wmh_conf2,
                                visuo_fa_conf2, bi_visuo_fa_conf2,
                                visuo_md_conf2, bi_visuo_md_conf2)
visuo_plot$model <- c('uni','bi')
visuo_plot <- mutate(visuo_plot, psig = case_when(p.value <= 0.001 ~ "***",
                                                  p.value <= 0.01 ~ "**",
                                                  p.value <= 0.05 ~ "*",
                                                  p.value > 0.05 ~ ""))
visuo_plot
#plot of visuo estimates and SEs
fig2e<-
  ggplot(visuo_plot, mapping=aes(x=measure,y=cor,
                                 fill=factor(model, levels=c('uni','bi'))))+
  geom_bar(stat="identity",position = "dodge",width = 1.5)+
  geom_errorbar(aes(ymin=(confidence1[,1]), ymax=(confidence2[,1])), 
                width=.3, position=position_dodge(1.5) )+
  scale_x_discrete(limits=visuo_plot$measure, labels=c('TBV','','GM','',
                                                       'NAWM','','WMH','',
                                                       'gFA','','gMD',''),
                   expand = expansion(add=c(1,0))) + 
  scale_y_continuous(limits = c(-0.3,0.35), breaks=c(-0.2,-0.1,0,0.1,0.2,0.3))+
  scale_fill_manual(values = c("#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3","#B2182B","#4393C3"), 
                    limits=visuo_plot$model,
                    labels=c('Single-order','Bifactor'),
                    guide = guide_legend(reverse = F))+
  labs(y='Effect size', x='', fill = 'Model')+
  theme_light() + 
  theme(plot.margin=margin(t=1,r=0.25,l=0.25,unit='cm'))+
  geom_text(label=visuo_plot$psig, size=4, 
            position=position_dodge(1.5),
            vjust= ifelse(visuo_plot$cor >= 0, 
                          yes=visuo_plot$cor-2.5*sign(visuo_plot$cor), 
                          no=visuo_plot$cor-4*sign(visuo_plot$cor)))
fig2e

# Figure 2: Model comparison of effect sizes ---------------------------------------------------------------
ggsave('fig2.png', plot=(
ggarrange(fig2b,fig2c,fig2d,fig2e,fig2a,
          ncol=3, nrow=2,
          labels = c('Processing speed','Crystallised ability',
                     'Memory','Visuospatial ability','g'),
          label.x = c(.15,.15,.35,.15,.5),
          label.y = c(1,1,1,1,1),
          #align='hv',
          font.label = list(size = 12, color = "black"),
          common.legend = TRUE, legend = 'bottom',
          theme(plot.margin=margin(t=1,r=0.25,l=0.25,unit='cm')))),
width= 30, units='cm', height = 20)


# % attenuation speed
100-( ((bi_speed_tbv$estimate/speed_tbv$estimate)*100) + ((bi_speed_gm$estimate/speed_gm$estimate)*100) + ((bi_speed_nawm$estimate/speed_nawm$estimate)*100) +
((bi_speed_wmh$estimate/speed_wmh$estimate)*100) + ((bi_speed_fa$estimate/speed_fa$estimate)*100) + ((bi_speed_md$estimate/speed_md$estimate)*100) ) / 6
(100* ((speed_tbv$estimate - bi_speed_tbv$estimate)/speed_tbv$estimate))
(100* ((speed_gm$estimate - bi_speed_gm$estimate)/speed_gm$estimate))
(100* ((speed_nawm$estimate - bi_speed_nawm$estimate)/speed_nawm$estimate))
(100* ((speed_wmh$estimate - bi_speed_wmh$estimate)/speed_wmh$estimate))
(100* ((speed_fa$estimate - bi_speed_fa$estimate)/speed_fa$estimate))
(100* ((speed_md$estimate - bi_speed_md$estimate)/speed_md$estimate))

# % attenuation cryst
100-( ((bi_cryst_tbv$estimate/cryst_tbv$estimate)*100) + ((bi_cryst_gm$estimate/cryst_gm$estimate)*100) + ((bi_cryst_nawm$estimate/cryst_nawm$estimate)*100) +
    ((bi_cryst_wmh$estimate/cryst_wmh$estimate)*100) + ((bi_cryst_fa$estimate/cryst_fa$estimate)*100) + ((bi_cryst_md$estimate/cryst_md$estimate)*100) ) / 6
a<- (100* ((cryst_tbv$estimate - -bi_cryst_tbv$estimate)/cryst_tbv$estimate))
b<- (100* ((cryst_gm$estimate - -bi_cryst_gm$estimate)/cryst_gm$estimate))
c<- (100* ((cryst_nawm$estimate - -bi_cryst_nawm$estimate)/cryst_nawm$estimate))
d<- (100* (((cryst_wmh$estimate) - bi_cryst_wmh$estimate)/(cryst_wmh$estimate)))
e<- (100* ((cryst_fa$estimate - -bi_cryst_fa$estimate)/cryst_fa$estimate))
f<- (100* ((cryst_md$estimate - bi_cryst_md$estimate)/cryst_md$estimate))

(a+b+c+d+e+f) /6
a
b
c
d
e
f


# % attenuation memory
100-( ((bi_memory_tbv$estimate/memory_tbv$estimate)*100) + ((bi_memory_gm$estimate/memory_gm$estimate)*100) + ((bi_memory_nawm$estimate/memory_nawm$estimate)*100) +
    ((bi_memory_wmh$estimate/memory_wmh$estimate)*100) + ((bi_memory_fa$estimate/memory_fa$estimate)*100) + ((bi_memory_md$estimate/memory_md$estimate)*100) ) / 6
100* ((memory_tbv$estimate - bi_memory_tbv$estimate)/memory_tbv$estimate)
100* ((memory_gm$estimate - bi_memory_gm$estimate)/memory_gm$estimate)
100* ((memory_nawm$estimate - bi_memory_nawm$estimate)/memory_nawm$estimate)
100* ((memory_wmh$estimate - bi_memory_wmh$estimate)/memory_wmh$estimate)
100* ((memory_fa$estimate - bi_memory_fa$estimate)/memory_fa$estimate)
100* ((memory_md$estimate - bi_memory_md$estimate)/memory_md$estimate)

# % attenuation visuo
100-( ((bi_visuo_tbv$estimate/visuo_tbv$estimate)*100) + ((bi_visuo_gm$estimate/visuo_gm$estimate)*100) + ((bi_visuo_nawm$estimate/visuo_nawm$estimate)*100) +
    ((bi_visuo_wmh$estimate/visuo_wmh$estimate)*100) + ((bi_visuo_fa$estimate/visuo_fa$estimate)*100) + ((bi_visuo_md$estimate/visuo_md$estimate)*100) ) / 6
100* ((visuo_tbv$estimate - bi_visuo_tbv$estimate)/visuo_tbv$estimate)
100* ((visuo_gm$estimate - bi_visuo_gm$estimate)/visuo_gm$estimate)
100* ((visuo_nawm$estimate - bi_visuo_nawm$estimate)/visuo_nawm$estimate)
100* ((visuo_wmh$estimate - bi_visuo_wmh$estimate)/visuo_wmh$estimate)
100* ((visuo_fa$estimate - bi_visuo_fa$estimate)/visuo_fa$estimate)
100* ((visuo_md$estimate - bi_visuo_md$estimate)/visuo_md$estimate)

# % attenuation g
100-( ((bi_g_tbv$estimate/g_tbv$estimate)*100) + ((bi_g_gm$estimate/g_gm$estimate)*100) + ((bi_g_nawm$estimate/g_nawm$estimate)*100) +
    ((bi_g_wmh$estimate/g_wmh$estimate)*100) + ((bi_g_fa$estimate/g_fa$estimate)*100) + ((bi_g_md$estimate/g_md$estimate)*100) ) / 6
100* ((g_tbv$estimate - bi_g_tbv$estimate)/g_tbv$estimate)
100* ((g_gm$estimate - bi_g_gm$estimate)/g_gm$estimate)
100* ((g_nawm$estimate - bi_g_nawm$estimate)/g_nawm$estimate)
100* ((g_wmh$estimate - bi_g_wmh$estimate)/g_wmh$estimate)
100* ((g_fa$estimate - bi_g_fa$estimate)/g_fa$estimate)
100* ((g_md$estimate - bi_g_md$estimate)/g_md$estimate)

# Create matlab dataframe -------------------------------------------------

#POST-CORRECTIONS VERSION - MATLABDF_2
matlabDF_2 = data.frame(lbc36no = data$lbc36no, sex = data$sex, agedays_w2 = data$agedays_w2,
                        speed = data$speed, cryst = data$cryst, memory = data$memory, visuo = data$visuo, g = data$g,
                        bi_speed = data$bi_speed, bi_cryst = data$bi_cryst, bi_memory = data$bi_memory, bi_visuo = data$bi_visuo, bi_g = data$bi_g,
                        icv = data$ICV_mm3_wX)
matlabDF_2 <- matlabDF_2[complete.cases(matlabDF_2),]
write.csv(matlabDF_2, file = 'matlabDF_2.csv', row.names = F)

##PRE-CORRECTIONS
## = data.frame(lbc36no = data$lbc36no, sex = data$sex, agedays_w2 = data$agedays_w2, 
                        ##speed = data$speed, cryst = data$cryst, memory = data$memory, visuo = data$visuo, g = data$g,
                        ##bi_speed = data$bi_speed, bi_cryst = data$bi_cryst, bi_memory = data$bi_memory, bi_visuo = data$bi_visuo, bi_g = data$bi_g,
                        ##icv = data$ICV_mm3_wX)
##matlabDF_1 <- matlabDF_1[complete.cases(matlabDF_1),]
##write.csv(matlabDF_1, file = 'matlabDF_1.csv', row.names = F)

#matlabDF = data.frame(lbc36no = data$lbc36no, sex = data$sex, agedays_w2 = data$agedays_w2,
#                      yrsedu_w1 = data$yrsedu_w1, hmsonum_w1 = data$hmsonum_w1, fathclass_w1 = data$fathclass_w1,
#                      alcunitwk_w2 = data$alcunitwk_w2, hibp_w2 = data$hibp_w2, diab_w2 = data$diab_w2,
#                      hichol_w2 = data$hichol_w2, cvdhist_w2 = data$cvdhist_w2, dement_w2 = data$dement_w2,
#                      smokcurr_w2 = data$smokcurr_w2, smokprev_w2 = data$smokprev_w2, 
#                      speed = data$speed, cryst = data$cryst, memory = data$memory, visuo = data$visuo, g = data$g,
#                      bi_speed = data$bi_speed, bi_cryst = data$bi_cryst, bi_memory = data$bi_memory, bi_visuo = data$bi_visuo, bi_g = data$bi_g)
#matlabDF <- matlabDF[complete.cases(matlabDF),]
#write.csv(matlabDF, file = 'matlabDF.csv', row.names = F)



# Figure 5: Comparison of overlapping vertices ----------------------------------------------------------------

#Individual modelling
totalvert <- c(152085, 41274, 39937, 183014, 177816)
vertices <- c(14068,162,642,24794,5565,138017,41112,39295,158220,172251)
type <- factor(c('unique','unique','unique','unique','unique','overlap','overlap','overlap','overlap','overlap'), levels = c('unique','overlap'))
domain <- factor(c('speed','cryst','memory','visuo','g'), levels = c('speed','cryst','memory','visuo','g'))
fig5 = data.frame(domain, vertices, type)

fig5A <- (ggplot(fig5, mapping=aes(x=domain, y=vertices, fill=type))+
            geom_bar(stat="identity",position = "stack",color='black')+ 
            ggtitle('Single-order')+
            scale_x_discrete(name = '', limits = c('speed','cryst','memory','visuo','g'), labels=c('Processing\nspeed','Crystallised\nability','Memory','Visuospatial\nability','g'))+
            scale_y_continuous(name = 'Number of significantly associated vertices', limits = c(-0,190000), breaks=c(0,20000,40000,60000,80000,100000,120000,140000,160000,180000), )+
            scale_fill_manual(values=c("#B2182B","#4393C3"), guide = guide_legend(reverse = T), name='', labels = c('Unique','Overlapping'))+
            theme_light()+
            theme(plot.title = element_text(hjust = 0.5)))+
            geom_text(aes(label = c("152085","","","","","","","","","")),
                      position = position_stack(1.3))+
            geom_text(aes(label = c("","41274","","","","","","","","")),
                      position = position_stack(28))+
            geom_text(aes(label = c("","","39937","","","","","","","")),
                      position = position_stack(7))+
            geom_text(aes(label = c("","","","183014","","","","","","")),
                      position = position_stack(1.16))+
            geom_text(aes(label = c("","","","","177816","","","","","")),
                      position = position_stack(1.7))
fig5A

#Bifactor modelling
totalvert <- c(20082,0,0,3659,177343)
vertices <- c(3175,0,0,469,19450,16907,0,0,3190,157893)
type <- c('unique','unique','unique','unique','unique','overlap','overlap','overlap','overlap','overlap')
domain <- factor(c('speed','cryst','memory','visuo','g'), levels = c('speed','cryst','memory','visuo','g'))
fig5 = data.frame(domain, vertices, type)

fig5B <- (ggplot(fig5, mapping=aes(x=domain, y=vertices, fill=type))+
            geom_bar(stat="identity",position = "stack",color='black')+ 
            ggtitle('Bifactor')+
            scale_x_discrete(name = '', limits = c('speed','cryst','memory','visuo','g'), labels=c('Processing\nspeed','Crystallised\nability','Memory','Visuospatial\nability','g'))+
            scale_y_continuous(name = '', limits = c(-0,190000), breaks=c(0,20000,40000,60000,80000,100000,120000,140000,160000,180000), )+
            scale_fill_manual(values=c("#B2182B","#4393C3"), guide = guide_legend(reverse = T), name='', labels = c('Unique','Overlapping'))+
            theme_light()+
            theme(plot.title = element_text(hjust = 0.5)))+
            geom_text(aes(label = c("20082","","","","","","","","","")),
                      position = position_stack(7.8))+
            geom_text(aes(label = c("","","","","","","0","0","","")),
                      nudge_y = 4000)+
            geom_text(aes(label = c("","","","3659","","","","","","")),
                      position = position_stack(17))+
            geom_text(aes(label = c("","","","","177343","","","","","")),
                      position = position_stack(9.35))
fig5B

ggsave('fig5.png', plot=(ggarrange(fig5A,fig5B, ncol=2, nrow=1,
                                   common.legend = TRUE, legend = 'bottom')),
       width= 40, units='cm', height = 20)

