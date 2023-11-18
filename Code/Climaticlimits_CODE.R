#Load packages
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape)
library(zoo)
library(Rmpfr)
library(gmp)


#Source IPM functions R script (one that varies at0 and z0)
source(".../Code/IPM_functions.R")

#Initial abundance
initial_abundance<- rep(0, n)
initial_abundance[1]<-1000

MT_z0<-0.00942460726973781
MT_at0<- 0.0145767684333848

at0 = c(0.00728838421, 0.01093257632, MT_at0)
z0 = c(0.01178075908, MT_z0)

#Create parameter dataset
Params2test<-crossing(at0, z0)
Params2test$Growth<- c( "50red", "50red", "25red", "25red", "MTG", "MTG")
Params2test$Survival<- c("MTS", "LowS", "MTS", "LowS", "MTS", "LowS")

################################################################################
#colours for plotting 
colorBlind8<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


#LAMBDAS FIXED TEMP
#Temp 0
LT0<-calc_fxd_temps_lambda(fxd_Temp = 0,  Parameter_data = Params2test)
#Temp2.5
LT2.5<-calc_fxd_temps_lambda(fxd_Temp = 2.5,  Parameter_data = Params2test)
#Temp 5
LT5<-calc_fxd_temps_lambda(fxd_Temp = 5,  Parameter_data = Params2test)
#Temp 7.5
LT7.5<-calc_fxd_temps_lambda(fxd_Temp = 7.5,  Parameter_data = Params2test)
#Temp 10
LT10<-calc_fxd_temps_lambda(fxd_Temp = 10,  Parameter_data = Params2test)
#Temp 12.5 
LT12.5<-calc_fxd_temps_lambda(fxd_Temp = 12.5,  Parameter_data = Params2test)
#Temp 15
LT15<-calc_fxd_temps_lambda(fxd_Temp = 15,  Parameter_data = Params2test)
#Temp 17.5
LT17.5<-calc_fxd_temps_lambda(fxd_Temp = 17.5,  Parameter_data = Params2test)
#Temp 20
LT20<-calc_fxd_temps_lambda(fxd_Temp = 20,  Parameter_data = Params2test)
#Temp 22.5
LT22.5<-calc_fxd_temps_lambda(fxd_Temp = 22.5,  Parameter_data = Params2test)
#Temp 25
LT25<-calc_fxd_temps_lambda(fxd_Temp = 25,  Parameter_data = Params2test)
#Temp 27.5
LT27.5<-calc_fxd_temps_lambda(fxd_Temp = 27.5,  Parameter_data = Params2test)
#Temp 30
LT30<-calc_fxd_temps_lambda(fxd_Temp = 30,  Parameter_data = Params2test)


fixd_temps_lambdas<-rbind(LT0, LT2.5, LT5, LT7.5, LT10, 
                          LT12.5, LT15, LT17.5 ,LT20, LT22.5,
                          LT25, LT27.5, LT30)


fixd_temps_lambdas$Temp_model<- "Mean"
fixd_temps_lambdas$Param_combo <- paste(fixd_temps_lambdas$Growth, fixd_temps_lambdas$Survival, sep="")



###########################################################################################
#MEDIANS FIXED TEMP
#Temp 0
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 0, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 0, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 0, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 0, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 0, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 0, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 0
PT0<-cbind(Params2test, peak, Temp)

#Temp 2.5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 2.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 2.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 2.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 2.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 2.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 2.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 2.5
PT2.5<-cbind(Params2test, peak, Temp)

#Temp 5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 5
PT_5<-cbind(Params2test, peak, Temp)

#Temp 7.5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 7.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 7.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 7.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 7.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 7.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 7.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 7.5
PT7.5<-cbind(Params2test, peak, Temp)

#Temp 10
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 10, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 10, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 10, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 10, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 10, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 10, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 10
PT10<-cbind(Params2test, peak, Temp)

#Temp 12.5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 12.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 12.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 12.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 12.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 12.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 12.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 12.5
PT12.5<-cbind(Params2test, peak, Temp)

#Temp 15
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 15, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 15, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 15, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 15, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 15, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 15, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 15
PT15<-cbind(Params2test, peak, Temp)

#Temp 17.5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 17.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 17.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 17.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 17.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 17.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 17.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 17.5
PT17.5<-cbind(Params2test, peak, Temp)

#Temp 20
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 20, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 20, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 20, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 20, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 20, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 20, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 20
PT20<-cbind(Params2test, peak, Temp)

#Temp 22.5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 22.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 22.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 22.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 22.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 22.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 22.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 22.5
PT22.5<-cbind(Params2test, peak, Temp)

#Temp 25
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 25, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 25, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 25, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 25, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 25, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 25, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 25
PT25<-cbind(Params2test, peak, Temp)

#Temp 27.5
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 27.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 27.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 27.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 27.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 27.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 27.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 27.5
PT27.5<-cbind(Params2test, peak, Temp)

#Temp 30
PT1<-aug_median_mass_fxdtemp(fxd_Temp = 30, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_fxdtemp(fxd_Temp = 30, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_fxdtemp(fxd_Temp = 30, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_fxdtemp(fxd_Temp = 30, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_fxdtemp(fxd_Temp = 30, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_fxdtemp(fxd_Temp = 30, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 30
PT30<-cbind(Params2test, peak, Temp)

fixd_temps_peaks<-rbind(PT0, PT2.5, PT_5, PT7.5, PT10, PT12.5, PT15, PT17.5, PT20,
                                     PT22.5, PT25, PT27.5, PT30)


fixd_temps_peaks$Param_combo <- paste(fixd_temps_peaks$Growth, fixd_temps_peaks$Survival, sep="")

################################################################################################
#FLUCTATING LAMBDAS


#Temp 0
LFT0<-calc_varied_temps_lambda(vary_Temp = 0,  Parameter_data = Params2test)
#Temp2.5
LFT2.5<-calc_varied_temps_lambda(vary_Temp = 2.5,  Parameter_data = Params2test)
#Temp 5
LFT5<-calc_varied_temps_lambda(vary_Temp = 5,  Parameter_data = Params2test) #run from here
#Temp 7.5
LFT7.5<-calc_varied_temps_lambda(vary_Temp = 7.5,  Parameter_data = Params2test)
#Temp 10
LFT10<-calc_varied_temps_lambda(vary_Temp = 10,  Parameter_data = Params2test)
#Temp 12.5 
LFT12.5<-calc_varied_temps_lambda(vary_Temp = 12.5,  Parameter_data = Params2test)
#Temp 15
LFT15<-calc_varied_temps_lambda(vary_Temp = 15,  Parameter_data = Params2test)
#Temp 17.5
LFT17.5<-calc_varied_temps_lambda(vary_Temp = 17.5,  Parameter_data = Params2test)
#Temp 20
LFT20<-calc_varied_temps_lambda(vary_Temp = 20,  Parameter_data = Params2test)
#Temp 22.5
LFT22.5<-calc_varied_temps_lambda(vary_Temp = 22.5,  Parameter_data = Params2test)
#Temp 25
LFT25<-calc_varied_temps_lambda(vary_Temp = 25,  Parameter_data = Params2test)
#Temp 27.5
LFT27.5<-calc_varied_temps_lambda(vary_Temp = 27.5,  Parameter_data = Params2test)
#Temp 30
LFT30<-calc_varied_temps_lambda(vary_Temp = 30,  Parameter_data = Params2test)


lambdas_fluctuating<-rbind(LFT0, LFT2.5, LFT5, LFT7.5, LFT10, 
                        LFT12.5, LFT15, LFT17.5 ,LFT20, LFT22.5, LFT25,
                        LFT27.5, LFT30)

lambdas_fluctuating$Temp_model<- "Fluctuating"

lambdas_fluctuating$Param_combo <- paste(lambdas_fluctuating$Growth, lambdas_fluctuating$Survival, sep="_")


########################################################################################
#Temp 0
PT1<-aug_median_mass_vary_temp(vary_Temp = 0, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 0, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 0, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 0, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 0, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 0, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 0
VPT0<-cbind(Params2test, peak, Temp)

#Temp 2.5
PT1<-aug_median_mass_vary_temp(vary_Temp = 2.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 2.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 2.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 2.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 2.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 2.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 2.5
VPT2.5<-cbind(Params2test, peak, Temp)

#Temp 5
PT1<-aug_median_mass_vary_temp(vary_Temp = 5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 5
VPT5<-cbind(Params2test, peak, Temp)

#Temp 7.5
PT1<-aug_median_mass_vary_temp(vary_Temp = 7.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 7.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 7.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 7.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 7.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 7.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 7.5
VPT7.5<-cbind(Params2test, peak, Temp)

#Temp 10
PT1<-aug_median_mass_vary_temp(vary_Temp = 10, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 10, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 10, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 10, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 10, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 10, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 10
VPT10<-cbind(Params2test, peak, Temp)

#Temp 12.5
PT1<-aug_median_mass_vary_temp(vary_Temp = 12.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 12.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 12.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 12.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 12.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 12.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 12.5
VPT12.5<-cbind(Params2test, peak, Temp)

#Temp 15
PT1<-aug_median_mass_vary_temp(vary_Temp = 15, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 15, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 15, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 15, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 15, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 15, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 15
VPT15<-cbind(Params2test, peak, Temp)

#Temp 17.5
PT1<-aug_median_mass_vary_temp(vary_Temp = 17.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 17.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 17.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 17.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 17.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 17.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 17.5
VPT17.5<-cbind(Params2test, peak, Temp)

#Temp 20
PT1<-aug_median_mass_vary_temp(vary_Temp = 20, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 20, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 20, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 20, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 20, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 20, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 20
VPT20<-cbind(Params2test, peak, Temp)

#Temp 22.5
PT1<-aug_median_mass_vary_temp(vary_Temp = 22.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 22.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 22.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 22.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 22.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 22.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 22.5
VPT22.5<-cbind(Params2test, peak, Temp)

#Temp 25
PT1<-aug_median_mass_vary_temp(vary_Temp = 25, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 25, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 25, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 25, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 25, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 25, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 25
VPT25<-cbind(Params2test, peak, Temp)

#Temp 27.5
PT1<-aug_median_mass_vary_temp(vary_Temp = 27.5, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 27.5, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 27.5, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 27.5, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 27.5, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 27.5, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 27.5
VPT27.5<-cbind(Params2test, peak, Temp)

#Temp 30
PT1<-aug_median_mass_vary_temp(vary_Temp = 30, at0 = Params2test$at0[1], z0 = Params2test$z0[1])
PT2<-aug_median_mass_vary_temp(vary_Temp = 30, at0 = Params2test$at0[2], z0 = Params2test$z0[2])
PT3<-aug_median_mass_vary_temp(vary_Temp = 30, at0 = Params2test$at0[3], z0 = Params2test$z0[3])
PT4<-aug_median_mass_vary_temp(vary_Temp = 30, at0 = Params2test$at0[4], z0 = Params2test$z0[4])
PT5<-aug_median_mass_vary_temp(vary_Temp = 30, at0 = Params2test$at0[5], z0 = Params2test$z0[5])
PT6<-aug_median_mass_vary_temp(vary_Temp = 30, at0 = Params2test$at0[6], z0 = Params2test$z0[6])

peak<-c(PT1, PT2, PT3, PT4, PT5, PT6)
Temp = 30
VPT30<-cbind(Params2test, peak, Temp)


vary_temps_peaks<-rbind(VPT0, VPT2.5, VPT5, VPT7.5, VPT10, VPT12.5, VPT15, VPT17.5, VPT20,
                        VPT22.5, VPT25, VPT27.5, VPT30)


vary_temps_peaks$Param_combo <- paste(vary_temps_peaks$Growth, vary_temps_peaks$Survival, sep="_")

############################################################################################
#rmax into into one data
both<-rbind(fixd_temps_lambdas, lambdas_fluctuating)
both$Temp_model<- factor(both$Temp_model, levels = c("Mean", "Fluctuating" ))


#medians add lambdas fixed
#fixd_temps_lambdasno30<-fixd_temps_lambdas[fixd_temps_lambdas$Temp < 27.5, ]
lambdas_fluctuatingno30<-lambdas_fluctuating[lambdas_fluctuating$Temp < 27.5, ]
#add lambdas fluct
fixd_temps_lambdas<-fixd_temps_lambdas[order(fixd_temps_lambdas$Param_combo),]
fixd_temps_peaks<-fixd_temps_peaks[order(fixd_temps_peaks$Param_combo),]

fixd_temps_peaks$lambda_1<- fixd_temps_lambdas$Lambda >= 1

lambdas_fluctuatingno30<-lambdas_fluctuatingno30[order(lambdas_fluctuatingno30$Param_combo),]
vary_temps_peaks<-vary_temps_peaks[order(vary_temps_peaks$Param_combo),]

vary_temps_peaks$lambda_1<- lambdas_fluctuatingno30$Lambda >= 1

############################################################################################
#PLOTTING paper

#Figure 2
both$Resource<- ifelse(both$Growth == "50red", "50% reduction", 
                              ifelse(both$Growth == "MTG", "Metabolically optimal", "25% reduction"
                              ))
both$Resource<- factor(both$Resource, levels = c("Metabolically optimal","25% reduction","50% reduction"))

#try plotting removing 30 stuff 
bothno30<-both[both$Temp < 27.5, ]

a<-ggplot(subset(bothno30, Survival == "MTS"),  aes(Temp, log(Lambda), group = Param_combo, colour = Resource, linetype = Temp_model))+
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = 3, xmax = 15,
           ymin = -Inf, ymax = Inf) +
  #annotate("rect", fill = "lightgrey", alpha = 0.4, 
  #         xmin = 0, xmax = 23,
  #         ymin = -Inf, ymax = Inf) +
  geom_line(stat="smooth", size =1.5, se = F)+
  theme_cowplot(font_size = 18)+
  scale_colour_manual(values = colorBlind8)+
  labs(colour = "Resource level", y = expression(r[max]), x = "Temperature", linetype = "Thermal regime feature")+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")+
  xlim(c(0,25))
  #ylim(c(-25,65))


ggplot(subset(bothno30, Survival == "MTS"),  aes(Temp, log(Lambda), group = Param_combo, colour = Resource, linetype = Temp_model))+
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = 3, xmax = 15,
           ymin = -Inf, ymax = Inf) +
  #annotate("rect", fill = "lightgrey", alpha = 0.4, 
  #         xmin = 0, xmax = 23,
  #         ymin = -Inf, ymax = Inf) +
  geom_point()+
  #geom_line(data = objectp2, aes(x , y, colour = colour, group = linetype), size =1.5)+
  theme_cowplot(font_size = 18)+
  scale_colour_manual(values = colorBlind8)+
  labs(colour = "Resource level", y = expression(r[max]), x = "Temperature", linetype = "Thermal regime feature")+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")+
  xlim(c(0,25))


#median mass figure
fixd_temps_peaks$Growth<- factor(fixd_temps_peaks$Growth, levels = c("MTG","25red","50red"))

b<-ggplot(subset(fixd_temps_peaks, lambda_1 == TRUE & Survival == "MTS"), aes(Temp, peak, colour = Growth, group = Param_combo))+
  geom_line(size= 1.5)+
  #scale_shape_manual(values=c(3, 15, 17))+
  theme_cowplot(font_size = 18)+
  scale_colour_manual(values = colorBlind8)+
  #theme(legend.position = "none", plot.background = element_rect(fill = "transparent", colour = NA),
   #     panel.background = element_rect(fill = "transparent", colour = NA))+
  labs(x = "Temperature", y = "Median mass in August")+
  geom_line(data =  subset(vary_temps_peaks, lambda_1 == TRUE & Survival == "MTS"), aes(Temp, peak, colour = Growth, group = Param_combo) ,size = 1.5, linetype = "dashed")+
  xlim(c(0,25))


key<-get_legend(a)

fig2<-plot_grid(a+theme(legend.position = "none"),b+theme(legend.position = "none"), nrow = 1, labels = "auto", label_size = 18)
plot_grid(fig2, leg, rel_widths = c(1,.3))

#with new plot c
plot_grid(a+theme(legend.position = "none"),b+theme(legend.position = "none"),paprc_2, nrow = 2, labels = "auto", label_size = 18)
plot_grid(key,templeg)

#Figure 3
scen1<-ggplot(subset(bothno30, Survival == "MTS"),  aes(Temp, log(Lambda), group = Param_combo, colour = Resource, linetype = Temp_model))+
  annotate("rect", fill = "grey", 
           xmin = 3, xmax = 15,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red", alpha = 0.5, 
           xmin = 3+1.8, xmax = 15+1.8,
          ymin = -Inf, ymax = Inf) +
  geom_line(stat="smooth", size =1.5, se = F)+
  theme_cowplot(font_size = 18)+
  scale_colour_manual(values = colorBlind8)+
  labs(title = "1.5", y =  expression(r[max]), x = "Temperature", linetype = "Temperature model")+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")+
  xlim(c(0,25))+
  coord_cartesian(ylim = c(-5,65))

scen4<-ggplot(subset(bothno30, Survival == "MTS"),   aes(Temp, log(Lambda), group = Param_combo, colour = Resource, linetype = Temp_model))+
  annotate("rect", fill = "grey", 
           xmin = 3, xmax = 15,
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "red3", alpha = 0.5, 
           xmin = 3+3.6, xmax = 15+3.6,
           ymin = -Inf, ymax = Inf) +
  geom_line(stat="smooth", size =1.5, se = F)+
  theme_cowplot(font_size = 18)+
  scale_colour_manual(values = colorBlind8)+
  labs(title = "3.6", y =  expression(r[max]), x = "Temperature", linetype = "Temperature model")+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")+
  xlim(c(0,25))+
  coord_cartesian(ylim = c(-5,65))

fig3<-plot_grid(scen1 +theme(legend.position = "none"), 
          scen4 +theme(legend.position = "none"), nrow = 1, labels = "auto", label_size = 18)
plot_grid(fig3, key, rel_widths = c(1,.3,1))


#Supp fig 1

bothno30$mort<- ifelse(bothno30$Survival == "LowS", "25% increase",  "Intrinsic mortality")
bothno30$mort<- factor(bothno30$mort, levels = c("Intrinsic mortality","25% increase"))

ggplot(bothno30,  aes(Temp, log(Lambda), group = Param_combo, colour = Resource, linetype = Temp_model))+
  annotate("rect", fill = "grey", alpha = 0.7, 
           xmin = 3, xmax = 15,
           ymin = -Inf, ymax = Inf) +
  #annotate("rect", fill = "lightgrey", alpha = 0.4, 
  #         xmin = 0, xmax = 23,
   #        ymin = -Inf, ymax = Inf) +
  geom_line(stat="smooth", size =1.5, se = F)+
  theme_cowplot(font_size = 18)+
  scale_colour_manual(values = colorBlind8)+
  labs(colour = "Resource level", y = expression(r[max]), x = "Temperature", linetype = "Thermal regime feature")+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")+
  xlim(c(0,25))+
  facet_wrap(~mort)+
  coord_cartesian(ylim = c(-5,65))






#r max and ssd


rmax<-subset(bothno30, Survival == "MTS")

ssds<-subset(fixd_temps_peaks, lambda_1 == TRUE & Survival == "MTS")

ssds_v<-subset(vary_temps_peaks, lambda_1 == TRUE & Survival == "MTS")

ssds$Temp_model<- "Mean"
ssds_v$Temp_model<- "Fluctuating"

ssds_both<-rbind(ssds, ssds_v)

test2<-merge(rmax, ssds)
test2$Growth<- factor(test2$Growth, levels = c("MTG","25red","50red"))



paprc<-ggplot(test2,  aes(y = log(Lambda), x = peak, group = Growth, colour = Temp))+
  #geom_line(stat="smooth", method = "loess", size = 1.5)+
  geom_line(size =2.5, lineend = "round" )+
  #geom_smoot( method = "loess", se = F)+ #aes(colour = ..Temp..), size = 1.5,)+
  #geom_point(size =3.5)+
  theme_cowplot(font_size = 18)+
  scale_colour_continuous(low = "skyblue", high = "tomato")+
  labs(y = expression(r[max]), x = "Median mass in August")+
  ylim(c(0,65))+
  #facet_wrap(~Temp_model)+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")



ggplot(test2,  aes(y = log(Lambda), x = peak, group = Growth, colour = Temp))+
  #geom_line(stat="smooth", method = "loess", size = 1.5)+
  geom_line( size =2.5, lineend = "round" )+
  #geom_smoot( method = "loess", se = F)+ #aes(colour = ..Temp..), size = 1.5,)+
  #geom_point(size =3.5 ,aes(colour = Temp))+
  theme_cowplot(font_size = 18)+
  scale_colour_continuous(low = "skyblue", high = "tomato")+
  labs(y = expression(r[max]), x = "Median mass in August")+
  ylim(c(0,65))+
  facet_wrap(~Temp_model)+
  geom_hline(yintercept=0, linetype = "longdash", colour = "darkgrey")


paprc_2<-plot_grid(paprc+theme(legend.position = "none"), labels = "c", label_size = 18, nrow = 1)
templeg<-get_legend(paprc)
  