#Load packages needed
library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape)
library(zoo)
library(Rmpfr)
library(gmp)


#Source IPM functions R script
source("Code/MIPM_functions.R")

#Initial abundance
initial_abundance<- rep(0, n)
initial_abundance[1]<-1000 #start with input of eggs and remove burn in after reaching SSD


#Create parameter data
at0 = c(MT_at0 * 0.5, MT_at0 * 0.75, MT_at0) #25 and 50% reduction from metabolically optimal
z0 = c(MT_z0 * 0.75, MT_z0) #25% from metabolically optimal
Params2test<-crossing(at0, z0)
Params2test$Growth<- c( "50red", "50red", "25red", "25red", "MTG", "MTG")
Params2test$Survival<- c("MTS", "LowS", "MTS", "LowS", "MTS", "LowS")


################################################################################

#LAMBDAS FIXED TEMP
#Temp 0
LT0<-calc_fxd_temps_lambda(fxd_Temp = 0,  Parameter_data = Params2test)
#Temp 2.5
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
#FLUCTUATING LAMBDAS

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
#combine into one data frame for plotting

#add column where lambda is greater than 1 (rmax greater than 0)
fixd_temps_lambdas<-fixd_temps_lambdas[order(fixd_temps_lambdas$Param_combo),]
fixd_temps_peaks<-fixd_temps_peaks[order(fixd_temps_peaks$Param_combo),]
fixd_temps_peaks$lambda_1<- fixd_temps_lambdas$Lambda >= 1

lambdas_fluctuating<-lambdas_fluctuating[order(lambdas_fluctuating$Param_combo),]
vary_temps_peaks<-vary_temps_peaks[order(vary_temps_peaks$Param_combo),]
vary_temps_peaks$lambda_1<- lambdas_fluctuating$Lambda >= 1


#rmax into into one data used for plotting
data2plot<-rbind(fixd_temps_lambdas, lambdas_fluctuating)

