#Atlantic salmon parameters

minmass<-0.03
maxmass<-20000

timescale = 1 #1 day

spawningmonth = 11 #November
spawningday = 1 #1st of the month

#Data needed for fecundity function
adult_mass_at_maturity = 1600 #g
fecundity_lm_intercept = 7.260878
fecundity_lm_slope = 1.012753

mean_egg_size = 0.125
egg_size_sd = 0.2

#Data needed for growth
#activation energy
E=0.45
#at0 varied
MT_at0<- 0.0145767684333848 #Initial parameter estimate
  
#Survival
T_pk = 9.05793 #t peak
power = 0.25 #quarter power
E_D = 2.10602 #deactivation energy
#z0 varied
MT_z0<-0.00942460726973781 #Initial parameter estimate
