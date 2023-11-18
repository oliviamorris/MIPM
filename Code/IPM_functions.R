#IPM functions for changing at0 and z0 parameters

#Global IPM params
# integration limits 
min.size<- log(0.03) 
max.size<- log(20000) #10000
# number of cells in the discretized kernel
n=300 
# boundary points (the edges of the cells defining the kernel)
b= seq(min.size, max.size, length.out = n+1)
# mesh points (midpoints of the cells)
y=0.5*(b[1:n]+b[2:(n+1)])
# width of the cells
h=y[2]-y[1] #same as (log(max.size) - log(min.size))/n



#Fecundity functions

#Number of eggs
f.yx <- function(xp, x, adult_mass = 1600, intercept = 7.260878, slope = 1.012753){  #m: mass(g), intercept and slope: from number of eggs produced per mass model
  x<-exp(x) 
  xp<- exp(xp)
  #create empty vector to populate
  egg <- rep(0, length(x))
  #sequence along mass to calculate number of eggs produced if over or equal to adult mass size
  for (i in seq_along(x)){
    if (x[i] >= adult_mass){
      m <- x[i]/1000     #convert mass as fecundity relationship is in kg and mass is in g
      egg_number <- exp(intercept)*(m^slope) 
      female_egg_number<- egg_number/2 #divide egg number by two to account for 1:1 sex ratio to model females
      egg[i]<- female_egg_number 
    } else {
      egg[i] <- 0
    }
  }
  return(egg)
}

#Size dist of egg sizes: input in log space
d.yx<- function(x, xp){ 
  dnorm(x, mean = log(0.125), sd = 0.2) 
} 



#Growth functions
k = 8.62*(10^(-5)) # gobal param: boltzmann's constant

#a. Development time from mass NEW
development_time <- function(m, Tc, at0, E=0.45, m0 = 0.03, M = 20000){ #m: mass(g), Tc: temp degrees celsius, t0: normalisation constant, a: slope, alpha=-E/kT20
  alpha = E/(k*(273^2))
  t <-  log( (1-( (m0/M)^(0.25) ))/ ( 1- ((m/M)^(0.25)) ) ) * (4 * (M^0.25) / ( at0 * exp(alpha* ( Tc/(1+(Tc/273))) )) ) 
  return(t)   #time in days to reach mass m
}

#b. Mass reached over developmental time NEW
development_mass <- function(t, Tc, at0, E = 0.45, m0 = 0.03, M = 20000){ #t: time in days from above function, Tc: temp degrees celsius, t0: normalisation constant, a: slope, alpha=-E/kT_0^2
  k <- 8.62*(10^(-5))
  alpha = E/(k*(273^2)) 
  a_temp<- at0 * exp(alpha* (Tc/(1+(Tc/273))) )
  m <- M * ( 1- (1-((m0/M)^0.25)) * (exp(-a_temp*t/(4*(M^0.25))) )) ^4
  return(m) #mass reached after amount of time t
}

#c. size next
sizeNext<-function(mass_vector, Tc, t, at0){ #Tc = temp degrees cels, t = time in days
  time <- development_time(mass_vector, Tc, at0)  + t
  mass_vector_next<- development_mass(time, Tc, at0)
  return(mass_vector_next)
}

#The growth function with normal distr and sd, input in log space
g.yx=function(x, xp, Tc , t, sd = h, at0) { 
  #exponentiate masses 
  mp<-exp(xp)
  #sizeNext under Tc after t days
  mass_t<-sizeNext(mp, Tc, t, at0)
  #Log sizeNext 
  log_mass_t<-log(mass_t)
  #create normal distribution 
  dnorm(x, mean = log_mass_t, sd =  sd) 
}



#U-shaped mortality
Inv_sharpe_schoolfield_edited<- function(Tc, m, T_pk, z0, E, E_D, power){ #K is specified as global param for all other eqns
  #Convert input temp to kelvin
  Temp <- Tc+273.15 
  
  #Reference temp
  T_ref<- 0+273.15 
  #T_pk has to be in kelvs too 
  T_pk_kel<-T_pk+273.15
  #Eqn
  Inv_S_S<- ((z0) * m^ (-power))* (( ( ( ( - E * exp( (E_D*(Temp-T_pk_kel))/(Temp*T_pk_kel*k) )  + E - E_D)  )* exp( - (( E*(Temp-T_ref))/ (Temp * T_ref * k)) ) ) /  (  E-E_D )))  
  
  return(Inv_S_S)
}

s.x<- function(x, Tc, z0, E = 0.45, T_pk = 9.05793, power = 0.25, E_D = 2.10602){
  m<- exp(x) #exponentiate masses
  #daily mortality rate
  Z<-Inv_sharpe_schoolfield_edited(Tc, m, T_pk, z0, E, E_D, power)
  #Convert to a daily survival probability
  S_prob<- exp(-(Z)) #using exponential CDF
  
  return(S_prob)
} 




#Start with normal day daily P
create_P_kernel<- function(y, h, Tc, t = 1, at0 = at0, z0 = z0){ #specify min and max mass, n number of mesh points, exp-mesh = yes/no?
  
  #Growth matrix G
  G = h*outer(y, y, g.yx, Tc, t, at0 = at0)  
  # survival S
  S = s.x(y, Tc, z0 = z0)
  #these give us our P matrix:
  P = G       
  for(i in 1:n) P[,i]=G[,i]*S[i]  
  
  return(P)
}

create_P_month<- function(y, h, month, Temp_vector, leapyear, at0 = at0, z0 = z0){  # month as a number e.g. Jan = 1, Feb = 2... etc, Tempvector data has be of year length and start from julian day 1 (01/01/YYY)
  #The Julian day at the beginning of each month 
  if (leapyear == "yes"){
    jday =  c(1,  32,  61,  92, 122, 153, 183, 214, 245, 275, 306, 336, 366)
  }else if (leapyear == "no") {
    jday =  c(1,  32,  60,  91, 121, 152, 182, 213, 244, 274, 305, 335, 365)
  } else (stop("check leapyear"))
  
  #Temperature vector depending on 
  if (month == 12){
    month_temps<- Temp_vector[(jday[month]):(jday[month+1])] #the final day in the vector
  }else  {
    month_temps<- Temp_vector[(jday[month]):(jday[month+1]-1)]
  }
  
  #Create survival matrices
  P_day <- vector("list",length(month_temps))
  for (i in seq_along(month_temps)){
    P_day[[i]]<- create_P_kernel (y, h, Tc = month_temps[i], at0 = at0, z0 = z0)
  }
  
  #Use Reduce function reduce to multiply all elements
  month_matrix<-Reduce("%*%", P_day)
  
  return(month_matrix)
}



#Reproduction day K= P+R
create_K_kernel<- function(y, h, Tc, t = 1, at0 = at0, z0 = z0){ #specify min and max mass, n number of mesh points, exp-mesh = yes/no?
  
  #Growth matrix
  G = h*outer(y,y,g.yx, Tc, t, at0 = at0)  
  # survival 
  S = s.x(y, Tc = Tc, z0 = z0)
  #these give us our P matrix:
  P = G          
  for(i in 1:n) P[,i]=G[,i]*S[i]  
  #Fecundity matrix
  F=outer(y, y, f.yx)   # reproduction kernel
  #Eggs
  D = h*outer(y,y,d.yx)
  #Number and size of eggs
  R= F*D
  #the whole kernel
  K = P+R 
  return(K)
}


#fecundity on one day of the month, later specified as 1st day here
create_K_month<- function(y, h, Temp_vector, leapyear, month, spawningday, at0 = at0, z0 = z0){  # month as a number e.g. Jan = 1, Feb = 2... etc, Tempvector data has be of year length and start from julian day 1 (01/01/YYY)
  #The Julian day at the beginning of each month 
  if (leapyear == "yes"){
    jday =  c(1,  32,  61,  92, 122, 153, 183, 214, 245, 275, 306, 336, 366)
  }else if (leapyear == "no") {
    jday =  c(1,  32,  60,  91, 121, 152, 182, 213, 244, 274, 305, 335, 365)
  } else (stop("check leapyear"))
  
  #Temperature vector depending on 
  if (month == 12){
    month_temps<- Temp_vector[(jday[month]):(jday[month+1])] #the final day in the vector
  }else  {
    month_temps<- Temp_vector[(jday[month]):(jday[month+1]-1)]
  }
  
  #Create survival matrices
  P_day <- vector("list",length(month_temps))
  for (i in seq_along(month_temps)){
    P_day[[i]]<- create_P_kernel(y, h, Tc = month_temps[i], at0 = at0, z0 = z0)
  }
  
  #Spawning day
  spawning_day_temp<-Temp_vector[(jday[month]+spawningday-1)]
  K_day<- create_K_kernel(y, h, Tc = spawning_day_temp, at0 = at0, z0 = z0)
  
  #Replace P matrix in P_day with K matrix created
  P_day[[spawningday]]<- K_day
  
  #Use Reduce function reduce to multiply all elements
  month_matrix<-Reduce("%*%", P_day)
  
  return(month_matrix)
}


#Create month matrices
project_monthmatrices<- function (y, h, n, Temp_vector_data, initial_abundance, spawningmonth = 11 , spawningday = 1, at0 = at0, z0 = z0){
  
  start_time<-Sys.time()

  years<- unique(Temp_vector_data$no.years.for.sim)
  
  #Convert input years to months
  total_period_projection_months<-max(years)* 12 
  #Create empty matrix to populate (abundance data)
  allabundances<- matrix(0, nrow = n, ncol = total_period_projection_months+1) #plus one for initial abundance
  #Insert initial abundance
  allabundances[,1]<- initial_abundance
  
  #create vector to put all months into
  allmonths<-vector("list", max(years))
  
  #time frame beginning with spawning month
  if (spawningmonth==1){
    time <- 1:12
  } else if (spawningmonth %in% 2:12) {
    time<- c(spawningmonth:12,1:(spawningmonth-1))
  }
  
  #leap year?
  for (i in seq_along(years)){
    
    data2subset<-subset(Temp_vector_data, no.years.for.sim == i) #novyears
    Temp_vector<- data2subset$meantemp 
    
    if (nrow(data2subset) == 366){
      leapyear <- "yes"
    } else if (nrow(data2subset) == 365){
      leapyear <- "no"
    } else {
      stop("wrong number of days in a year in temperature data")
    }
    
    #create months list 
    months<- list()
    
    #spawning month as first month (IS THIS CORRECT OR SHOULD I AUTOMATE IT TO BE EASILY SPECIFIED?)
    months[[1]] <-create_K_month(y, h, Temp_vector = Temp_vector, leapyear = leapyear, month = spawningmonth, spawningday = spawningday, at0 = at0, z0 = z0)
    #rest of months
    for (j in 2:12){
      months[[j]] <- create_P_month(y, h, month = time[j] , Temp_vector = Temp_vector, leapyear = leapyear, at0 = at0, z0 = z0)
    }
    
    #months into a year
    allmonths[[i]]<- months
    
  }
  
  #Un list allmonths to loop over
  allmonths<-lapply(rapply(allmonths, enquote, how="unlist"), eval)
  
  #Project loop
  for (i in 2:(total_period_projection_months+1)){
    allabundances[,i]<- allmonths[[i-1]] %*% allabundances[,i-1]
  }

  end_time<- Sys.time()
  #How long function takes to run
  print(end_time - start_time)
  
  return(allabundances)
}



#DATA
#function that adds burn in time by repeating first year of data n times
add_burn_in<- function(data, n_burnin_years){
  #subset first years data
  firstyear_data<- subset(data, novyears == 1)
  #repeat first year data
  firstyear_data_burnin<-cbind(firstyear_data[rep(1:nrow(firstyear_data), n_burnin_years), ], no.years.for.sim = rep(1:n_burnin_years, each = nrow(firstyear_data)))
  #Add no of years for sim column
  data$no.years.for.sim<- data$novyears+n_burnin_years
  #combine the burn in data to the actual data set
  Data2Run<-rbind(firstyear_data_burnin, data)
  
  return(Data2Run) #returns dataframe added in with burn in added on
}


# 
# ###### a modification from projectmonthmatrices
# #Create month matrices
# get_month_matrices<- function (y, h, n, Temp_vector_data, spawningmonth = 11 , spawningday = 1, at0 = at0, z0 = z0){
#   
#   start_time<-Sys.time()
#   
#   years<- unique(Temp_vector_data$no.years.for.sim)
#   
#   #create vector to put all months into
#   allmonths<-vector("list", max(years))
#   
#   #time frame beginning with spawning month
#   if (spawningmonth==1){
#     time <- 1:12
#   } else if (spawningmonth %in% 2:12) {
#     time<- c(spawningmonth:12,1:(spawningmonth-1))
#   }
#   
#   #leap year?
#   for (i in seq_along(years)){
#     
#     data2subset<-subset(Temp_vector_data, no.years.for.sim == i) #novyears
#     Temp_vector<- data2subset$meantemp #Temp_vector_data needs to have meantemp column
#     
#     if (nrow(data2subset) == 366){
#       leapyear <- "yes"
#     } else if (nrow(data2subset) == 365){
#       leapyear <- "no"
#     } else {
#       stop("wrong number of days")
#     }
#     
#     #create months list 
#     months<- list()
#     
#     #spawning month as first month (IS THIS CORRECT OR SHOULD I AUTOMATE IT TO BE EASILY SPECIFIED?)
#     months[[1]] <-create_K_month(y, h, Temp_vector = Temp_vector, leapyear = leapyear, month = spawningmonth, spawningday = spawningday, at0 = at0, z0 = z0)
#     #rest of months
#     for (j in 2:12){
#       months[[j]] <- create_P_month(y, h, month = time[j] , Temp_vector = Temp_vector, leapyear = leapyear, at0 = at0, z0 = z0)
#     }
#     
#     #months into a year
#     allmonths[[i]]<- months
#     
#   }
#   
#   #Un list allmonths to loop over
#   allmonths<-lapply(rapply(allmonths, enquote, how="unlist"), eval)
#   #allmonths<-unlist(allmonths, recursive = F)
#   
#   end_time<- Sys.time()
#   #How long function takes to run
#   print(end_time - start_time)
#   
#   return(allmonths) #allmonths
# }
# 
# 
# 
# 
# #ouput from above goes into this to get matrix out
# 
# Project_abundances_matrix<- function(allmonths, Temp_vector_data, initial_abundance){
#   
#   #From temperature data input
#   years<- unique(Temp_vector_data$no.years.for.sim)
#   #Convert input years to months
#   total_period_projection_months<-max(years)* 12 
#   #Create empty matrix to populate (abundance data)
#   allabundances<- matrix(0, nrow = n, ncol = total_period_projection_months+1) #plus one for initial abundance
#   #Insert initial abundance
#   allabundances[,1]<- initial_abundance
#   
#   
#   #Project loop
#   for (i in 2:(total_period_projection_months+1)){
#     allabundances[,i]<- allmonths[[i-1]] %*% allabundances[,i-1]
#   }
#   
#   return(allabundances)
# }



# 
# #NEW WITH PROPORTIONAL W_ij IN RELATION TO WHERE DATA IS
# #2. Function that uses project_matrices() and spits out just the proportions in each august in each year
# SAD_props_count_data_month_nogrouping<- function(abundancemodeldata, countdata, burnindates, modeldates, startdate, month2select){
#   
#   #start_time<-Sys.time()
#   
#   #IPM over data time frame called abundance model data
#   #abundancemodeldata<- project_monthmatrices(y = y, h = h, n = n, initial_abundance = initial_abundance, Temp_vector_data = Data2Run)
#   
#   #Convert matrix to data frame with dates and abundances to subset August (make a date argument in function?)
#   #Save length of abundance_data (long format)
#   length<- ncol(abundancemodeldata)
#   #Convert into data frame
#   abundancemodeldata<-data.frame(abundancemodeldata)
#   #Add log mass column
#   abundancemodeldata$logmass<- y
#   #Add exp mass column
#   abundancemodeldata$mass<- exp(y)
#   
#   #Wide to long format
#   abundance_data<-gather(abundancemodeldata, col, abundance, 1:length, factor_key = T)
#   
#   #Number of data divisions = n in global environment
#   #Add dates to data
#   abundance_data_df_allmonthyears_noinit<-abundance_data[-c(1:n), ]
#   abundance_data_df_allmonthyears_noinit$Date <-c(rep(burnindates, each = n), rep(modeldates, each = n))
#   abundance_data_df_allmonthyears_noinit$Date <-as.yearmon(abundance_data_df_allmonthyears_noinit$Date)
#   abundance_data_df_allmonthyears_noinit_monthyear<-abundance_data_df_allmonthyears_noinit  %>% group_by(mass, abundance) %>%
#     mutate(month = format(Date, "%m"), year = format(Date, "%Y" ))
#   
#   #Subset the august dates
#   month_abundance_subsetted_year<- subset(abundance_data_df_allmonthyears_noinit_monthyear, month == month2select)
#   
#   #Remove burnin
#   months_years_to_test<-subset(month_abundance_subsetted_year, Date >= startdate) #specify date?
#   
#   #Create one data frame to run analyses on
#   data<-cbind(countdata, months_years_to_test)
#   #Add number of years to loop over
#   data$noyears<-as.numeric(factor(data$year))
#   
#   #make count data 0s for fish that are catchable, between exp(-1) ~0.36g and exp(3)
#   data$count[is.na(data$count) & data$mass > exp(-1) & data$mass < exp(3)] <- 0
#   
#   #Calculate proportions from abundances, proportions only include rows where there is data (excludes where count is NA)
#   #Sum
#   data_total<-data %>% group_by(year) %>%  mutate(total = sum(abundance * !is.na(count))) 
#   #Divide proportions by total
#   data_props<-data_total %>% group_by(year, mass) %>%  mutate(prop = (abundance * !is.na(count))/total)
#   
#   #Replace prop 0s with NA to remove
#   data_props$prop[data_props$prop==0] <- NA
#   data4LL<-data_props %>% drop_na(prop)
#   
#   #return Data to run log likelihood on
#   return(data4LL)
# }
# 
# 
# 
# 
# #trying diff param space
# paramsearch_juveloglike_lambda_slope_intercept<- function(Paramdata, Temp_vector_data, initial_abundance, countdata, n_burnin_years = 40, burnindates, modeldates, startdate, month2select, at0 = at0, z0= z0, spawningmonth){
#   
#   #Create empty data
#   out <- matrix(nrow = nrow(Paramdata), ncol = 2, NA)
#   #and empty list for slope calc (size frequency dist)
#   out_list<-list()
#   #Loop for every parameter combination
#   for (i in 1:nrow(Paramdata)) {
#     #Print count for looking at how far through we are
#     #print(paste(i)) 
#     at0<- Paramdata$at0[i] 
#     z0<- Paramdata$z0[i]
#     #print(paste(at0, z0))
#     #Create IPM over data time frame called abundance model data
#     #abundancemodeldata<- project_monthmatrices(y = y, h = h, n = n, initial_abundance = initial_abundance, Temp_vector_data = Data2Run)
#     
#     #Get month matrices over which to calculate lambda
#     allmonthsdata<-get_month_matrices(y = y, h = h, n = n, Temp_vector_data = Temp_vector_data, at0 = at0, z0 = z0, spawningmonth = spawningmonth)
#     #Project these as usual
#     abundancemodeldata<-Project_abundances_matrix(allmonths = allmonthsdata, Temp_vector_data = Temp_vector_data, initial_abundance = initial_abundance)
#     #Turn into proportions and combine with empirical data
#     SAD_props_count<-SAD_props_count_data_month_nogrouping(abundancemodeldata = abundancemodeldata, countdata = countdata, burnindates = burnindates, modeldates = modeldates, startdate = startdate, month2select = month2select)
#     #Fill data frame with log likelihood
#     out[i,1] <- overall_loglik(data = SAD_props_count)
#     
#     #remove burnin for calculating lambda
#     burnin2remove<-1:(n_burnin_years*12)
#     allmonthsdata_noburnin<-allmonthsdata[-burnin2remove]
#     #Fill data with Lambda (no burn in)
#     out[i,2] <- get_lambda(allmonths = allmonthsdata_noburnin)
#     #Slope
#     out_list[[i]] <- calc_slope_allyears(abundance_data = abundancemodeldata, burnindates = burnindates, modeldates = modeldates, startdate = startdate)
#   }
#   
#   
#   data.out <- as.data.frame(out)
#   names(data.out) <- paste(c("JuveLogLikelihood","Lambda"))
#   data.out <- cbind(Paramdata, data.out)
#   
#   #slope
#   #Calculate slope and add to dataframe
#   slope<- c()
#   for (i in 1:length(out_list)){
#     #slope[i]<-coef(out_list[[i]])[2]
#     slope[i]<-tryCatch(coef(out_list[[i]])[2], error=function(err) NA)  
#   }
#   data.out$Slope<- slope
#   
#   #Add intercept to dataframe
#   intercept<- c()
#   for (i in 1:length(out_list)){
#     #intercept[i]<-coef(out_list[[i]])[1]
#     intercept[i]<-tryCatch(coef(out_list[[i]])[1], error=function(err) NA)  
#   }
#   data.out$Intercept<- intercept
#   
#   
#   return(data.out)
#   
# }
# 

#Convert to dataframe for plotting
matrix2dataframe<- function(abundancemodeldata, countdata, burnindates, modeldates, startdate){
  
  #start_time<-Sys.time()
  
  #IPM over data time frame called abundance model data
  #abundancemodeldata<- project_monthmatrices(y = y, h = h, n = n, initial_abundance = initial_abundance, Temp_vector_data = Data2Run)
  
  #Convert matrix to data frame with dates and abundances to subset August (make a date argument in function?)
  #Save length of abundance_data (long format)
  length<- ncol(abundancemodeldata)
  #Convert into data frame
  abundancemodeldata<-data.frame(abundancemodeldata)
  #Add log mass column
  abundancemodeldata$logmass<- y
  #Add exp mass column
  abundancemodeldata$mass<- exp(y)
  
  #Wide to long format
  abundance_data<-gather(abundancemodeldata, col, abundance, 1:length, factor_key = T)
  
  #Number of data divisions = n in global environment
  #Add dates to data
  abundance_data_df_allmonthyears_noinit<-abundance_data[-c(1:n), ]
  abundance_data_df_allmonthyears_noinit$Date <-c(rep(burnindates, each = n), rep(modeldates, each = n))
  abundance_data_df_allmonthyears_noinit$Date <-as.yearmon(abundance_data_df_allmonthyears_noinit$Date)
  abundance_data_df_allmonthyears_noinit_monthyear<-abundance_data_df_allmonthyears_noinit  %>% group_by(mass, abundance) %>%
    mutate(month = format(Date, "%m"), year = format(Date, "%Y" ))
  
  #Remove burnin
  abundance_data_df_allmonthyears_noinit_monthyear<-subset(abundance_data_df_allmonthyears_noinit_monthyear, Date >= startdate) #specify date?
  
  
  return(abundance_data_df_allmonthyears_noinit_monthyear)
  
}




################################################################################
#Functions needed for running manuscript code

#so amplitude needs to be a function of meantemp
amplitude_func<- function(meantemp, x = 6.2285, y = -0.1227){
  amp<-( y * meantemp) + x
  return(amp)
}


#Create function to make fixed temp dataframe (fixed for 60 years for now)
create_fixed_temp_df<- function(fixdtemp){
  meantemp<-  rep(fixdtemp,365*70) 
  fixed_temp_data<-as.data.frame(meantemp)
  fixed_temp_data$no.years.for.sim<-rep(1:70, each =365)
  return(fixed_temp_data)
}

#fixed temps lambda
calc_fxd_temps_lambda<- function(fxd_Temp, Parameter_data){
  
  #Lambda and temp
  Lambda<- c()
  Temp<-fxd_Temp
  
  #create temp data frame
  fixed_temp_data = create_fixed_temp_df(fixdtemp = fxd_Temp)
  head(fixed_temp_data)
  
  allmonthsdata_0_PDT1<-get_month_matrices(z0 = Parameter_data$z0[1] , at0 = Parameter_data$at0[1],  y = y, h = h, n = n, Temp_vector_data = fixed_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT1_noburnin<-allmonthsdata_0_PDT1[-burnin2remove]
  Lambda[1]<-get_lambda(allmonths = allmonthsdata_0_PDT1_noburnin)
  
  allmonthsdata_0_PDT2<-get_month_matrices(z0 = Parameter_data$z0[2] , at0 = Parameter_data$at0[2],  y = y, h = h, n = n, Temp_vector_data = fixed_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT2_noburnin<-allmonthsdata_0_PDT2[-burnin2remove]
  Lambda[2]<-get_lambda(allmonths = allmonthsdata_0_PDT2_noburnin)
  
  allmonthsdata_0_PDT3<-get_month_matrices(z0 = Parameter_data$z0[3] , at0 = Parameter_data$at0[3],  y = y, h = h, n = n, Temp_vector_data = fixed_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT3_noburnin<-allmonthsdata_0_PDT3[-burnin2remove]
  Lambda[3]<-get_lambda(allmonths = allmonthsdata_0_PDT3_noburnin)
  
  allmonthsdata_0_PDT1<-get_month_matrices(z0 = Parameter_data$z0[4] , at0 = Parameter_data$at0[4],  y = y, h = h, n = n, Temp_vector_data = fixed_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT1_noburnin<-allmonthsdata_0_PDT1[-burnin2remove]
  Lambda[4]<-get_lambda(allmonths = allmonthsdata_0_PDT1_noburnin)
  
  allmonthsdata_0_PDT2<-get_month_matrices(z0 = Parameter_data$z0[5] , at0 = Parameter_data$at0[5],  y = y, h = h, n = n, Temp_vector_data = fixed_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT2_noburnin<-allmonthsdata_0_PDT2[-burnin2remove]
  Lambda[5]<-get_lambda(allmonths = allmonthsdata_0_PDT2_noburnin)
  
  allmonthsdata_0_PDT3<-get_month_matrices(z0 = Parameter_data$z0[6] , at0 = Parameter_data$at0[6],  y = y, h = h, n = n, Temp_vector_data = fixed_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT3_noburnin<-allmonthsdata_0_PDT3[-burnin2remove]
  Lambda[6]<-get_lambda(allmonths = allmonthsdata_0_PDT3_noburnin)
  
  
  LambdaTempdf<-cbind(Parameter_data, Lambda, Temp)
  
  return(LambdaTempdf)
} #FOR SIX PARAMS

#median
weighted.median <- function(x, w) {
  w <- w[order(x)]
  x <- x[order(x)]
  
  prob <- cumsum(w)/sum(w)
  ps <- which(abs(prob - .5) == min(abs(prob - .5)))
  return(x[ps])
}

#august sixe distribution median mass
aug_median_mass_fxdtemp<- function(fxd_Temp, at0, z0){
  
  #create temp data frame
  fixed_temp_data <- create_fixed_temp_df(fixdtemp = fxd_Temp)
  #head(fixed_temp_data)
  
  abundance_data<-project_monthmatrices(y = y, h = h, n = n, initial_abundance = initial_abundance, Temp_vector_data = fixed_temp_data, at0 = at0, z0 = z0, spawningmonth = 11)
  #plot(colSums(abundance_data), type ="l" , log = "y")
  #Save length of abundance_data (long format)
  length<- ncol(abundance_data)
  #Into data frame
  abundance_data<-data.frame(abundance_data)
  abundance_data$logmass<- y
  #Wide to long format
  abundance_data<-gather(abundance_data, col, abundance, 1:length , factor_key = T)
  #Remove first 40 years
  abundance_data_noinit<-abundance_data[-c(1:144300), ]
  #add month and years for 30 years
  months.<-rep(c(11,12, 1:10), each = 300)
  abundance_data_noinit$month<-rep(months., 30)
  abundance_data_noinit$years<- rep(1:30, each = 300*12) 
  DF_fixed<-abundance_data_noinit
  DF_fixed$month<- factor(DF_fixed$month, levels = c("11", "12", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10") )
  
  #Select August
  DF_fixed_aug<-subset(DF_fixed, month == "8")
  
  #Stable Age distribution as graph?
  
  #need to time average over 'model' time 
  #group by mass calculate mean abundance
  avprop_overTS<-DF_fixed_aug %>% group_by(logmass) %>% summarise(average_ts = mean(abundance))
  
  avprop_overTS$total<-sum(avprop_overTS$average_ts)
  avprop_overTS$averageprop_ts<- avprop_overTS$average_ts / avprop_overTS$total
  
  
  #plot
  p<-ggplot(avprop_overTS, aes(logmass, averageprop_ts))+
    geom_line()+
    labs(title = "Model: Time averaged proportions", y = "Average proportion over years")+
    theme_cowplot()
  print(p)
  
  #aug_median_mass<-avprop_overTS$logmass[which.min(abs(avprop_overTS$averageprop_ts - median(avprop_overTS$averageprop_ts)))]
  aug_median_mass<-weighted.median(avprop_overTS$logmass, avprop_overTS$averageprop_ts)
  
  #aug_mean_mass<- weighted.mean(avprop_overTS$logmass, avprop_overTS$averageprop_ts)
  
  #stretchiness<- aug_mean_mass - aug_median_mass
  
  return(aug_median_mass)
}

#remove 40 years (480 momnths) burn in time
burnin2remove<-1:(40*12)

#amplitudes relationship to temp
amplitude_func<- function(meantemp, x = 5.4265, y = -0.0163){   #x = 6.2285, y = -0.1227
  amp<-( y * meantemp) + x
  return(amp)
}

#fluctuations lambdas
calc_varied_temps_lambda<- function(vary_Temp, Parameter_data){
  
  #flucs
  sine_x<- seq(1, 365*70, by = 1)
  sine_data<-as.data.frame(sine_x)
  sine_y<- amplitude_func(vary_Temp)*sin((2*pi*sine_x)/365-2.5)+vary_Temp
  sine_data$sine_y <- sine_y
  
  #plot(sine_x,sine_y,type="l", ylab = "Temp", xlab = "Days")
  sine_temp_data<-as.data.frame(sine_y)
  #need nov.years and no.years.for.sim
  #burn in uses nov years to calculate no years for sim so just need no years for sim in thsi bit
  sine_temp_data$no.years.for.sim<-rep(1:70, each =365)
  sine_temp_data$meantemp<- sine_y
  
  
  #Lambda and temp
  Lambda<- c()
  Temp<-vary_Temp
  
  #Lambda[1]
  allmonthsdata_0_PDT1<-get_month_matrices(z0 = Parameter_data$z0[1] , at0 = Parameter_data$at0[1],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT1_noburnin<-allmonthsdata_0_PDT1[-burnin2remove]
  Lambda[1]<-get_lambda(allmonths = allmonthsdata_0_PDT1_noburnin)
  
  allmonthsdata_0_PDT2<-get_month_matrices(z0 = Parameter_data$z0[2] , at0 = Parameter_data$at0[2],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT2_noburnin<-allmonthsdata_0_PDT2[-burnin2remove]
  Lambda[2]<-get_lambda(allmonths = allmonthsdata_0_PDT2_noburnin)
  
  allmonthsdata_0_PDT3<-get_month_matrices(z0 = Parameter_data$z0[3] , at0 = Parameter_data$at0[3],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT3_noburnin<-allmonthsdata_0_PDT3[-burnin2remove]
  Lambda[3]<-get_lambda(allmonths = allmonthsdata_0_PDT3_noburnin)
  
  allmonthsdata_0_PDT4<-get_month_matrices(z0 = Parameter_data$z0[4] , at0 = Parameter_data$at0[4],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT4_noburnin<-allmonthsdata_0_PDT4[-burnin2remove]
  Lambda[4]<-get_lambda(allmonths = allmonthsdata_0_PDT4_noburnin)
  
  allmonthsdata_0_PDT4<-get_month_matrices(z0 = Parameter_data$z0[5] , at0 = Parameter_data$at0[5],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT4_noburnin<-allmonthsdata_0_PDT4[-burnin2remove]
  Lambda[5]<-get_lambda(allmonths = allmonthsdata_0_PDT4_noburnin)
  
  allmonthsdata_0_PDT4<-get_month_matrices(z0 = Parameter_data$z0[6] , at0 = Parameter_data$at0[6],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  allmonthsdata_0_PDT4_noburnin<-allmonthsdata_0_PDT4[-burnin2remove]
  Lambda[6]<-get_lambda(allmonths = allmonthsdata_0_PDT4_noburnin)
  
  # allmonthsdata_0_PDT4<-get_month_matrices(z0 = Parameter_data$z0[7] , at0 = Parameter_data$at0[7],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  # allmonthsdata_0_PDT4_noburnin<-allmonthsdata_0_PDT4[-burnin2remove]
  # Lambda[7]<-get_lambda(allmonths = allmonthsdata_0_PDT4_noburnin)
  # 
  # allmonthsdata_0_PDT4<-get_month_matrices(z0 = Parameter_data$z0[8] , at0 = Parameter_data$at0[8],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  # allmonthsdata_0_PDT4_noburnin<-allmonthsdata_0_PDT4[-burnin2remove]
  # Lambda[8]<-get_lambda(allmonths = allmonthsdata_0_PDT4_noburnin)
  # 
  # allmonthsdata_0_PDT4<-get_month_matrices(z0 = Parameter_data$z0[9] , at0 = Parameter_data$at0[9],  y = y, h = h, n = n, Temp_vector_data = sine_temp_data, spawningmonth = 11)
  # allmonthsdata_0_PDT4_noburnin<-allmonthsdata_0_PDT4[-burnin2remove]
  # Lambda[9]<-get_lambda(allmonths = allmonthsdata_0_PDT4_noburnin)
  # 
  LambdaTempdf<-cbind(Parameter_data, Lambda, Temp)
  
  return(LambdaTempdf)
} #FOR SIX PARAMS


aug_median_mass_vary_temp<- function(vary_Temp, at0, z0){
  
  #create temp data frame
  
  # #0-25
  # sine_x<- seq(1, 365*70, by = 1)
  # sine_data<-as.data.frame(sine_x)
  # sine_y<- amplitude_func(vary_Temp)*sin((2*pi*sine_x)/365-2.5)+vary_Temp
  # sine_data$sine_y <- sine_y
  # 
  # #plot(sine_x,sine_y,type="l", ylab = "Temp", xlab = "Days")
  # sine_temp_data<-as.data.frame(sine_y)
  # #need nov.years and no.years.for.sim
  # #burn in uses nov years to calculate no years for sim so just need no years for sim in thsi bit
  # sine_temp_data$no.years.for.sim<-rep(1:70, each =365)
  # sine_temp_data$meantemp<- sine_y
  
  #flucs
  sine_x<- seq(1, 365*70, by = 1)
  sine_data<-as.data.frame(sine_x)
  sine_y<- amplitude_func(vary_Temp)*sin((2*pi*sine_x)/365-2.5)+vary_Temp
  sine_data$sine_y <- sine_y
  
  #plot(sine_x,sine_y,type="l", ylab = "Temp", xlab = "Days")
  sine_temp_data<-as.data.frame(sine_y)
  #need nov.years and no.years.for.sim
  #burn in uses nov years to calculate no years for sim so just need no years for sim in thsi bit
  sine_temp_data$no.years.for.sim<-rep(1:70, each =365)
  sine_temp_data$meantemp<- sine_y
  
  abundance_data<-project_monthmatrices(y = y, h = h, n = n, initial_abundance = initial_abundance, Temp_vector_data = sine_temp_data, at0 = at0, z0 = z0, spawningmonth = 11)
  #plot(colSums(abundance_data), type ="l" , log = "y")
  #Save length of abundance_data (long format)
  length<- ncol(abundance_data)
  #Into data frame
  abundance_data<-data.frame(abundance_data)
  abundance_data$logmass<- y
  #Wide to long format
  abundance_data<-gather(abundance_data, col, abundance, 1:length , factor_key = T)
  #Remove first 40 years
  abundance_data_noinit<-abundance_data[-c(1:144300), ]
  #add month and years for 30 years
  months.<-rep(c(11,12, 1:10), each = 300)
  abundance_data_noinit$month<-rep(months., 30)
  abundance_data_noinit$years<- rep(1:30, each = 300*12) 
  DF_fixed<-abundance_data_noinit
  DF_fixed$month<- factor(DF_fixed$month, levels = c("11", "12", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10") )
  
  #Select August
  DF_fixed_aug<-subset(DF_fixed, month == "8")
  
  #Stable Age distribution as graph?
  
  #need to time average over 'model' time 
  #group by mass calculate mean abundance
  avprop_overTS<-DF_fixed_aug %>% group_by(logmass) %>% summarise(average_ts = mean(abundance))
  
  avprop_overTS$total<-sum(avprop_overTS$average_ts)
  avprop_overTS$averageprop_ts<- avprop_overTS$average_ts / avprop_overTS$total
  
  
  #plot
  p<-ggplot(avprop_overTS, aes(logmass, averageprop_ts))+
    geom_line()+
    labs(title = "Model: Time averaged proportions", y = "Average proportion over years")+
    theme_cowplot()
  print(p)
  
  #aug_median_mass<-avprop_overTS$logmass[which.min(abs(avprop_overTS$averageprop_ts - median(avprop_overTS$averageprop_ts)))]
  aug_median_mass<-weighted.median(avprop_overTS$logmass, avprop_overTS$averageprop_ts)
  
  return(aug_median_mass)
}

#################################################################################

