##################################################################################
############################ CODING BLOCK ########################################
##################################################################################
# Run-Time Environment:   R version 3.5.1
# Author:			            M.J. van Esdonk, LACDR
# Short title:			      Pharmacokinetic simulations in R 
# Version:			          V.1.0
# Model based on:
##################################################################################
##################################################################################
###
rm(list = ls(all.names = TRUE))
###

version <- 1.0


###############################################
############ Load libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)


###############################################
### Dosing

oral <- T

dose <- 100 # mg


###############################################
## Parameter estimates 
ka  <- 1 # Absorption rate constant
cl  <- 1 # Clearance
vd  <- 1 # Central distribution volume
vd2 <- 1 # Peripheral distribution volume
q1  <- 1 # Inter-compartmental clearance (set to 0 for 1-cmt model)




###############################################
## Insert etas and sigmas

etaka  <- 0.09
etacl  <- 0.09
etavd  <- 0.09
etavd2 <- 0.09
etaq1  <- 0.09




#################################################
## Insert residual error

sigmaprop <- 0.01 # Proportional error
sigmaadd  <- 0 # Additive error




###############################3
## Simulation info

nsamples <- 1000 ### Number of simulated individuals
sim_time <- 12 ## Time of simulation



# set probabilities of ribbon in figure
minprobs=0.1
maxprobs=0.9




###################################################
############# Set Dosing objects

if(oral){
  ## Oral dose
  Administration <-  as.data.frame(ev(ID=1:nsamples,ii=24, cmt=1, addl=9999, amt=dose, rate = 0,time=0)) 
}else{
  ## IV BOLUS
  Administration <-  as.data.frame(ev(ID=1:nsamples,ii=24, cmt=2, addl=9999, amt=dose, rate = 0,time=0)) 
}



## Sort by ID and time
data <- Administration[order(Administration$ID,Administration$time),]



######################################################################
### Load in model for mrgsolve
mod <- mread_cache("popPK")



## Specify the omegas and sigmas matrices
omega <- cmat(etaka,
              0, etacl,
              0,0,etavd,
              0,0,0,etavd2,
              0,0,0,0,etaq1)


sigma <- cmat(sigmaprop,
              0,sigmaadd)



## Set parameters in dataset
data$TVKA <- ka
data$TVCL <- cl
data$TVVC <- vd
data$TVVP1 <- vd2
data$TVQ1 <- q1




#################################################################
###### Perform simulation
out <- mod %>%
  data_set(data) %>%
  omat(omega) %>%
  smat(sigma) %>%
  mrgsim(end=sim_time,delta=sim_time/100, obsonly=TRUE) # Simulate 100 observations, independent of the total simulated time


### Output of simulation to dataframe
df <- as.data.frame(out)




#################################################################  
## Calculate summary statistics

sum_stat <- df %>%
  group_by(time) %>%
  summarise(Median_C=median(DV),
            Low_percentile=quantile(DV,probs=minprobs),
            High_percentile=quantile(DV,probs=maxprobs)
  )  




################################### Graphical output

ggplot(sum_stat, aes(x=time,y=Median_C)) +
  
  ## Add ribbon for variability
  geom_ribbon(aes(ymin=Low_percentile, ymax=High_percentile, x=time), alpha = 0.15, linetype=0)+
  
  ## Add median line
  geom_line(size=2) +
  
  # scale_y_log10()+
  
  # Set axis and theme
  ylab(paste("Concentration",sep=""))+
  xlab("Time after dose (h)")+
  theme_bw()+
  
  # Remove legend
  theme(legend.position="none")


###### Monte Carlo simulation
  
#Every curve has an equation. The equation for the PK curve can be found via excel trendline
#The equation for this curve was:
#y = -0.0015x6 + 0.0601x5 - 0.9303x4 + 7.0138x3 - 25.866x2 + 37.139x + 7.482
#Now that we have the equation, we can perform monte carlo simulation

f <- function(x) -0.0015*x^6 + 0.0601*x^5 - 0.9303*x^4 + 7.0138*x^3 - 25.866*x^2 + 37.139*x + 7.482

#Plotting the curve again via the equation

curve(f, from =0, to =8) 
library(DescTools)  
Shade(f, breaks = c(0,8), col = "red") # Shading the required region

###1st Monte carlo estimation

query_points <- runif(n=2000, min = 0, max = 8) #sampling 2000 points between 0 & 8
8* mean(f(query_points)) #Applying our function f to the random points and compute the corresponding mean
#here 8 is the length of interval

#Now lets see if the number of random points we sampled affects the accuracy of our estimate

area_estimates <- vector (length = 10000)
for (i in 1:10000) {
  query_points <- runif(n = i,min = 0, max = 8)
  area_estimates[i] <- 8* mean(f(query_points))
}
plot(area_estimates)
abline(h= 95, col= "red", lwd =2)

