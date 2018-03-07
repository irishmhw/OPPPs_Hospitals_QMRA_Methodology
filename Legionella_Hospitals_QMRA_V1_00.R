#================================================================================================|
# ------------------------------------------- HEADER ------------------------------------------- |
# Legionella_Hospitals_QMRA_V1_00.R is an R source code developed based on legionell.draft.18.R  |
# developed during the 2017 QMRA III summer institute. Significant advancements and repairs to   |
# non-working components have been accomplished.                                                 |
# Legionella_Hospitals_QMRA_V1_00.R is written by Alexis L. Mraz MS and Mark H. Weir Ph.D. both  |
# working jointly for the Applied Research Center (ARC) of NSF International, Ann Arbor, MI, USA |
# and with the College of Public Health, The Ohio State University, Columbus, OH, USA            |
# Funding for this sourec code development is provided by NSF International, Ann Arbor, MI, USA  |
# QMRAIII group members include: Courtney Gardner, Rubayat Jamal, Emily Kumpel, Alexis Mraz,     |
# Joyjit Saha, and Amanda Wilson; Mentored by: Jade Mitchell Ph.D., MSU, Lansing, MI, USA, and   | 
# Kerry Hamilton Ph.D. Drexel University, Phila, PA, USA                                         | 
#------------------------------------------------------------------------------------------------|
# Legionella_Hospitals_QMRA_V1_00.R models the risks of infection, illness and assocaited DALYs  |
# from exposure to Legionella pneumophila in healthcare environments                             |
# The QMRA model this source code is for has been subsequently published in _____________________|
# REFERENCE DOI   |
#================================================================================================|


#Function Form - Exploring Interventions


rm(list=ls()) #clearing environment --> MHW not particularly necessary and removes any data you need or want in workspace, better at the end. 

#USING PACKAGE TRUNCNORM 
#library(truncnorm)
#======================================================================================================================|
# Set working directory and install required packages if not installed, truncnorm() is from original code, does not    |
# with older versions of R
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#setwd("U:/Publications/Legionella_2016_Case_Study/R_Codes/Non_Confounding_RScript/10000_Iter/Pop_testing/")          #|
if("truncnorm" %in% rownames(installed.packages()) == FALSE) {install.packages("truncnorm")} else{require(truncnorm)} #| 
source("tnormal.r") # truncnorm seems to be having problems working with older versions of R. This will always work   #|
#======================================================================================================================|



set.seed(24) #setting seed

######################################################################################

# Code for Function:

######################################################################################

#Explanation of arguments in function:
   # l.plumbing: length of pipe from the water heater to the shower pipe (m)
   # l.shower: length of pipe that may experience stagnation in between showers (m)
   # duration.shower: the length of the showering event (in minutes)
   # Temperature.waterheater: the temperature at which the water heater is set (deg C)
   # Desired.temp: temperature at which shower is desired (deg C)
   # scalding guard: 1 or 0 (1 = yes, scalding guard installed; 0 = no, scalding guard not installed)

#legionella.model<-function(l.plumbing=443,l.shower,duration.shower=15,Temperature.waterheater=59,desired.temp=43,scalding.guard=0){
WHtempoptions <- c(49,54,59)
lengths <- c(50, 200, 400)

#Loop for using various water heater temperatures (WHtempoptions) and various pipe lengths from shower to water
#heater (lengths) -> (3 temp x 3 pipe lengths) = 9 different scenarios

for(a in 1:3){ 
  Temperature.waterheater <- WHtempoptions[a] # water temperature loop

    for(b in 1:3){ #shower pipe length loop
  l.plumbing=lengths[b]
  
duration.shower=15 #duration of shower in minutes
loop.3<-duration.shower*60 #duration of shower in seconds

desired.temp=43 #desired temperature of shower water (in deg C)

  ##############################################################
  
  
  #######  PART 1: WATER HEATER
  
  
  ##############################################################
  residence.time<-2.1*60*60 #amount of time water is in heater before being distributed (in seconds)
  loops.wh<-residence.time #assuming residence time in water heater #NEED TO ASK EMILY AND COURTNEY ABOUT THIS ASSUMPTION
  matrix.wh<-matrix(nrow=loops.wh,ncol=3) #matrix for the output from the water heater portion of the model 
  # # rows = # seconds of residence time; 
  #columns = 2 (column 1: temperature in waterheater in non-growth region; column 2: temperature in water heater in growth region; column 3: concentration of L. pneumophila)
  p.1<-0.14 #proportion that is staying at original temp
  p.2<-0.86 #proportion undergoing cooling
  
  C.0<-rtruncnorm(1,0,Inf,10,8)     #incoming concentration of legionella (copies/mL) 
  C.r<-0    #incoming concentration of legionella from recirculating pipe (CFU/mL)
  C.hwh.2<-0.89/1000 # incoming concentration of legionella from municipal water -> median value from Borella et al., 2004 #cfu/mL
  
  for (i in 1:loops.wh){
    
    #Decay rate
    r.d<-runif(1,2.8*10^-6,8.3*10^-6) #min and max values found in literature - assumed a normal distribution
    
    #Temperature of Water - for non-growth region
    matrix.wh[i,1]<-Temperature.waterheater
    
    #Temperature of Water - for growth region
    matrix.wh[i,2]<-46 #assuming 46 deg C for now - (we acknowledge that this is a function of the set temperature of the water heater)
      #To be addressed! :) 
    
    #change in concentration based on growth rate (temperature dependent)
    if(matrix.wh[i,2]>=48){
      g.rate<-0 #/s
    }else if (matrix.wh[i,2]<48 & matrix.wh[i,2]>=40){
      g.rate<-0.162 #/s
    }else if (matrix.wh[i,2]<40 & matrix.wh[i,2]>=35){
      g.rate<-.048 #/s
    }else if (matrix.wh[i,2]<35 & matrix.wh[i,2]>=25){
      g.rate<-0.00139 #/s
    }else{
      g.rate<-0 #/s
    }
    
    
    #Legionella concentration
    
    ########
    
    #C.ng<- concentration of legionella in the non growth region of the water heater
    C.ng<-(C.0+C.r+C.hwh.2)-(C.0+C.r+C.hwh.2)*r.d #non growth concentration in water heater
    
    #no account for growth - concentration change per second is a function of original concentration (C.0),
    #incoming concentration from the recirculation pipe (C.r), incoming concentration from municipal water (C.hwh.2),
    # and loss from decay for the sum of these concentrations ((C.0+C.r+C.hwh.2)*r.d)
    
    ########
    
    #C.g<- concentraiton fo legionella in growth region of water heater
    C.g<-C.ng*exp(g.rate)-(C.ng*r.d) #growth concentration in water heater 
    
    #Account for growth and decay - concentration change per second is equal to the concentration in the non-growth section
    #times an exponential growth rate (C.ng*exp(g.rate)) minus any expected decay ((C.ng*r.d))
    #Assuming concentrations in the two sections of the water heater are equal. This is why the changes in 
    #C.g are related to growth and decay of C.ng concentration
    
    ########
    
    #total legionella concentration in water heater - a weighted sum of legionella concentration
    #in growth and non growth regions in water heater. We are assuming 14% of the water heater is not
    #at the "set temperature" therefore allow for growth to occur
    
    matrix.wh[i,3]<-(0.86*C.ng)+(0.14*C.g) #combined concentration in water heater 
      #Assuming 14% of water heater is not at "set temperature," therefore allowing for growth to occur 
    
  }
  print("Part 1 = check")
  ##############################################################
  
  
  #######  PART 2: PREMISE PLUMBING
  
  
  ##############################################################

loop.1<-l.plumbing #looping per 1m distance as opposed to per second in this part of the model, 
  #because equations used to calculate change in concentration were a function of distance from the water heater
  #Therefore, the length of this loop is in terms of length of plumbing in meters
  
#matrix for output from second part of model 
# number of rows = number of meters
# number of columns = 3: column 1 = distance from water heater; column 2 = temperature of water in pipe at that distance 
# from water heater; column 3 = concentration of Legionella 
matrix<-matrix(nrow=loop.1,ncol=3)

for(i in 2:loop.1){
  
  C.0.1<-matrix.wh[loops.wh,3] 
  #defining intial concentration of Legionella in premise plumbing to be
  #equal to final concentration in
  #water heater after holding time
  
  T.air<-21 #(deg C)
  T.w<-Temperature.waterheater #(deg C)
  c<-4.1868 #specific heat of water kJ/(kg x K)
  vel<-0.16 #velocity (m/s)
  r<-.0127/2 #radius (m)
  UA<-.346 #UA flow: W/(m x K)
  p<-1000 #density of water 
  proportion.infectious<-10^-4 #Acknowledging that only a proportion of Legionella may be infectious (coming from amoeba)
  
  #Decay rate
  r.d<-runif(1,2.8*10^-6,8.3*10^-6) #min and max values found in literature - using a uniform dist
  
  #initial concentration
  matrix[1,3]<-C.0.1 # assigning first concentration in loop to be equal to initial concentration in premise plumbing
  
  #distance from water heater
  matrix[i,1]<-i # distance from water heater in meters = # in the loop (starting at 2m since the concentration at 1st
  # meter is equal to final concentration measured in water heater)
  
  #temperature at that distance (deg C)
  matrix[i,2]<-T.air+(T.w-T.air)*exp(-(UA*i)/(c*vel*pi*r^2*p*1000))
  
  #change in concentration based on growth rate (temperature dependent)
  #Couldn't find a growth curve - only found ranges. So, based on the range of temperatures, a different
  #growth rate was assigned
  if(matrix[i,2]>=48){
    g.rate<-0 #/s
  }else if (matrix[i,2]<48 & matrix[i,2]>=40){
    g.rate<-0.162 #/s
  }else if (matrix[i,2]<40 & matrix[i,2]>=35){
    g.rate<-.048 #/s
  }else if (matrix[i,2]<35 & matrix[i,2]>=25){
    g.rate<-0.00139 #/s
  }else{
    g.rate<-0 #/s
  }
  #concentration of Legionella per meter away from the water heater
  matrix[i,3]<-matrix[i-1,3]+(matrix[i-1,3]*proportion.infectious*exp(g.rate*(1/vel)))-(matrix[i-1,3]*r.d*(1/vel))
  #concentration at current time step is equal to concentration at previous time step (matrix[i-1,3) plus
  #any growth that occurs ((matrix[i-1,3]*proportion.infectious*exp(g.rate*(1/vel)))) minus any decay that occurs
  #(matrix[i-1,3]*r.d*(1/vel))
}

#Reminder - output for premise plumbing (part 2) of model: 
# 1st column in matrix -> distance from water heater (m) at that moment in the loop
# 2nd column in the matrix -> temperature of water (deg C) at that moment in the loop 
# 3rd column in the matrix -> concentration of legionella (cfu/mL) in the pipe at that moment in the loop 
print("Part 2 = Check")



##############################################################


#######  PART 3: STAGNATION


##############################################################

#loop.2<-(24-duration.shower)*60*60 #assuming that the time of stagnation in the pipe is equal to 24 hours minus showering event duration

loop.2<-(24-15/60)*60*60  #Assume stagnation time is equal to 24 hours minus showering time event duration
# 15/60 = showering duration in hours (15 minutes or 1/4 hour)
# multiply (24-15/60) by 3600 to convert to seconds, since we will be tracking concentration change per second

matrix.2<<-matrix(nrow=loop.2,ncol=2)
#number of rows = number of seconds of assumed stagnation
#number of columns = 2: 1st column: temperature in pipe per second; 2nd column: concentration of Legionella

matrix.2[1,2]<-matrix[l.plumbing,3] #Initial concentration for this portion of the model
#This was the concentration at 433 m in the premise plumbing from the water heater

matrix.2[1,1]<-matrix[l.plumbing,2] #Temperature at end of premise-plumbing pipe is the intial temperature in the
#stagnation scenario

UA.2<-0.222 #not sure what this parameter is, but Emily might know -Amanda 11/14

Q.shower<-((0.425*pi*r^2)*1000*1000) #flow rate of shower (mL/s)

l.shower<-(Q.shower/(1000*1000))*15/(pi*r^2) #m  #not sure what this parameter is, but Emily might know -Amanda 11/14

#Determining cp for no flow:
water <- 4.1868*1000*pi*r^2 #not sure what this parameter is, but Emily might know -Amanda 11/14
pipe<-0.39*8960*((r)^2- ((r-2*.00124))^2) #not sure what this parameter is, but Emily might know -Amanda 11/14
insulation<-1.3*15*pi*(((r+.013)^2)-((r/2)^2)) #not sure what this parameter is, but Emily might know -Amanda 11/14

cp.noflow<-sum(water,pipe,insulation)*1000 #not sure what this parameter is, but Emily might know -Amanda 11/14


for(i in 2:loop.2){
  
  r.d<-runif(1,2.8*10^-6,8.3*10^-6) #inactivation distribution - assumed uniform distribution using min and max values from literature
  
  #Temperature
  #matrix.2[i,1]<-T.air+(matrix.2[i-1,1]-T.air)*exp(-(UA.2*i)/(cp.noflow))
  matrix.2[i,1]<-T.air+(matrix[l.plumbing,2] -T.air)*exp(-(UA.2*i)/(cp.noflow))
  #change in temperature per second
  
  #Defining growth rate: Growth rates associated with ranges of temperature
  if(matrix.2[i,1]>=48){
    g.rate<-0 #/s
  }else if (matrix.2[i,1]<48 & matrix.2[i,1]>=40){
    g.rate<-0.162 #/s
  }else if (matrix.2[i,1]<40 & matrix.2[i,1]>=35){
    g.rate<-.048 #/s
  }else if (matrix.2[i,1]<35 & matrix.2[i,1]>=25){
    g.rate<-0.00139 #/s
  }else{
    g.rate<-0 #/s
  }
  #Concentration of Legionella (CFU/ml)
  s<-1 #change in time is 1 second per step
  matrix.2[i,2]<-matrix.2[i-1,2]+(matrix.2[i-1,2]*proportion.infectious*exp(g.rate)*s)-(matrix.2[i-1,2]*r.d*s)
  
}


print ("Part 3 = Check")
#Reminder - output for part 3 stagnation pipe portion of model:
# 1st column in matrix.2 -> temperature (deg C) of water in pipe at that moment in the loop
# 2nd column in matrix.2 -> concentration of Legionella (deg C) in that pipe at moment in loop


##############################################################


#######  PART 4: SHOWERING EVENT


##############################################################


#Time length of shower 
loop.3<-duration.shower*60 #duration shower is given in minutes in fuction as form of convenience; this is converted to seconds in model so we can loop by seconds

#Setting up Matrix for Temp of water, Concentration of Legionella coming out of shower,
#concentration of Legionella in the room, and CFU inhaled
matrix.3<-matrix(nrow=loop.3,ncol=3)
# number of rows: duration of shower * 60 = time in seconds
#number of columns: 3: column 1: concentration of Legionella coming out of shower per second; 
#column 2: concentration of legionella in the room; column 3: time in shower

#Initial Conditions
initial.temp<-matrix.2[loop.2,1] #Initial Temperature of Water (deg C)
matrix.3[1,2]<-matrix.2[loop.2,2] #Initial Concentration of Legionella coming out of shower (cfu/mL)
matrix.3[1,3]<-0 #Initial concentration of Legionella in the room, respirable (assumption)

air.part<-10^-5 #Perkins et al., 2011
volume.shower<-4.06 #m^3 (based on ADA length and width; floor-to-floor height based on simple online calculator)

#Loop for Monte Carlo Simulation, Calculating Dose, and Calculating Associated Probabilities
#of Infection

for (i in 2:loop.3){
  
  ###############################################
  # Temperature, Concentration
  ###############################################
  
  #Time
  matrix.3[i,1]<-i
  
  
  
  
  
  
  #Concentration of Legionella in cold water line
  conc.cold<-2 #cfu/mL
  
  #Determining Amount of Mixing
  T.m<-desired.temp #Desired temperature for shower
  T.cold<-T.air #Assuming cold air that may be mixed will have same temperature as ambient temperature 
  T.h<-matrix[l.plumbing,2]

  P<-(T.m-T.cold)/(T.h-T.cold) #determines whether mixing is needed or not

  if(P>=1){
    P.effect<-1 #If P >= 1, we assume that all water is coming from hot water pipe 
    #Concentration coming out of shower head
    if (i<=(l.shower*pi*r^2)/(Q.shower*P.effect*(1/1000)*(1/1000))){
      matrix.3[i,2]<-matrix.3[1,2] # (cfu/ml) (we are assuming we have already accounted for proportion that may be infectious in previous steps)
      
    }else{ 
      matrix.3[i,2]<-matrix[l.plumbing,3] # (cfu/ml) (we are assuming we have already accounted for proportion that may be infectious in previous steps)
    }
  }else{ #if P< 1, we assume that the % of inflowing water coming from cold water pipe = 1 - P
    P.effect<-P
    #Concentration coming out of shower head
    if (i<=(l.shower*pi*r^2)/(Q.shower*P.effect*(1/1000)*(1/1000))){
      matrix.3[i,2]<-(matrix.3[1,2]*P) + (conc.cold*(P-1))  # (cfu/mL) (we are assuming we have already accounted for proportion that may be infectious in previous steps)
      #assuming concentration coming out of shower head is equal to concentration in stagnant pipe * the portion of water coming from
      #stagnant pipe (matrix.3[1,1]*P) plus the concentration in municipal water being mixed with hot water * the portion of water
      #coming from municipal water ((conc.cold*(P-1))) Amount of mixing is a function of what the water tempermature is set at
      
      }else{
      matrix.3[i,2]<-(matrix[l.plumbing,3]*P) + (conc.cold*(P-1)) # (cfu/mL) (we are assuming we have already accounted for proportion that may be infectious in previous steps)
      #after the volume of water in the stagnant pipe has been "used up," hot water will be coming directly from the premise plumbing
      #and so the concentration per second will now be equal to the final concentration of legionella in the premise plumbing
      #as it arrives at the shower * the portion of hot water being used (matrix[l.plumbing,3]*P) plus the concentration in municipal
      #water * the portion of municipal water being used (conc.cold*(P-1))
      }
  }
  
  #The if else statement above is communicating that the first portion of water
  #coming out of the stagnat pipe will contribute a different concentration of Legionella.
  #Once that water is "used up," water will be pulled from the piping system, where the
  #concentration of legionella will be assumed to be that calculatd for that distance
  #(433 m) from the water heater with growth/decay in the circulating pipe. If p < 1, a particular volume of cold water
  #may be mixed, diluting the concentration of legionella coming from the shower head. This is accounted for. 
  
  
  #Respirable Concentration in the Room
  
  matrix.3[i,3]<-matrix.3[i,2]*air.part*1000 #CFU/m^3
  
  #Q is in units of mL/second
  #air.part is the proportion of water coming out of the shower that is aerosolized
  #to a size of water droplet small enough to be inhaled
  #volume of shower should be in same units as inhalation for next step
  
  
  #matrix.3 1st column -> temperature of water (deg C) ###### STILL NEEDS TO BE INCORPORATED IN CODE ########
  #matrix.3 2nd column -> concentration of legionella in water coming from shower (cfu/mL)
  #matrix.3 3rd column -> concentration in room that is respirable (cfu/m^3)
  
}

print ("Part 4 = Check")
# * Note about Inhalation Rates - Source of Inhalation Rates: Exposure Factors Handbook 2011
#Inhalation rates are in units of m^3 / day. We are converting these values to m^3 / minute and
#then multiplying by the number of minutes given for "duration.shower"


#Distributions are on in table on page 6-24 of the handbook. Due to age ranges of adults being broken down into smaller
#sections in handbook than what was desired in model, parameters for adult, male and female, distributions were created by
#taking 100 random samples from each adult distribution (23 - <30, 30 - <40, 40 - <65)  (assumed normal) and fitting a normal distribution to those 300 values. 
#FitdistRplus was used to fit these values to a normal distribution. The code for this process is not shown here.

#Fraction factor 
#F1 - values from S.D. Perkins, R. Byerns, K. Cowen, D. Lorch, M. Shaw, M. Taylor, T. Nichols/Assessing exposures to shower aerosolized Brevundimonas diminuta and Pseudomonas aeruginosa/Applied and Environment Microbiology in Review (2011)
#F2 - values from R.B. Schlesinger/R.O. McClellan, R.F. Henderson (Eds.), Concepts in inhalation toxicology, Hemisphere Publishing Corp, New York (1989), pp. 163-192
F1.5<-.75 #Fraction of total aerosolized organisms in aerosols of size 1-5 um
F2.5<-.2 #Fraction of aerosols of size range 1-5 um deposited in alveoli
F1.6<-.09 #Fraction of total aerosolized organisms in aerosols of size 5-6 um
F2.6<-.1 #Fraction of aerosols of size range 5-6 um deposited in alveoli
F1.10<-.14#Fraction of total aerosolized organisms in aerosols of size 6-10 um
F2.10<-.01 #Fraction of aerosols of size range 6-10 um deposited in alveoli

fraction.factor<-(F1.5*F2.5)+(F1.6*F2.6)+(F1.10*F2.10)
#########################################################################

#Calculating Dose and Prob. of Infection for 1,000 people - Child - Male

#########################################################################

#Creating vectors to store cumulative dose (Dose.child.male) and probability of infection (Prob.infect.child.male) and dose per second (dose.moment)
Dose.child.male<<-numeric(1000)
Prob.infect.child.male<<-numeric(1000)
Dose.moment<<-numeric(1000) #Used for each loop for children, adults, and elderly - male and female

#Calculating Inhaled CFU per second per person.. then taking cumulative dose for
#whole showering event

#Exponential Dose-Response Model Recommended by CAMRA wiki
#k value for exponential dose-response:
k<-5.99*10^(-2) #50% percentile for k value


for(j in 1:1000){
  #inhalation rate
  inhalation.c.m<-rnorm(1,17.23,3.67)*(1/24)*(1/60) #Children - Male (11-<23 years old) - Exposure Factors Handbook
  for (i in 1:loop.3){
    #Dose per second
    Dose.moment[i]<-inhalation.c.m*(matrix.3[i,3])*fraction.factor
    # dose per second = inhalation rate (m^3/sec) * concentration of Legionella (CFU/mL) * fraction inhaled to alveolar region
  }
  Dose.child.male[j]<-sum(Dose.moment)
  #summed dose over full time of shower (loop.3) 

  Prob.infect.child.male[j]<-1-exp(-k*Dose.child.male[j])
  #Exponential dose-response
}
#hist(Prob.infect.child.male)

#########################################################################

#Calculating Dose and Prob. of Infection for 1,000 people - Child - Female

#########################################################################

Dose.child.female<<-numeric(1000)
Prob.infect.child.female<<-numeric(1000)

#Calculating Inhaled CFU per second per person.. then taking cumulative dose for
#whole showering event

for(j in 1:1000){
  
  inhalation.c.f<-rnorm(1,13.28,2.6)*(1/24)*(1/60) #Children - Female (11-<23 years old)
  
  for (i in 1:loop.3){
    Dose.moment[i]<-inhalation.c.f*matrix.3[i,3]*fraction.factor  # converting mL to L for concentration 
    
  }
  Dose.child.female[j]<-sum(Dose.moment)
  Prob.infect.child.female[j]<-1-exp(-k*Dose.child.female[j])
}


#########################################################################

#Calculating Dose and Prob. of Infection for 1,000 people - Adult - Male

#########################################################################

Dose.adult.male<<-numeric(1000)
Prob.infect.adult.male<<-numeric(1000)

#Calculating Inhaled CFU per second per person.. then taking cumulative dose for
#whole showering event

for(j in 1:1000){
  inhalation.a.m<-rnorm(1,16.94,2.615)*(1/24)*(1/60)        #Adults - Male (23-64 years)
  for (i in 1:loop.3){
    Dose.moment[i]<-inhalation.a.m*matrix.3[i,3]*fraction.factor
    
  }
  Dose.adult.male[j]<-sum(Dose.moment)
  Prob.infect.adult.male[j]<-1-exp(-k*Dose.adult.male[j])
}

#########################################################################

#Calculating Dose and Prob. of Infection for 1,000 people - Adult - Female

#########################################################################


Dose.adult.female<<-numeric(1000)
Prob.infect.adult.female<<-numeric(1000)

#Calculating Inhaled CFU per second per person.. then taking cumulative dose for
#whole showering event

for(j in 1:1000){
  
  inhalation.a.f<-rnorm(1,13.14,2.04)*(1/24)*(1/60)   #Adults - Female (23-64 years)
  for (i in 1:loop.3){
    Dose.moment[i]<-inhalation.a.f*matrix.3[i,3]*fraction.factor
  }
  Dose.adult.female[j]<-sum(Dose.moment)
  Prob.infect.adult.female[j]<-1-exp(-k*Dose.adult.female[j])
}

#########################################################################

#Calculating Dose and Prob. of Infection for 1,000 people - Elderly - Male

#########################################################################

Dose.elderly.male<<-numeric(1000)
Prob.infect.elderly.male<<-numeric(1000)

#Calculating Inhaled CFU per second per person.. then taking cumulative dose for
#whole showering event

for(j in 1:1000){
  
  inhalation.e.m<-rnorm(1,12.96,2.48)*(1/24)*(1/60) #Elderly - Male (65+ years)
  
  for (i in 1:loop.3){
    Dose.moment[i]<-inhalation.e.m*matrix.3[i,3]*fraction.factor 
    
  }
  Dose.elderly.male[j]<-sum(Dose.moment)
  Prob.infect.elderly.male[j]<-1-exp(-k*Dose.elderly.male[j])
}

#########################################################################

#Calculating Dose and Prob. of Infection for 1,000 people - Elderly - Female

#########################################################################

Dose.elderly.female<-numeric(1000)
Prob.infect.elderly.female<-numeric(1000)

#Calculating Inhaled CFU per second per person.. then taking cumulative dose for
#whole showering event

for(j in 1:1000){
  
  inhalation.e.f<-rnorm(1,9.80,2.17)*(1/24)*(1/60) #Elderly - Female (65+ years)
  for (i in 1:loop.3){
    Dose.moment[i]<-inhalation.e.f*matrix.3[i,3]*fraction.factor  
    #print(Dose.moment[i])
    #print(inhalation.e.f)
    #print(matrix.3[i,3])
  }
  Dose.elderly.female[j]<-sum(Dose.moment)
  Prob.infect.elderly.female[i]<-1-exp(-k*Dose.elderly.female[j])
}

#############################################################
##Calculating the illness incidence based on the infection rates
##Calculating the probability of illness given infection based on demographic information for age and sex
#prob_Disease(given infection) =demographic rates of incidence/(exposure_factor)???(probability of illness given exposure) )  ??????pop of demographic/???Pop???_USA 	

###Probability of Illness given infection Male Children
Dem.ill.child.male <- 0.002583
Exp.fact <- 1
Pop.child.male <-53606769
Pop.US <- 308745538
Prob.ill.child.male <- ((Dem.ill.child.male/Exp.fact)*Prob.infect.child.male)*(Pop.child.male/Pop.US)


###Probability of Illness given infection Female Children
Dem.ill.child.female <- 0.001476
Exp.fact <- 1
Pop.child.female <-51246786
Pop.US <- 308745538
Prob.ill.child.female <- ((Dem.ill.child.female/Exp.fact)*Prob.infect.child.female)*(Pop.child.female/Pop.US)

###Probability of Illness given infection Male Adults
Dem.ill.adult.male <- 0.010332
Exp.fact <- 1
Pop.adult.male <-80811597
Pop.US <- 308745538
Prob.ill.adult.male <- ((Dem.ill.adult.male/Exp.fact)*Prob.infect.adult.male)*(Pop.adult.male/Pop.US)

###Probability of Illness given infection Female Adults
Dem.ill.adult.female <- 0.005904
Exp.fact <- 1
Pop.adult.female <-82812402
Pop.US <- 308745538
Prob.ill.adult.female <- ((Dem.ill.adult.female/Exp.fact)*Prob.infect.adult.female)*(Pop.adult.female/Pop.US)

###Probability of Illness given infection Male Elderly
Dem.ill.elderly.male <- 0.108486
Exp.fact <- 1
Pop.elderly.male <-17362960
Pop.US <- 308745538
Prob.ill.elderly.male <- ((Dem.ill.elderly.male/Exp.fact)*Prob.infect.elderly.male)*(Pop.elderly.male/Pop.US)

###Probability of Illness given infection Female Elderly
Dem.ill.elderly.female <- 0.061992
Exp.fact <- 1
Pop.elderly.female <-22905024
Pop.US <- 308745538
Prob.ill.elderly.female <- ((Dem.ill.elderly.female/Exp.fact)*Prob.infect.elderly.female)*(Pop.elderly.female/Pop.US)



##############################################################


#######  PART 4: FINAL SUMMARY STATISTICS AND FIGURES

# summary(Prob.infect.child.male)
# summary(Prob.ill.child.male)
# summary(Prob.infect.child.female)
# summary(Prob.ill.child.female)
# summary(Prob.infect.adult.male)
# summary(Prob.ill.adult.male)
# summary(Prob.infect.adult.female)
# summary(Prob.ill.adult.female)
# summary(Prob.infect.elderly.male)
# summary(Prob.ill.elderly.male)
# summary(Prob.infect.elderly.female)
# summary(Prob.ill.elderly.female)

#plot(density(Prob.ill.adult.male),main = "Risk to Adults", xlab = "Risk of Illness", col="blue",ylim=c(0,2*10^13))
#lines(density(Prob.ill.adult.female),col="red")


#plot(density(Prob.ill.child.male),main = "Risk to Children", xlab = "Risk of Illness", col="blue",ylim=c(0,2*10^13.5))
#lines(density(Prob.ill.child.female),col="red")


#plot(density(Prob.ill.elderly.male),main = "Risk to the Elderly", xlab = "Risk of Illness", col="blue",ylim=c(0,2*10^13))
#lines(density(Prob.ill.elderly.female),col="red")

##############################################################
#Probability of illness for elderly when altering starting temperature (1.n) or pipe length (n.1)
##############################################################
#Starting temperature = 49C and pipe length varies

##############################################################
##############################################################
if(a==1 & b==1){
  Prob.ill.elderly.male.1.1 <- Prob.ill.elderly.male
  Prob.ill.elderly.female.1.1 <- Prob.ill.elderly.female
  Dose.ill.elderly.male.1.1<-Dose.elderly.male
  Dose.ill.elderly.female.1.1<-Dose.elderly.female}
if(a==1 & b==2){
  Prob.ill.elderly.male.1.2<-Prob.ill.elderly.male
  Prob.ill.elderly.female.1.2<-Prob.ill.elderly.female
  Dose.ill.elderly.male.1.2<-Dose.elderly.male
  Dose.ill.elderly.female.1.2<-Dose.elderly.female}
if(a==1 & b==3){
  Prob.ill.elderly.male.1.3<-Prob.ill.elderly.male
  Prob.ill.elderly.female.1.3<-Prob.ill.elderly.female
  Dose.ill.elderly.male.1.3<-Dose.elderly.male
  Dose.ill.elderly.female.1.3<-Dose.elderly.female}


#concentration.1<-matrix.3[,3]
#dose.m.1<-Dose.elderly.male
#dose.f.1<-Dose.elderly.female


if(a==2 & b==1){
  Prob.ill.elderly.male.2.1<-Prob.ill.elderly.male
  Prob.ill.elderly.female.2.1<-Prob.ill.elderly.female
  Dose.ill.elderly.male.2.1<-Dose.elderly.male
  Dose.ill.elderly.female.2.1<-Dose.elderly.female}
if(a==2 & b==2){
  Prob.ill.elderly.male.2.2<-Prob.ill.elderly.male
  Prob.ill.elderly.female.2.2<-Prob.ill.elderly.female
  Dose.ill.elderly.male.2.2<-Dose.elderly.male
  Dose.ill.elderly.female.2.2<-Dose.elderly.female}
if(a==2 & b==3){
  Prob.ill.elderly.male.2.3<-Prob.ill.elderly.male
  Prob.ill.elderly.female.2.3<-Prob.ill.elderly.female
  Dose.ill.elderly.male.2.3<-Dose.elderly.male
  Dose.ill.elderly.female.2.3<-Dose.elderly.female}

if(a==3 & b==1){
  Prob.ill.elderly.male.3.1<-Prob.ill.elderly.male
  Prob.ill.elderly.female.3.1<-Prob.ill.elderly.female
  Dose.ill.elderly.male.3.1<-Dose.elderly.male
  Dose.ill.elderly.female.3.1<-Dose.elderly.female}
if(a==3 & b==2){
  Prob.ill.elderly.male.3.2<-Prob.ill.elderly.male
  Prob.ill.elderly.female.3.2<-Prob.ill.elderly.female
  Dose.ill.elderly.male.3.2<-Dose.elderly.male
  Dose.ill.elderly.female.3.2<-Dose.elderly.female}
if(a==3 & b==3){
  Prob.ill.elderly.male.3.3<-Prob.ill.elderly.male
  Prob.ill.elderly.female.3.3<-Prob.ill.elderly.female
  Dose.ill.elderly.male.3.3<-Dose.elderly.male
  Dose.ill.elderly.female.3.3<-Dose.elderly.female}

##############################################################
###Probability for illness in Adults under varying pipe length and temperature
#Temperature 49C piple length varies

if(a==1 & b==1){
  Prob.ill.adult.male.1.1 <- Prob.ill.adult.male
  Prob.ill.adult.female.1.1 <- Prob.ill.adult.female
  Dose.ill.adult.male.1.1<-Dose.adult.male
  Dose.ill.adult.female.1.1<-Dose.adult.female}
if(a==1 & b==2){
  Prob.ill.adult.male.1.2<-Prob.ill.adult.male
  Prob.ill.adult.female.1.2<-Prob.ill.adult.female
  Dose.ill.adult.male.1.2<-Dose.adult.male
  Dose.ill.adult.female.1.2<-Dose.adult.female}
if(a==1 & b==3){
  Prob.ill.adult.male.1.3<-Prob.ill.adult.male
  Prob.ill.adult.female.1.3<-Prob.ill.adult.female
  Dose.ill.adult.male.1.3<-Dose.adult.male
  Dose.ill.adult.female.1.3<-Dose.adult.female}

#Starting temperature = 54C and pipe length varies

if(a==2 & b==1){
  Prob.ill.adult.male.2.1 <- Prob.ill.adult.male
  Prob.ill.adult.female.2.1 <- Prob.ill.adult.female
  Dose.ill.adult.male.2.1<-Dose.adult.male
  Dose.ill.adult.female.2.1<-Dose.adult.female}
if(a==2 & b==2){
  Prob.ill.adult.male.2.2<-Prob.ill.adult.male
  Prob.ill.adult.female.2.2<-Prob.ill.adult.female
  Dose.ill.adult.male.2.2<-Dose.adult.male
  Dose.ill.adult.female.2.2<-Dose.adult.female}
if(a==2 & b==3){
  Prob.ill.adult.male.2.3<-Prob.ill.adult.male
  Prob.ill.adult.female.2.3<-Prob.ill.adult.female
  Dose.ill.adult.male.2.3<-Dose.adult.male
  Dose.ill.adult.female.2.3<-Dose.adult.female}

#Starting temperature = 59C and pipe length varies

if(a==3 & b==1){
  Prob.ill.adult.male.3.1 <- Prob.ill.adult.male
  Prob.ill.adult.female.3.1 <- Prob.ill.adult.female
  Dose.ill.adult.male.3.1<-Dose.adult.male
  Dose.ill.adult.female.3.1<-Dose.adult.female}
if(a==3 & b==2){
  Prob.ill.adult.male.3.2<-Prob.ill.adult.male
  Prob.ill.adult.female.3.2<-Prob.ill.adult.female
  Dose.ill.adult.male.3.2<-Dose.adult.male
  Dose.ill.adult.female.3.2<-Dose.adult.female}
if(a==3 & b==3){
  Prob.ill.adult.male.3.3<-Prob.ill.adult.male
  Prob.ill.adult.female.3.3<-Prob.ill.adult.female
  Dose.ill.adult.male.3.3<-Dose.adult.male
  Dose.ill.adult.female.3.3<-Dose.adult.female}

##############################################################
###Probability for illness in Children under varying pipe length and temperature
#Temperature 49C piple length varies

if(a==1 & b==1){
  Prob.ill.child.male.1.1 <- Prob.ill.child.male
  Prob.ill.child.female.1.1 <- Prob.ill.child.female
  Dose.ill.child.male.1.1<-Dose.child.male
  Dose.ill.child.female.1.1<-Dose.child.female}
if(a==1 & b==2){
  Prob.ill.child.male.1.2<-Prob.ill.child.male
  Prob.ill.child.female.1.2<-Prob.ill.child.female
  Dose.ill.child.male.1.2<-Dose.child.male
  Dose.ill.child.female.1.2<-Dose.child.female}
if(a==1 & b==3){
  Prob.ill.child.male.1.3<-Prob.ill.child.male
  Prob.ill.child.female.1.3<-Prob.ill.child.female
  Dose.ill.child.male.1.3<-Dose.child.male
  Dose.ill.child.female.1.3<-Dose.child.female}

#Starting temperature 54C, pipe length varies

if(a==2 & b==1){
  Prob.ill.child.male.2.1<-Prob.ill.child.male
  Prob.ill.child.female.2.1<-Prob.ill.child.female
  Dose.ill.child.male.2.1<-Dose.child.male
  Dose.ill.child.female.2.1<-Dose.child.female}
if(a==2 & b==2){
  Prob.ill.child.male.2.2<-Prob.ill.child.male
  Prob.ill.child.female.2.2<-Prob.ill.child.female
  Dose.ill.child.male.2.2<-Dose.child.male
  Dose.ill.child.female.2.2<-Dose.child.female}
if(a==2 & b==3){
  Prob.ill.child.male.2.3<-Prob.ill.child.male
  Prob.ill.child.female.2.3<-Prob.ill.child.female
  Dose.ill.child.male.2.3<-Dose.child.male
  Dose.ill.child.female.2.3<-Dose.child.female}

#Starting temperature 59C, pipe length varies
if(a==3 & b==1){
  Prob.ill.child.male.3.1<-Prob.ill.child.male
  Prob.ill.child.female.3.1<-Prob.ill.child.female
  Dose.ill.child.male.3.1<-Dose.child.male
  Dose.ill.child.female.3.1<-Dose.child.female}
if(a==3 & b==2){
  Prob.ill.child.male.3.2<-Prob.ill.child.male
  Prob.ill.child.female.3.2<-Prob.ill.child.female
  Dose.ill.child.male.3.2<-Dose.child.male
  Dose.ill.child.female.3.2<-Dose.child.female}
if(a==3 & b==3){
  Prob.ill.child.male.3.3<-Prob.ill.child.male
  Prob.ill.child.female.3.3<-Prob.ill.child.female
  Dose.ill.child.male.3.3<-Dose.child.male
  Dose.ill.child.female.3.3<-Dose.child.female}
#############Emily
count <- ifelse(a==1 & b==1, 1, ifelse(a==1 & b==2, 2, ifelse(a==1 & b==3, 3,ifelse(a==2 & b==1, 4,
                                                                                    ifelse(a==2 & b==2, 5,ifelse(a==2 & b==3, 6, ifelse(a==3 & b==1, 7,ifelse(a==3 & b==2, 8,
                                                                                                                                                              ifelse(a==3 & b==3, 9,0)))))))))

finalcolnames <-  c("H1L1", "H1L2", "H1L3", "H2L1", "H2L2", "H2L3", "H3L1", "H3L2", "H3L3")
if(a==1 & b==1){
  rctime <- as.data.frame(matrix(NA, ncol=9, nrow=max(lengths)))
  colnames(rctime) <- finalcolnames
  rctemp <- as.data.frame(matrix(NA, ncol=9, nrow=max(lengths)))
  colnames(rctemp) <- finalcolnames
  rcconc <- as.data.frame(matrix(NA, ncol=9, nrow=max(lengths)))
  colnames(rcconc) <- finalcolnames
  
  sptime <- as.data.frame(matrix(NA, ncol=9, nrow=nrow(matrix.2)))
  colnames(sptime) <-finalcolnames
  sptemp <- as.data.frame(matrix(NA, ncol=9, nrow=nrow(matrix.2)))
  colnames(sptemp) <- finalcolnames
  spconc <- as.data.frame(matrix(NA, ncol=9, nrow=nrow(matrix.2)))
  colnames(spconc) <- finalcolnames
  
  showertime <- as.data.frame(matrix(NA, ncol=9, nrow=nrow(matrix.3)))
  colnames(showertime) <-finalcolnames
  showertemp <- as.data.frame(matrix(NA, ncol=9, nrow=nrow(matrix.3)))
  colnames(showertemp) <- finalcolnames
  showerconc <- as.data.frame(matrix(NA, ncol=9, nrow=nrow(matrix.3)))
  colnames(showerconc) <- finalcolnames
}

rctime[,count] <- c(matrix[,1], rep(NA, nrow(rctime)-nrow(matrix)))
rctemp[,count] <- c(matrix[,2], rep(NA, nrow(rctemp)-nrow(matrix)))
rcconc[,count] <- c(matrix[,3], rep(NA, nrow(rcconc)-nrow(matrix)))
sptemp[,count] <- matrix.2[,1]
spconc[,count] <- matrix.2[,2]

showertime[,count] <- matrix.3[,1]
showertemp[,count] <- matrix.3[,2]
showerconc[,count] <- matrix.3[,3]
    }}
################################################# RUN FROM HERE UP













sptime[,1:9] <- seq(1, nrow(sptime))

#Plot of concentration along pipe (matrix)
pdf("RC-Plots.pdf", width=8, height=8)
par(mfrow=c(3, 3), mar=c(2, 3, 1, 3), oma=c(4, 3.5, 2, 3))
for(j in 1:9){
  plot(rctime[,j], rctemp[,j], xlab="Distance along pipe (m)", ylab="Temperature (C)", las=1, type="l", col="red")
  if(j==1){mtext("49C:", 2, las=2, oma=TRUE, line=3.5, at=max(rctemp[which(!is.na(rctemp[,j])),j]))}
  if(j==4){mtext("54C:", 2, las=2, oma=TRUE, line=3.5, at=max(rctemp[which(!is.na(rctemp[,j])),j]))}
  if(j==7){mtext("59C:", 2, las=2, oma=TRUE, line=3.5, at=max(rctemp[which(!is.na(rctemp[,j])),j]))}
  if(j<4){mtext(c("50 m", "200 m", "400 m")[j], 3, oma=TRUE, line=1)}
  par(new=T)
  plot(rctime[,j], rcconc[,j], las=2, ylab="Concentration (CFU/mL)", xlab="", axes=FALSE, type="l")
  axis(4, las=2) #,  labels=round(rcconc[,j], 1))
}
mtext("Temperature (C)", 2, outer=TRUE, line=2, cex=1.2)
mtext("Distance along pipe (x)", 1, outer=TRUE, line=1.2, cex=1.2)
mtext("Concentration (CFU/mL)", 4, outer=TRUE, line=1.2, cex=1.2)
legend(200, min(rcconc[,j]-.1), c("Temperature", "Concentration"), col=c("red", "black"), lty=1, xpd=NA, bty='n')
dev.off()


#Plot of concentration in stagnant pipe
pdf("SP-Plots.pdf", width=8, height=8)
par(mfrow=c(3, 3), mar=c(2, 3, 1, 3), oma=c(4, 3.5, 2, 3))
for(j in 1:9){
  plot(sptime[,j]/60/60, sptemp[,j], xlab="Time (h)", ylab="Temperature (C)", las=1, type="l", col="red")
  if(j==1){mtext("49C:", 2, las=2, oma=TRUE, line=3.5, at=max(rctemp[which(!is.na(rctemp[,j])),j]))}
  if(j==4){mtext("54C:", 2, las=2, oma=TRUE, line=3.5, at=max(rctemp[which(!is.na(rctemp[,j])),j]))}
  if(j==7){mtext("59C:", 2, las=2, oma=TRUE, line=3.5, at=max(rctemp[which(!is.na(rctemp[,j])),j]))}
  if(j<4){mtext(c("50 m", "200 m", "400 m")[j], 3, oma=TRUE, line=1)}
  par(new=T)
  plot(sptime[,j]/60/60, log10(spconc[,j]), las=2, ylab="Concentration (CFU/mL)", xlab="", axes=FALSE, type="l")
  axis(4, las=2) #,  labels=round(rcconc[,j], 1))
}
mtext("Temperature (C)", 2, outer=TRUE, line=2, cex=1.2)
mtext("Time (s)", 1, outer=TRUE, line=1.2, cex=1.2)
mtext("Concentration (CFU/mL)", 4, outer=TRUE, line=1.2, cex=1.2)
legend(200, 7.9, c("Temperature", "Concentration"), col=c("red", "black"), lty=1, xpd=NA, bty='n')
dev.off()












ymax <- 5*10^-3
#pdf("Results.pdf", width=8, height=5) 
par(mfrow=c(3,3), mar=1*c(1, 2, 1, 3), oma=c(1, 3, 1, 1))
boxplot(Prob.ill.elderly.male.1.1,ylab="Probability of Illness",xlim=c(.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.1.1,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.1.2,ylab="Probability of Illness",xlim=c(.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.1.2,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.1.3,ylab="Probability of Illness",xlim=c(.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.1.3,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.2.1,ylab="Probability of Illness",xlim=c(0.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.2.1,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.2.2,ylab="Probability of Illness",xlim=c(0.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.2.2,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.2.3,ylab="Probability of Illness",xlim=c(0.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.2.3,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.3.1,ylab="Probability of Illness",xlim=c(0.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.3.1,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.3.2,ylab="Probability of Illness",xlim=c(0.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.3.2,add=TRUE,at=2, axes=FALSE)

boxplot(Prob.ill.elderly.male.3.3,ylab="Probability of Illness",xlim=c(0.5, 2.5), ylim=c(0,ymax), cex.axis=1.2, las=2)
boxplot(Prob.ill.elderly.female.3.3,add=TRUE,at=2, axes=FALSE)
dev.off()



# 
# 
# boxplot(Prob.ill.elderly.male.1,ylab="Probability of Illness",xlim=c(1,6.5),ylim=c(0,0.7*10^-4))
# boxplot(Prob.ill.elderly.female.1,add=TRUE,at=2)
# boxplot(Prob.ill.elderly.male.2,add=TRUE,at=3)
# boxplot(Prob.ill.elderly.female.2,add=TRUE,at=4)
# boxplot(Prob.ill.elderly.male.3,add=TRUE,at=5)
# boxplot(Prob.ill.elderly.female.3,add=TRUE,at=6)
# 
# 








##############################################################


#######  PART 5: SENSITIVITY ANALYSIS


##############################################################

####Spearman's Correlation for association between exponential of beginning L. pneumophila concentration and Probability of Illness###
expC.0 <- exp(matrix.wh[1:1000,3]) #Concentration of Legionella
cor.test(expC.0, Prob.ill.adult.female.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.female.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.adult.male.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.male.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.child.female.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.male.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.1.1,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.2.1,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.3.1,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.1.2,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.2.2,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.3.2,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.1.3,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.2.3,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)

cor.test(expC.0, Prob.ill.elderly.female.3.3,
         method = c("pearson"),
         exact = NULL, conf.level = 0.95)
###Correlation test between inhalation rates and probability of illness

cor.test(Dose.ill.adult.male.1.1, Prob.ill.adult.male.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.1.2, Prob.ill.adult.male.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.1.3, Prob.ill.adult.male.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.2.1, Prob.ill.adult.male.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.2.2, Prob.ill.adult.male.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.2.3, Prob.ill.adult.male.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.3.1, Prob.ill.adult.male.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.3.2, Prob.ill.adult.male.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.male.3.3, Prob.ill.adult.male.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.1.1, Prob.ill.adult.female.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.1.2, Prob.ill.adult.female.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.1.3, Prob.ill.adult.female.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.2.1, Prob.ill.adult.female.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.2.2, Prob.ill.adult.female.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.2.3, Prob.ill.adult.female.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.3.1, Prob.ill.adult.female.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.3.2, Prob.ill.adult.female.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.adult.female.3.3, Prob.ill.adult.female.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.1.1, Prob.ill.child.male.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.1.2, Prob.ill.child.male.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.1.3, Prob.ill.child.male.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.2.1, Prob.ill.child.male.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.2.2, Prob.ill.child.male.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.2.3, Prob.ill.child.male.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.3.1, Prob.ill.child.male.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.3.2, Prob.ill.child.male.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.male.3.3, Prob.ill.child.male.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.1.1, Prob.ill.child.female.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.1.2, Prob.ill.child.female.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.1.3, Prob.ill.child.female.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.2.1, Prob.ill.child.female.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.2.2, Prob.ill.child.female.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.2.3, Prob.ill.child.female.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.3.1, Prob.ill.child.female.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.3.2, Prob.ill.child.female.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.child.female.3.3, Prob.ill.child.female.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.1.1, Prob.ill.elderly.male.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.1.2, Prob.ill.elderly.male.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.1.3, Prob.ill.elderly.male.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.2.1, Prob.ill.elderly.male.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.2.2, Prob.ill.elderly.male.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.2.3, Prob.ill.elderly.male.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.3.1, Prob.ill.elderly.male.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.3.2, Prob.ill.elderly.male.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.male.3.3, Prob.ill.elderly.male.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.1.1, Prob.ill.elderly.female.1.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.1.2, Prob.ill.elderly.female.1.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.1.3, Prob.ill.elderly.female.1.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.2.1, Prob.ill.elderly.female.2.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.2.2, Prob.ill.elderly.female.2.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.2.3, Prob.ill.elderly.female.2.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.3.1, Prob.ill.elderly.female.3.1,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.3.2, Prob.ill.elderly.female.3.2,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)

cor.test(Dose.ill.elderly.female.3.3, Prob.ill.elderly.female.3.3,
         method = c("spearman"),
         exact = NULL, conf.level = 0.95)