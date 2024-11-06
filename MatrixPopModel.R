##################################################################################################
########################## Coding Script #########################################################
############### by Aissa Thiombane and Laura Emrich ##############################################
##################################################################################################

####Loading the packages##########################################################################
library(popbio) #useful things for population biology
library(pracma) #practical maths
library(tidyverse) #tidy coding and plots

####Loading the Data##############################################################################
rawdata <- read.table("All_data.csv", header=TRUE, sep = ";", dec=",") #run data: name of file, including the first row with names, data is separated by semicolon, the values contain commas

afterwar <- rawdata %>%   #selecting and redefining the data we will use primarily 
  filter((Year >= 1945) & (Year < 2002)) %>% #only containing the years from 1945 and later, and without 2002
  view() #load table
  
####Defining the Values used in our Matrix########################################################
nbj <- 0.4 #successful annual breeds for first time breeders
nba <- 0.9 #successful annual breeds for adults
SR <- 0.5 #sex ratio
  
Sa <- afterwar$Adult_survival[1] #adult survival at [year 1 -> 1945]
Sj <- afterwar$Juvenile_survival[1] #juvenile survival at [year 1 -> 1945]
bs <- afterwar$Brood_size[2] #brood size for year at t+1 [year 2 -> 1946]

#The universal matrix we will be using from owl paper Altwegg et al. 2006, with 2 rows and 2 columns
M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) 

##################################################################################################
#### Model 1 (WITH Variation) ####################################################################
##################################################################################################

par(mfrow=c(1,2)) #display plots in x rows and y columns

#Loading the Data 
afterwar <- rawdata %>%  #selecting and redefining the data we will use primarily 
filter(Year >= 1945)  #only showing the years from 1945 and later

#Defining our Values
nbj <- 0.4 # % of successful annual breeds for first time breeders
nba <- 0.9 # % of successful annual breeds for adults
SR <- 0.5 # sex ratio
N <- c(118,82) # number of individuals at the beginning of time step (juveniles, adults)
Npop <- rep(0, length(afterwar$Year) - 1) # empty matrix for size of population each year -1
growth <- rep(0, length(afterwar$Year) - 1) # empty matrix for growth each year -1

#Creating the loop for every year (x57)
for(i in 1:(length(afterwar$Year) -1)){ #for every i in the amount of years -1 
  Sa <- afterwar$Adult_survival[i] #adult survival [year i]
  Sj <- afterwar$Juvenile_survival[i] #juvenile survival [year i]
  bs <- afterwar$Brood_size[i+1] #brood size for t+1
  M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) # Create a Matrix 
  Nprev <- N #save N before the projection
  N <-    M %*% N #the new population size is the matrix * the old N
  Npop[i] <- sum(N) #the population size is the sum of N at every point i (The sum of juveniles and adults)
  growth[i] <- sum(N)/sum(Nprev) #"lambda" for geometric growth rate
}

####Plotting the population size and the growth rate over 57 years
plot(Npop~c(1945:2001), type="l", #take the years from 1945-2001 from Npop, as a line
     main = "Model 1\nPopulation size",  #title for plot
     ylab = "Population size", xlab = "Year", cex.lab = 1.3, cex.axis = 1.3, #with axis labels, cex changing font size
     lwd = 2) # change thickness of line

plot(growth~c(1945:2001), type="l", #take the years from 1945-2001 from growth, as a line
     main = "Model 1\nGrowth", #add title
     ylab = "Growth rate", xlab = "Year", cex.lab = 1.3, cex.axis = 1.3, lwd = 2) #with axis labels, cex changing font size, thicker line width
abline(h=1, lty=2, col = "firebrick") #add a horizontal dashed line in red, to see turning point at 1

#Calculating the long term growth rate 
long_term_gr1 <- nthroot((Npop[57]/200), 57) #Calculating the long term growth rate ((Npop at end/population size in beginning), 57th root)


##################################################################################################
#### Model 2 (WITHOUT Variation) #################################################################
##################################################################################################

#Loading the Data
afterwar <- rawdata %>% #selecting and redefining the data we will use primarily
  filter((Year >= 1945)) #only showing the years from 1945 and later

#Defining our Values
nbj <- 0.4 # % of successful annual breeds for first time breeders
nba <- 0.9 # % of successful annual breeds for adults
SR <- 0.5 # sex ratio
N <- c(118,82) # number of individuals at the beginning of time step
Npop <- rep(0, length(afterwar$Year) - 1) # empty matrix for size of population each year
growth <- rep(0, length(afterwar$Year) - 1) # empty matrix for growth rate each year

#Creating the loop for every year, using the means of survival rates and brood size, so M is always the same
for(i in 1:(length(afterwar$Year) -1)){ #for every i in the amount of years -1 
  Sa <- mean(afterwar$Adult_survival, na.rm = TRUE) #constant survival through mean, no variation, stripping of NA values
  Sj <- mean(afterwar$Juvenile_survival, na.rm = TRUE) #constant survival through mean, no variation, stripping of NA values
  bs <- mean(afterwar$Brood_size) #constant brood size through mean, no variation
  M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) #Create a Matrix 
  Nprev <- N #save N before the projection
  N <-    M %*% N #the new population is the matrix-product of M * N
  Npop[i] <- sum(N) #the population size is the sum of N at every point i (The sum of juveniles and adults)
  growth[i] <- sum(N)/sum(Nprev) #"lambda" for geometric growth rate
}
Npop #print Npop

######Plotting the growth rate and the population size over 57 years
plot(Npop~c(1945:2001),type="l", main = "Model 2\nPopulation size", #plot population size from 1945 - 2001, as a line, add a title
     ylab = "Population size", xlab = "Year", cex.lab = 1.3, cex.axis = 1.3, lwd = 2) #name axes, change font size and line width

plot(growth~c(1945:2001), main = "Model 2\nGrowth", ylim = c(0.9, 1.2), #plot growth rate from 1945 - 2001, as a line, add a title, change y-axis range
     type="l", ylab = "Growth rate", xlab = "Year", cex.lab = 1.3, cex.axis = 1.3, lwd = 2) #name axes, change font size and line width
abline(h=1, lty=2,  col = "firebrick") #add a horizontal dashed line in red, to see turning point at 1

#Calculating the long term growth rate 
long_term_gr2 <- nthroot((Npop[57]/200), 57) #Calculating the long term growth rate ((Npop at end/population size in beginning), 57th root)


##################################################################################################
#### Model 3: Projection for 60 Years for 10k times###############################################
##################################################################################################

set.seed(26) #set a seed to get the same random number sampled. seed 26

#Loading the Data
afterwar <- rawdata %>% #selecting and redefining the data we will use primarily
  filter((Year >= 1945)) #only showing the years from 1945 and later

#Defining our Values
Nsimul <- 10000 #number of simulations, iterations 
Nend <- rep(0, Nsimul)  #empty matrix for population size at end of simulation
nbj <- 0.4 # % of successful annual breeds for first time breeders
nba <- 0.9 # % of successful annual breeds for adults
SR <- 0.5 #sex ratio
y <- 60 #number of years to project

for(j in 1:Nsimul){ #For every j, repeat simulation Nsimul times
 #giving parameters names for iteration j 
  Npop <- rep(0, length(afterwar$Year) - 1 + y) #empty matrix for size of population each year
  growth <- rep(0, length(afterwar$Year) - 1 + y) #empty matrix for growth each year
  N <- c(118,82) #number of individuals at the beginning of time step
  
  for(i in 1:(length(afterwar$Year) -1)){ #for every year -1
   #same as model 1 
    Sa <- afterwar$Adult_survival[i] #adult survival [year i]
    Sj <- afterwar$Juvenile_survival[i] #juvenile ""
    bs <- afterwar$Brood_size[i+1] #brood size for t+1
    M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) #create universal matrix
    Nprev <- N #save N before the projection
    N <-    M %*% N #the new population is the matrix * the old N
    Npop[i] <- sum(N) #the population size is the sum of N at every point i
    growth[i] <- sum(N)/sum(Nprev) #"lambda" for geometric growth rate
  } 
  #future projection
  for(i in length(afterwar$Year):(length(afterwar$Year)-1 +y)){ #for every future year from 58 to year 57 + y (here:Year 58 to 118)
    a <- sample(1:57, size = 1, replace = TRUE) #sample a random number between 1 and 57, 1 number of values taken out, replacing values after sampling
    Sa <- afterwar$Adult_survival[a] #take the adult survival from the random matrix number a
    Sj <- afterwar$Juvenile_survival[a] #take the juvenile survival from the random matrix number a
    bs <- afterwar$Brood_size[a+1] #take the brood size from the random matrix number a
    M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) #Create a Matrix for the random number/Year
    Nprev <- N #save N before the projection
    N <-  M %*% N #the new population is the matrix * the old N
    Npop[i] <- sum(N) #the population size is the sum of N at every point i
    growth[i] <- sum(N)/sum(Nprev) #"lambda" for geometric growth rate
    
  }
  realgrowthsim <- growth #save model 3 growth under this name for 
  Nend[j] <- Npop[length(afterwar$Year)-1 + y] #saving the last Npop in every simulation, for j simulations
  print(j) #show where you are in the simulation 
}

e3 <- Nend[Nend < 10] #when the end population is smaller than 10, we say the species go extinct
length(e3) #count the number of Extinctions

#Comparing the actual End-Pop-size to the one we calculated from the long term growth rate
exp(mean(log(Nend))) #Calculating the (geometric) mean end Population size from our 10k simulations
200*long_term_gr^118 #calculating the end Population size with the long term Growth rate (200 is the start Pop size at 1945 and 118 is the total number of years for which we calculated the Population size and the growth rate)
 
#Histogram of Nend for model 3
hist(log(Nend), col = "mediumorchid1", #frequency of Nend model 3, taking log to see values near 0, change colour
     main = "Final population size\n Model 3", #title
     xlab = "log(Final population size)", #x-axis label
     cex.lab = 1.3, cex.axis = 1.3, #font size
    ylim = c(0, 2500), xlim = c(-15, 20)) #manually adjust axes area
abline(v = mean(log(Nend)), col = "black", lwd = 2) #mark the mean with a vertical line, thicker
text(x = 11, y = 2500, "mean", col = "black") #label the mean at point x and y, name it
abline(v = log(10), lt = 2, col ="mediumblue", lwd = 2) #mark the extinction point at <10 individuals
text(x = -2.5, y = 2250, "extinct", col = "mediumblue") #label the mean at point x and y, name it

#calculate the mean of log(Nend) for model 3
meanhist3 <- mean(log(Nend))
meanhist3


##################################################################################################
#### Model 4 Projection for 60 Years for 10k times with increased Probability of harsh winters####
##################################################################################################

#Loading the Data
afterwar <- rawdata %>% #selecting and redefining the data we will use primarily
  filter((Year >= 1945)) #only showing the years from 1945 and later

#Defining our Values
Nsimul <- 10000 #number of simulations 
Nend <- rep(0, Nsimul)  #population size at end of simulation
nbj <- 0.4 # % of successful annual breeds for first time breeders
nba <- 0.9# % of successful annual breeds for adults
SR <- 0.5 #sex ratio
y <- 60 #number of years to project
w <- 2 #probability of sampling a bad year (here doubled)
p <- rep(1, (length(afterwar$Year) - 1)) #probability matrix with all 1 (100 %)
p[3] <- w #change probability of 3rd Year (1947) in p to w
p[8] <- w #change probability of 8th Year (1952) in p to w
p[18] <- w #change probability of 18th Year (1962) in p to w

for(j in 1:Nsimul){ #For every j, repeat simulation Nsimul times
  #giving parameters names for iteration j
  Npop <- rep(0, length(afterwar$Year) - 1 + y) #empty matrix for size of population each year
  growth <- rep(0, length(afterwar$Year) - 1 + y) #empty matrix for growth
  N <- c(118,82) #number of individuals at the beginning of time step
  
  for(i in 1:(length(afterwar$Year) -1)){ #for every year -1
  #same as model 3  
    Sa <- afterwar$Adult_survival[i] #adult survival [year i]
    Sj <- afterwar$Juvenile_survival[i] #juvenile " "
    bs <- afterwar$Brood_size[i+1] #brood size for year i+1
    M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) #create universal matrix
    Nprev <- N #save N before the projection
    N <-    M %*% N #the new population is the matrix * the old N
    Npop[i] <- sum(N) #the population size is the sum of N at every point i
    growth[i] <- sum(N)/sum(Nprev) #"lambda" for geometric growth rate
  } 
  #same as model 3, except adding probability in a
  for(i in length(afterwar$Year):(length(afterwar$Year)-1 + y)){ #for every future year (Year 57 to 118)
    a <- sample(1:57, size = 1, replace = TRUE, prob = p) #sample a random number between 1 and 57, number of values taken out, replacing values after sampling, with probabilities p
    Sa <- afterwar$Adult_survival[a] #take the adult survival from the random matrix number a
    Sj <- afterwar$Juvenile_survival[a] #take the juvenile survival from the random matrix number a
    bs <- afterwar$Brood_size[a+1] #take the brood size from the random matrix number a+1
    M <- matrix(c(bs * nbj * SR * Sj, Sj, bs * nba * SR * Sa, Sa), nrow = 2, ncol = 2) #Create a Matrix for the random number/Year
    Nprev <- N #save N before the projection
    N <-  M %*% N #the new population is the matrix * the old N
    Npop[i] <- sum(N) #the population size is the sum of N at every point i
    growth[i] <- sum(N)/sum(Nprev) #"lambda" for geometric growth rate
    
  }
  
  Nend[j] <- Npop[length(afterwar$Year)-1 +y] #saving the last Npop in simulation, for j simulations
  print(j) #show where you are in the simulation 
}

e4 <- Nend[Nend < 10] #when the end population is smaller than 10, we say the species goes extinct
length(e4) #count the number of Extinctions

#Histogram of Nend for model 4
hist(log(Nend), col = "mediumorchid1", #making a logarithmic scale of Nend, giving it a colour
     main = "Final population size\n Model 4", #title
     xlab = "log(Final population size)", #naming x-axis
     cex.lab = 1.3, cex.axis = 1.3, #change font size 
    ylim = c(0, 2500), xlim = c(-15, 20)) #manually adjust axes areas
abline(v = mean(log(Nend)), col = "black", lwd = 2) #mark the mean with a vertical line, thicker line
text(x = 7, y = 2500, "mean", col = "black") #label the mean at point x and y, name it
abline(v = log(10), lt = 2, col ="mediumblue", lwd = 2) #mark the extinction point at <10 individuals, thicker line
text(x = -2.5, y = 2250, "extinct", col = "mediumblue") #label the mean at point x and y, name it

#show means of histograms for Nend
meanhist4 <- mean(log(Nend)) 
meanhist4

#both growth plots together
plot(growth~c(1945:(2002-1 +y)), col="green3", type="l", #growth rate of model 4
     ylim=c(0,max(growth,realgrowthsim)), #define y-axis limit
     main = "Model 3 & 4\nCombined growth rate", ylab = "Growth rate", xlab = "Year", #adding title and axes-labels
     cex.lab = 1.3, cex.axis = 1.3, lwd = 2) #changing font size and line width
lines(realgrowthsim~c(1945:(2002-1 +y)), lwd = 2) #growthrate of model 3
abline(v=2002, col = "royalblue3") #add a vertical line at 2002, to see where future projection starts
abline(h=1, lty=2, col = "firebrick") #add dashed horizontal line at 1, to better see when lamda is below or above 1
mtext("Start future prediction", col = "royalblue", side = 3, cex = 0.75) #label of blue line at top

#Calculating the (geometric) mean end Population size from our 10k simulations
exp(mean(log(Nend))) 

#Extinctions count
length(e3) #show number of extinctions for model 3
length(e4) #show number of extinctions for model 4
