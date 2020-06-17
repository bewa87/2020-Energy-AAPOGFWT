#######################################################################
# Code written by Benjamin Wacker, 2020                               #
#                                                                     #
# Goal: Wind Speed Data Analysis by Assuming Two-Parameter            #
#       Weibull Distribution, Four-Parameter Kappa Distribution       #
#       and Five-Parameter Wakeby Distribution                        #
#######################################################################

library(lmomco)
library(fitdistrplus)

Aktuell <- Sys.time()

################################
# Step 1: Power Curve Modeling #
################################

# Step 1.1: Data Preparation for Power Curve Data

# DEFINE YOUR LOADING PATH FOR POWER CURVE DATA!

PCData     <- read.table("YOUR PATH/vestas_data_powercurve_R.txt", header=TRUE, sep=";")
vel_log    <- PCData[,1]
pow_log    <- PCData[,2]

for (j in seq(1,length(vel_log),1)){
  while (pow_log[j]==0.0){
    j <- j + 1
  }
  break
}
start_cubic <- j

for (j in seq(start_cubic,length(vel_log),1)){
  while (pow_log[j]!=pow_log[j+1]){
    j <- j+1
  }
  break
}
end_cubic <- j

for (j in seq(end_cubic,length(vel_log),1)){
  while (pow_log[j]!=0.0){
    j <- j+1
  }
  break
}
end_power <- j

vel_cubic    <- vel_log[start_cubic:end_cubic]
pow_cubic    <- pow_log[start_cubic:end_cubic]

# Step 1.2: Cubic Spline Interpolation

cubic_power     <- splinefun(vel_cubic,pow_cubic)

# Step 1.3: Final Power Curve

Powercurve <- function(x){
  (0.0 <= x & x < vel_log[start_cubic])*0.0 + (vel_log[start_cubic] <= x & x < vel_log[end_cubic])*cubic_power(x) + (vel_log[end_cubic] <= x & x < vel_log[end_power])*3075 + (x >= vel_log[end_power])*0
}

######################################################
# Step 2: Helper functions for Weibull Distributions #
######################################################

Weibullcurve         <- function(x,wei_A_fitdist,wei_k_fitdist){
  (wei_k_fitdist/wei_A_fitdist) * ((x/wei_A_fitdist)^(wei_k_fitdist-1)) * exp( -(x/wei_A_fitdist)^(wei_k_fitdist) )
}

Weibulldiff_v        <- function(x,wei_A_fitdist,wei_k_fitdist){
  ( ( wei_k_fitdist-1 ) * wei_k_fitdist * exp(-(x/wei_A_fitdist)^(wei_k_fitdist)) * (x/wei_A_fitdist)^(wei_k_fitdist-2) )/(wei_A_fitdist^2) - ( wei_k_fitdist^2 * exp(-(x/wei_A_fitdist)^(wei_k_fitdist)) * (x/wei_A_fitdist)^(2*wei_k_fitdist - 2) )/(wei_A_fitdist^2)
}

Weibulldiff_A        <- function(x,wei_A_fitdist,wei_k_fitdist){
  ( wei_k_fitdist^2 * exp(-(x/wei_A_fitdist)^(wei_k_fitdist)) * (x/wei_A_fitdist)^(wei_k_fitdist-1) * ( ( (x/wei_A_fitdist)^(wei_k_fitdist) )-1 ) )/(wei_A_fitdist^2)
}

Weibulldiff_k        <- function(x,wei_A_fitdist,wei_k_fitdist){
  ( exp(-(x/wei_A_fitdist)^(wei_k_fitdist)) * (x/wei_A_fitdist)^(wei_k_fitdist-1) )/(wei_A_fitdist) + ( wei_k_fitdist * exp(-(x/wei_A_fitdist)^(wei_k_fitdist)) * (x/wei_A_fitdist)^(wei_k_fitdist-1) * log(x/wei_A_fitdist) )/(wei_A_fitdist) - ( wei_k_fitdist * exp(-(x/wei_A_fitdist)^(wei_k_fitdist)) * (x/wei_A_fitdist)^(2*wei_k_fitdist-1) * log(x/wei_A_fitdist) )/(wei_A_fitdist)
}

#####################################
# Step 3: Loop for DWD Station Data #
#####################################

# DEFINE PATH OF SORTED FULL WIND SPEED STATION DATA AND ONE SORTED FILE OF STATION META DATA BY STATION ID!

OFiles     <- list.files(path="YOUR PATH OF WIND SPEED STATION DATA", pattern="*.txt", full.names=TRUE, recursive=FALSE)
SFiles     <- list.files(path="YOUR PATH OF STATION META DATA", pattern="*.txt", full.names=TRUE, recursive=FALSE)

# Build an Array for Data File Paths

length     <- length(OFiles)
Files      <- array(c(OFiles, SFiles), dim=c(length,2,1))

######################################################################
# Step 4: Create Matrix                                              #
######################################################################

# Pre-Step: Initialize Final Array for Statistical Analysis

GoalMatrix  <- array(0, dim=c(length+1,21,1))
GoalMatrix2 <- array(0, dim=c(length+1,11,1))

# Step 4.1: Create First Row with Headings

GoalMatrix[1,1,1]   <- "StationID"
GoalMatrix[1,2,1]   <- "Location"
GoalMatrix[1,3,1]   <- "Longitude"
GoalMatrix[1,4,1]   <- "Latitude"
GoalMatrix[1,5,1]   <- "StationHeight"
GoalMatrix[1,6,1]   <- "ReadingHeight"
GoalMatrix[1,7,1]   <- "NumberOfData"
GoalMatrix[1,8,1]   <- "MeanValue"
GoalMatrix[1,9,1]   <- "StandardDeviation"
GoalMatrix[1,10,1]  <- "MinVelocity"
GoalMatrix[1,11,1]  <- "MaxVelocity"
GoalMatrix[1,12,1]  <- "Weibull_k"
GoalMatrix[1,13,1]  <- "Weibull_A"
GoalMatrix[1,14,1]  <- "Semi-empirical_Pow"
GoalMatrix[1,15,1]  <- "Estimated_Pow_Weibull"
GoalMatrix[1,16,1]  <- "Error_Weibull"
GoalMatrix[1,17,1]  <- "Difference_Weibull"
GoalMatrix[1,18,1]  <- "Estimated_Pow_Kappa"
GoalMatrix[1,19,1]  <- "Difference_Kappa"
GoalMatrix[1,20,1]  <- "Estimated_Power_Wakeby"
GoalMatrix[1,21,1]  <- "Difference_Wakeby"

GoalMatrix2[1,1,1]  <- "StationID"
GoalMatrix2[1,2,1]  <- "Weibull_k"
GoalMatrix2[1,3,1]  <- "Weibull_A"
GoalMatrix2[1,4,1]  <- "Semi-empirical_Pow"
GoalMatrix2[1,5,1]  <- "Estimated_Pow_Weibull"
GoalMatrix2[1,6,1]  <- "Error_Weibull"
GoalMatrix2[1,7,1]  <- "Difference_Weibull"
GoalMatrix2[1,8,1]  <- "Estimated_Pow_Kappa"
GoalMatrix2[1,9,1]  <- "Difference_Kappa"
GoalMatrix2[1,10,1] <- "Estimated_Power_Wakeby"
GoalMatrix2[1,11,1] <- "Difference_Wakeby"

# Step 4.2: Loop to fill GoalMatrix

for(j in 1:length){
  
# Statistical Values
  
A                    <- read.table(Files[j,1,1], header=TRUE, sep=";")
B                    <- read.table(Files[j,2,1], header=TRUE, sep=";")
GoalMatrix[j+1,1,1]  <- B[1,1,1]                        #Station ID
GoalMatrix[j+1,2,1]  <- as.character(B[2,2,1])          #Name of Station
GoalMatrix[j+1,3,1]  <- B[dim(B)[1],3,1]                #Longitude
GoalMatrix[j+1,4,1]  <- B[dim(B)[1],4,1]                #Latitude
GoalMatrix[j+1,5,1]  <- B[dim(B)[1],5,1]                #Height of Actual Station (Last Entry)
GoalMatrix[j+1,6,1]  <- B[dim(B)[1],6,1]                #Reading Height of Actual Station (Last Entry)

# Preparation of Wind Speeds

C                    <- A[,4]                           #Wind Speeds
C                    <- C[which(C>=0.0)]
C_Weibull            <- C[which(C>=0.1)]                #Wind Speeds for Weibull
GoalMatrix[j+1,7,1]  <- length(C)
GoalMatrix[j+1,8,1]  <- mean(C)
GoalMatrix[j+1,9,1]  <- sd(C)
GoalMatrix[j+1,10,1] <- min(C)
GoalMatrix[j+1,11,1] <- max(C)

# Weibull Distribution

wei                  <- fitdist(C_Weibull, "weibull", method="mle")
wei_k_fitdist        <- getElement(wei$estimate, "shape")
wei_A_fitdist        <- getElement(wei$estimate, "scale")
wei_sigma_k_fitdist  <- getElement(wei$sd, "shape")
wei_sigma_A_fitdist  <- getElement(wei$sd, "scale")
GoalMatrix[j+1,12,1] <- wei_k_fitdist
GoalMatrix[j+1,13,1] <- wei_A_fitdist
GoalMatrix[j+1,14,1] <- sum(Powercurve(C))*24*365/(as.numeric(GoalMatrix[j+1,7,1])*1000000)

N_weibull              <- length(C_Weibull)
P_hourly_weibull       <- sum( 0.1 * Powercurve(seq(0,vel_log[end_power],0.1)) * Weibullcurve(seq(0,vel_log[end_power],0.1),wei_A_fitdist,wei_k_fitdist) )
Delta_P_hourly_weibull <- sum( Powercurve(seq(0.1,vel_log[end_power],0.1)) * Weibullcurve(seq(0.1,vel_log[end_power],0.1),wei_A_fitdist,wei_k_fitdist) * ( abs(Weibulldiff_v(seq(0.1,vel_log[end_power],0.1),wei_A_fitdist,wei_k_fitdist))*0.1 + abs(Weibulldiff_A(seq(0.1,vel_log[end_power],0.1),wei_A_fitdist,wei_k_fitdist))*abs(wei_sigma_A_fitdist) + abs(Weibulldiff_k(seq(0.1,vel_log[end_power],0.1),wei_A_fitdist,wei_k_fitdist))*abs(wei_sigma_k_fitdist) )) + sum( abs(cubic_power(seq(vel_log[start_cubic],vel_log[end_cubic]),deriv=1))*0.1*Weibullcurve(seq(vel_log[start_cubic],vel_log[end_cubic]),wei_A_fitdist,wei_k_fitdist) ) - sum( abs(cubic_power(seq(vel_log[start_cubic],vel_log[end_cubic]),deriv=1))*0.1*( abs(Weibulldiff_v(seq(vel_log[start_cubic],vel_log[end_cubic]),wei_A_fitdist,wei_k_fitdist))*0.1 + abs(Weibulldiff_A(seq(vel_log[start_cubic],vel_log[end_cubic]),wei_A_fitdist,wei_k_fitdist))*abs(wei_sigma_A_fitdist) + abs(Weibulldiff_k(seq(vel_log[start_cubic],vel_log[end_cubic]),wei_A_fitdist,wei_k_fitdist))*abs(wei_sigma_k_fitdist) ))
P_annual_weibull       <- 24*365*P_hourly_weibull/1000000
Delta_P_annual_weibull <- 24*365*Delta_P_hourly_weibull/1000000
GoalMatrix[j+1,15,1]   <- P_annual_weibull
GoalMatrix[j+1,16,1]   <- Delta_P_annual_weibull
GoalMatrix[j+1,17,1]   <- abs(as.numeric(GoalMatrix[j+1,14,1])-as.numeric(GoalMatrix[j+1,15,1]))

# Kappa Distribution and Wakeby Distribution

lmom_wak               <- lmoms(C,7)
par_kap                <- parkap(lmom_wak,checklmom=TRUE)
par_wak                <- parwak(lmom_wak,checklmom=TRUE)

N_kappa                <- length(C)
P_hourly_kappa         <- sum(0.1*Powercurve(seq(vel_log[start_cubic],vel_log[end_power],0.1))*pdfkap(seq(vel_log[start_cubic],vel_log[end_power],0.1),par_kap))
P_annual_kappa         <- 24*365*P_hourly_kappa/1000000
GoalMatrix[j+1,18,1]   <- P_annual_kappa
GoalMatrix[j+1,19,1]   <- abs(as.numeric(GoalMatrix[j+1,14,1])-as.numeric(GoalMatrix[j+1,18,1]))

N_wakeby               <- length(C)
P_hourly_wakeby        <- sum(0.1*Powercurve(seq(vel_log[start_cubic],vel_log[end_power],0.1))*pdfwak(seq(vel_log[start_cubic],vel_log[end_power],0.1),par_wak))
P_annual_wakeby        <- 24*365*P_hourly_wakeby/1000000
GoalMatrix[j+1,20,1]   <- P_annual_wakeby
GoalMatrix[j+1,21,1]   <- abs(as.numeric(GoalMatrix[j+1,14,1])-as.numeric(GoalMatrix[j+1,20,1]))

# Interesting Results in other Matrix

GoalMatrix2[j+1,1,1]  <- B[1,1,1]
GoalMatrix2[j+1,2,1]  <- GoalMatrix[j+1,12,1]
GoalMatrix2[j+1,3,1]  <- GoalMatrix[j+1,13,1]
GoalMatrix2[j+1,4,1]  <- GoalMatrix[j+1,14,1]
GoalMatrix2[j+1,5,1]  <- GoalMatrix[j+1,15,1]
GoalMatrix2[j+1,6,1]  <- GoalMatrix[j+1,16,1]
GoalMatrix2[j+1,7,1]  <- GoalMatrix[j+1,17,1]
GoalMatrix2[j+1,8,1]  <- GoalMatrix[j+1,18,1]
GoalMatrix2[j+1,9,1]  <- GoalMatrix[j+1,19,1]
GoalMatrix2[j+1,10,1] <- GoalMatrix[j+1,20,1]
GoalMatrix2[j+1,11,1] <- GoalMatrix[j+1,21,1]

# Data Cleaning for next loop

A                      <- NULL
B                      <- NULL
C                      <- NULL
C_Weibull              <- NULL
wei                    <- NULL
wei_k_fitdist          <- NULL
wei_A_fitdist          <- NULL
wei_sigma_k_fitdist    <- NULL
wei_sigma_A_fitdist    <- NULL
N_weibull              <- NULL
P_hourly_weibull       <- NULL
Delta_P_hourly_weibull <- NULL
P_annual_weibull       <- NULL
Delta_P_annual_weibull <- NULL
lmom_wak               <- NULL
par_kap                <- NULL
par_wak                <- NULL
N_kappa                <- NULL
P_hourly_kappa         <- NULL
P_annual_kappa         <- NULL
N_wakeby               <- NULL
P_hourly_wakeby        <- NULL
P_annual_wakeby        <- NULL
}

######################################################################
# Step 5: Save Results in Seperate File                              #
#         Results_01: Contains all Results                           #
#         Results_02: Contains solely Station Informations and       #
#                     Parameters/Results for Weibull, Results for    #
#                     Kappa and Wakeby Distributions                 #
######################################################################

write.table(GoalMatrix, file="Results_01.txt", sep=";")
write.table(GoalMatrix2, file="Results_02.txt", sep=";")
Ende <- Sys.time()