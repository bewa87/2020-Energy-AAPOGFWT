#######################################################################
# Code written by Benjamin Wacker, 2020                               #
#                                                                     #
# Goal: Wind Speed Data Analysis by Assuming Two-Parameter            #
#       Weibull Distribution and Power Curve Modelling with           #
#       Logistic Regression Function                                  #
#######################################################################

library(EnvStats)
library(neldermead)

Aktuell <- Sys.time()

#######################################################################
# Step 1: Get path files (Files) for Power Curve Data, Original Data  # 
#         and Station Data                                            #
#######################################################################

# Get Power Curve and prepare Logistic Regression (Vestas112-Wind Turbine) - domain for regression must be defined
# CHOOSE YOUR PATH FOR POWER CURVE DATA

PCData     <- read.table(".../Daten/vestas_data_powercurve_R.txt", header=TRUE, sep=";")
vel_log    <- PCData[,1]
pow_log    <- PCData[,2]
vel_log    <- vel_log[6:28]
pow_log    <- pow_log[6:28]

# Extract files in Original Data Path from produkt_ff_stunde_XXXXX and
# in Station Data Path from Metadaten_Geraete_Windgeschwindigkeit_XXXXX
# CHOOSE YOUR PATH FOR HOURLY WIND SPEED DATA AND CHOOSE YOUR PATH
# FOR META STATION DATA

OFiles     <- list.files(path=".../Daten/Original-Data", pattern="*.txt", full.names=TRUE, recursive=FALSE)
SFiles     <- list.files(path=".../Daten/Station-Data", pattern="*.txt", full.names=TRUE, recursive=FALSE)

# Build an Array for Data File Paths

length     <- length(OFiles)
Files      <- array(c(OFiles, SFiles), dim=c(length,2,1))

######################################################################
# Step 2: Create Matrix                                              #
######################################################################

# Pre-Step: Initialize Final Array for Statistical Analysis

GoalMatrix  <- array(0, dim=c(length+1,15,1))
GoalMatrix2 <- array(0, dim=c(length+1,5,1))

# Step 2.1: Create First Row with Headings

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
GoalMatrix[1,12,1]  <- "ShapePark_Env"
GoalMatrix[1,13,1]  <- "ScaleParA_Env"
GoalMatrix[1,14,1]  <- "Semi-empirical_Pow"
GoalMatrix[1,15,1]  <- "Estimated_Pow"

GoalMatrix2[1,1,1]  <- "StationID"
GoalMatrix2[1,2,1]  <- "ShapePark_Env"
GoalMatrix2[1,3,1]  <- "ShapeParA_Env"
GoalMatrix2[1,4,1]  <- "Semi-empirical_Pow"
GoalMatrix2[1,5,1]  <- "Estimated_Pow"

# Step 2.2: Logistic Regression for Power Curve

# Definition of Logistic Regression Function

CostLogReg <- function(x){
  y <- sum((x[1]/(x[2]+x[3]*exp(x[4]*vel_log+x[5]))-pow_log)^2)
}

# Solution of Least-Squares Problem

sol        <- fminsearch(CostLogReg, c(3.0,0.0,1.0,-1.0,1.0))
opt_par    <- neldermead.get(sol,"xopt")

# Final Logistic Regression Funktion

PowerCurve <- function(x){
  P <- (0 <= x & x <3)*0.0+(3 <= x & x < 13)*(opt_par[1]/(opt_par[2]+opt_par[3]*exp(opt_par[4]*x+opt_par[5])))+(13 <= x & x < 25.5)*3110+(x >= 25.5)*0
}

plot(PCData[,1],PowerCurve(PCData[,1]),xlim=c(0.0,40.0),ylim=c(0.0,3200.0))

# Step 2.3: Loop to fill GoalMatrix

for(j in 1:length){

# Statistical Values

A                    <- read.table(Files[j,1,1], header=TRUE, sep=";")
B                    <- read.table(Files[j,2,1], header=TRUE, sep=";")
GoalMatrix[j+1,1,1]  <- B[1,1,1]
GoalMatrix[j+1,2,1]  <- as.character(B[2,2,1])          #Name of Station
GoalMatrix[j+1,3,1]  <- B[dim(B)[1],3,1]                #Longitude
GoalMatrix[j+1,4,1]  <- B[dim(B)[1],4,1]                #Latitude
GoalMatrix[j+1,5,1]  <- B[dim(B)[1],5,1]                #Height of Actual Station (Last Entry)
GoalMatrix[j+1,6,1]  <- B[dim(B)[1],6,1]                #Reading Height of Actual Station (Last Entry)
C                    <- A[,4]
#C                    <- subset(A[,4],1999010100<=(A$MESS_DATUM) & (A$MESS_DATUM)<=2016123123)
C                    <- C[which(C>=0.1)]
GoalMatrix[j+1,7,1]  <- length(C)
GoalMatrix[j+1,8,1]  <- mean(C)
GoalMatrix[j+1,9,1]  <- sd(C)
GoalMatrix[j+1,10,1] <- min(C)
GoalMatrix[j+1,11,1] <- max(C)
C_new_80             <- C #C*(80.0/(as.numeric(GoalMatrix[j+1,5,1])+as.numeric(GoalMatrix[j+1,6,1])))^0.2 #C
wei_new_80_Env       <- eweibull(C_new_80, method="mle")
wei_k_Env            <- getElement(wei_new_80_Env$parameters, "shape")
wei_A_Env            <- getElement(wei_new_80_Env$parameters, "scale")
GoalMatrix[j+1,12,1] <- wei_k_Env
GoalMatrix[j+1,13,1] <- wei_A_Env
GoalMatrix[j+1,14,1] <- sum(PowerCurve(C))*24*365/(as.numeric(GoalMatrix[j+1,7,1])*1000000)
GoalMatrix[j+1,15,1] <- sum(PowerCurve(seq(0,25,by=0.1))*(wei_k_Env/wei_A_Env)*(seq(0,25,by=0.1)/wei_A_Env)^(wei_k_Env-1)*exp(-(seq(0,25,by=0.1)/wei_A_Env)^(wei_k_Env)))*(0.1*24*365/1000000)

GoalMatrix2[j+1,1,1] <- B[1,1,1]
GoalMatrix2[j+1,2,1] <- GoalMatrix[j+1,12,1]
GoalMatrix2[j+1,3,1] <- GoalMatrix[j+1,13,1]
GoalMatrix2[j+1,4,1] <- GoalMatrix[j+1,14,1]
GoalMatrix2[j+1,5,1] <- GoalMatrix[j+1,15,1]

A                    <- NULL
B                    <- NULL
C                    <- NULL
C_new_80             <- NULL
wei_new_80_Env       <- NULL
wei_k_Env            <- NULL
wei_A_Env            <- NULL
}

######################################################################
# Step 3: Save Results in Seperate File                              #
#         Results_01: Contains all Results                           #
#         Results_02: Contains solely Station Informations and       #
#                     Parameters for Weibull Distribution            #
######################################################################

write.table(GoalMatrix, file="Results_01.txt", sep=";")
write.table(GoalMatrix2, file="Results_02.txt", sep=";")
Ende <- Sys.time()
