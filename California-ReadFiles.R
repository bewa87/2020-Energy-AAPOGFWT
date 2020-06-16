# GOAL: WIND SPEED LOCATION DATA FOR FULL AVAILABLE TIME RANGE!

location          <- 69015093121

# DEFINE YOUR PATH OF YEARLY LOCATION DATA!

dir_name_help     <- "YOUR PATH"
dir_name          <- paste(dir_name_help,location,"/","Originaldaten/",sep="")

OFiles            <- list.files(path=dir_name, pattern="*", full.names=TRUE, recursive=FALSE)

original_table1   <- read.table(OFiles[1])
vel_help1         <- original_table1$V9
vel_help1         <- vel_help1[which(vel_help1>=0.0)]
vel_help1         <- vel_help1/10
vel_help1         <- t(vel_help1)
vel               <- vel_help1

for (j in 2:length(OFiles)){
  original_table2 <- read.table(OFiles[j])
  vel_help2       <- original_table2$V9
  vel_help2       <- vel_help2[which(vel_help2>=0.0)]
  vel_help2       <- vel_help2/10
  vel_help2       <- t(vel_help2)
  vel             <- cbind(vel,vel_help2)
}

name_velocity     <- paste(location,".","txt",sep="")
velocity          <- t(vel)
write.table(velocity, file=name_velocity, sep=";")