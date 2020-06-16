require(RCurl)

location   <- 72381523161
start      <- 1973
ende       <- 2019
sixnumber  <- 723815
fivenumber <- 23161

url        <- "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/"

filenames  <- array(0,dim=c(ende-start+1,1))

for (j in start:ende){
  filenames[j-start+1,1] <- paste(url,j,"/",sixnumber,"-",fivenumber,"-",j,".","gz",sep="")
}

# DEFINE YOUR PATH WHERE YOU WANT TO SAVE THE YEARLY STATION FILES!

for (j in start:ende){
  download.file(filenames[j-start+1,1],paste("YOUR PATH",location,"/",sixnumber,"-",fivenumber,"-",j,".","gz",sep=""))
}