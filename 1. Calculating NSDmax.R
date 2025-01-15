library(dplyr)
library(terra)


setwd("C:/Users/wewang/Desktop/Feedback")


#### This function calculates the rate of spread (m/min) for each fire day
ros<-function(A1,A0=NULL,hours=NULL,area.unit="m"){  
  if (is.null(hours)){hours<-rep(4,length(A1))}   
  if (is.null(A0)&length(A1)>1) {A0<-rep(0,length(A1));for (i in 1:(length(A1)-1)){A0[i+1]<-A0[i]+A1[i]}} 
  if (is.null(A0)&length(A1)==1) {A0<-rep(0,length(A1))}
  if (area.unit=="ha") {A1<-10000*A1;A0<-10000*A0}
  ros <-(2*((sqrt(A1+A0)-sqrt(A0))/1.772454))/(hours*60)
  ros
}


################## MAIN #######################
nbac <- terra::vect("Data/Fire_Growth/NBAC_Downloaded20240702/nbac_1972_2023_20240530_shp/nbac_1972_2023_20240530.shp")
Ecozones <- terra::vect("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones.shp")
Ecozones <- terra::project(Ecozones, nbac)
years <- seq(2001, 2023, 1)
ROS_outputs <- data.frame()
for (yi in years) {
  print(paste("Year:",yi,"/ Whole range:", min(years), "-", max(years)))
  growth_yi <- read.csv(paste0("Data/Fire_Growth/Fire_Growth_2001_2023/fire_growth_",yi,"_region.csv"))
  growth_yi <- subset(growth_yi, ADJ_HA>=50)
  growth_yi$FireID <- sapply(strsplit(growth_yi$FIRE_ID_YR, "_"), function(x) as.numeric(x[1]))
  growth_yi$Year <- sapply(strsplit(growth_yi$FIRE_ID_YR, "_"), function(x) as.numeric(x[2]))
  
  #### add ecozone based on NBAC locations
  fires <- unique(growth_yi$FireID)
  nbac_fires <- nbac[(nbac$YEAR==yi) & (nbac$NFIREID %in% fires), ]
  nbac_fires <- terra::aggregate(nbac_fires, by="NFIREID")
  zone_fires <- terra::relate(terra::centroids(nbac_fires), Ecozones, relation="within", pairs=TRUE)
  matid <- match(growth_yi$FireID, nbac_fires$NFIREID[zone_fires[,1]])
  growth_yi$ECOZONE <- Ecozones$ECOZONE[zone_fires[matid,2]]
  growth_yi$ZONE_NAME <- Ecozones$ZONE_NAME[zone_fires[matid,2]]
  growth_yi$ZONE_NOM <- Ecozones$ZONE_NOM[zone_fires[matid,2]]
  growth_yi <- na.omit(growth_yi)
  
  ### add ROS and SD
  growth_yi[, c("ROS", "SD")] <- NA
  for (fi in unique(growth_yi$FIRE_ID_YR)) {
    fi_rows <- which(growth_yi$FIRE_ID_YR==fi)
    growth_fi <- growth_yi[fi_rows,]
    growth_fi <- growth_fi[order(growth_fi$loc_jd), ]
    growth_yi$ROS[fi_rows] <- ros(growth_fi$aob, area.unit="ha")
  }
  growth_yi$SD <- ifelse(growth_yi$ROS>=1, 1, 0)
  ROS_outputs <- rbind(ROS_outputs, growth_yi)
}
write.table(ROS_outputs, file= paste0("AllFiresLT50ha_", min(years),"-",max(years),"_ROS_SD.csv"), sep=",", append=T, col.names=T, row.names=F)

### Summarize number of spread days (NSD) for each fire event
NSD_outputs <- ROS_outputs %>% group_by(Year, FireID, ECOZONE, ZONE_NAME, ZONE_NOM) %>% summarise(NSD=sum(SD), Fire_aob=sum(aob), Fire_AJDha=unique(ADJ_HA))
write.table(NSD_outputs, file= paste0("AllFiresLT50ha_", min(years),"-",max(years),"_NSD.csv"), sep=",", append=T, col.names=T, row.names=F)

### Summarize maximum NSD among all fire events for each year in each ecozone
NSDmax_outputs <- NSD_outputs %>% group_by(Year, ECOZONE, ZONE_NAME, ZONE_NOM) %>% summarise(NSD_max=max(NSD), AAB_aob=sum(Fire_aob), 
                                                                                             AAB_AJD=sum(Fire_AJDha), MFS_AJD=max(Fire_AJDha), ANF_AJD=n())
write.table(NSDmax_outputs, file= paste0("AllFiresLT50ha_", min(years),"-",max(years),"_NSDmax_AAB-MFS-ANF.csv"), sep=",", append=T, col.names=T, row.names=F)
