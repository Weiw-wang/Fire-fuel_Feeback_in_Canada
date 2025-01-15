library(dplyr)
library(tidyverse)
library(terra)

setwd("C:/Users/wewang/Desktop/Feedback")


#### This function generates the potential spread-event day distribution table
psed<-function(dat,n=10000,wt=c(0,0,0,0.020,0.016,0.50,0.25,0.20,0.02,0,0,0)){
  foo<-data.frame(PSED=1:200) 
  PSED<-PSED0<-NULL
  for (k in 1:n){
    koo0<-subset(dat,id==sample(unique(dat$id),1))
    koo1<-subset(koo0,yr==sample(unique(koo0$yr),1))
    koo2<-subset(koo1,mon==sample(1:12,1,prob=wt[1:12]))
    if (nrow(koo2)==0){PSED0<-0}else
      if(length(unique(koo2$spl))==1&unique(koo2$spl)[1]==0){PSED0<-0} else {
        n0<-which(row.names(koo1)==sample(row.names(koo2[koo2$spl!=0,]),1))
        foo0<-koo1[n0:nrow(koo1),]
        foo1<-subset(foo0,spl==0|spl==2)
        x<-rle(foo1$spl)
        x0<-data.frame(sed=x$values,freqs=x$lengths)
        PSED0<-ifelse(nrow(x0)>0&x0[1,1]==2,x0[1,2],0)}
    PSED<-c(PSED,PSED0)
  }
  PSED<-PSED[PSED>0]
  a<-tapply(PSED,PSED,length)
  a<-(a/sum(a))*100
  a<-data.frame(PSED=names(a),percent=as.numeric(a))
  names(a)[2]<-"pct"
  foo<-merge(foo,a,by.x="PSED",by.y="PSED",all.x=TRUE,all.y=FALSE)
  foo<-foo[!is.na(foo$pct),]
  foo
}


################## MAIN ############################
### Estimate fire counts for each month in each ecozone
fire_events <- read.csv("Results/1. Linking AAB and NSDmax/AllFiresLT50ha_2001-2023_ROS_SD.csv")
date_growth <- as.Date(fire_events$loc_jd - 1, origin = paste0(fire_events$Year, "-01-01"))
month_growth <- format(date_growth, "%m")
fire_events$Month <- as.integer(month_growth)
wt0<-tapply(fire_events$FireID,list(fire_events$ECOZONE,fire_events$Month),length)
wt1<-wt0/apply(wt0,1,function(.x)sum(na.omit(.x)))
wt1[is.na(wt1)]<-0
wt <- data.frame(matrix(0, nrow = nrow(wt1), ncol = 12))
rownames(wt) <- rownames(wt1); colnames(wt) <- 1:12
wt[rownames(wt1), colnames(wt1)] <- wt1

### period; climate changes
years <- seq(2001, 2023, 1)
GCMs <- c("CanESM2", "CSIRO-Mk3-6-0", "HadGEM2-ES")
RCPs <- c("rcp26", "rcp45", "rcp85")

### variables used to generate psd
wis <- c("dmc", "fwi")

### baseline weather
dat_base <- read.csv("Data/FWI_weather/1.PSD/observed_baseline/Baseline.csv")
colnames(dat_base) <- tolower(colnames(dat_base))
used_cols <- c("long", "lat", "yr", "mon", "day", wis)
dat_base <- dat_base[, used_cols]
dat_base$id <- paste(sprintf("%.3f", dat_base$long), sprintf("%.3f", dat_base$lat), sep="_")

## assign ecozone
Ecozones <- terra::vect("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones.shp")
Ecozones <- terra::project(Ecozones, "+proj=longlat +datum=WGS84")
ws_id <- unique(dat_base$id)
dat_points <- data.frame(id=ws_id, long=sapply(strsplit(ws_id, "_"), function(x) as.numeric(x[1])), 
                         lat=sapply(strsplit(ws_id, "_"), function(x) as.numeric(x[2])))
ws_points <- terra::vect(dat_points, geom=c("long", "lat"), crs="+proj=longlat +datum=WGS84")
ws_zones <- terra::relate(ws_points, Ecozones, relation="within", pairs=TRUE)

PSD<-data.frame()
for (gi in GCMs) {
  for (ri in RCPs) {
    print(paste(gi, ":", ri))
    dat_gri <- read.csv(paste0("Data/FWI_weather/3. FuturePSD&SD/fwi/",gi,"_2020s_",ri,".csv"))
    colnames(dat_gri) <- tolower(colnames(dat_gri))
    dat_gri <- dat_gri[, used_cols]
    dat_gri$id <- paste(sprintf("%.3f", dat_gri$long), sprintf("%.3f", dat_gri$lat), sep="_")
    dat_all <- rbind(dat_base, dat_gri)
    dat_all <- subset(dat_all, yr>=min(years) & yr<=max(years))
    
    ## assign ecozone
    matid <- match(dat_all$id, dat_points$id[ws_zones[,1]])
    dat_all$ECOZONE <- Ecozones$ECOZONE[ws_zones[matid,2]]
    dat_all$ZONE_NAME <- Ecozones$ZONE_NAME[ws_zones[matid,2]]
    dat_all$ZONE_NOM <- Ecozones$ZONE_NOM[ws_zones[matid,2]]
    dat_all <- na.omit(dat_all)
    
    for(zi in unique(dat_all$ECOZONE)) {
      dat_zonei <- subset(dat_all, ECOZONE==zi)
      wt_zi <- wt[which(row.names(wt)==zi),]
      dat_zonei$spl<-0
      dat_zonei$spl<-ifelse(dat_zonei$dmc>20,1,dat_zonei$spl)
      dat_zonei$spl<-ifelse(dat_zonei$dmc>20&dat_zonei$fwi>=19,2,dat_zonei$spl)
      
      ### psd distribution
      dat_zonei<-dat_zonei[with(dat_zonei,order(id,yr,mon,day)),]
      psed_out<-psed(dat_zonei,wt=wt_zi,n=50000)
      psed_out$ECOZONE <- zi
      psed_out$ZONE_NAME  <- unique(dat_zonei$ZONE_NAME)
      psed_out$ZONE_NOM <- unique(dat_zonei$ZONE_NOM)
      psed_out$GCM <- gi
      psed_out$RCP <- ri
      
      PSD<-rbind(PSD, psed_out)
    }
  }
}

write.csv(PSD, file=paste0("PSD_simulations_",min(years),"-",max(years),"_3GCMs_3RCPs.csv"),row.names=F) 

