library(dplyr)
library(tidyverse)
library(terra)

setwd("C:/Users/wewang/Desktop/Feedback")


##### This function calculates maxPSD for all grids in one year
Cal_maxPSD <- function(dat){
  foo <- data.frame()
  for (zi in unique(dat$ECOZONE)){
    dat1<-subset(dat,ECOZONE==zi)
    ZONE_NAME_zi <- unique(dat1$ZONE_NAME)
    ZONE_NOM_zi <- unique(dat1$ZONE_NOM)
    dat1$spl<-0
    dat1$spl<-ifelse(dat1$dmc>20,1,dat1$spl)
    dat1$spl<-ifelse(dat1$dmc>20&dat1$fwi>=19,2,dat1$spl)
    foo1<-subset(dat1,spl==0|spl==2)
    for (ii in unique(foo1$id)) {
      foo2 <- subset(foo1, id==ii)
      foo2<-foo2[with(foo2,order(id,yr,mon,day)),]
      x<-rle(foo2$spl)
      x0<-data.frame(sed=x$values,freqs=x$lengths)
      x1<-subset(x0,sed==2)
      x1_maxPSD <- ifelse(nrow(x1)>0, max(x1$freqs), 0)
      foo <- rbind(foo, data.frame(ECOZONE = zi, ZONE_NAME=ZONE_NAME_zi, ZONE_NOM=ZONE_NOM_zi, 
                                   id=ii, PSDmax=x1_maxPSD))
    }
  }
  foo
}


####################### MAIN #########################
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

foo_out <- data.frame()
## baseline
matid_b <- match(dat_base$id, dat_points$id[ws_zones[,1]])
dat_base$ECOZONE <- Ecozones$ECOZONE[ws_zones[matid_b,2]]
dat_base$ZONE_NAME <- Ecozones$ZONE_NAME[ws_zones[matid_b,2]]
dat_base$ZONE_NOM <- Ecozones$ZONE_NOM[ws_zones[matid_b,2]]
dat_base <- na.omit(dat_base)
for (yi in unique(dat_base$yr)){
  datb_yi <- subset(dat_base, yr==yi)
  foob_yi <- Cal_maxPSD(datb_yi)
  foob_yi$GCM <- "Baseline"
  foob_yi$RCP <- "Baseline"
  foob_yi$Period <- "Baseline"
  foob_yi$yr <- yi
  foo_out <- rbind(foo_out, foob_yi)
}

## future
GCMs <- c("CanESM2", "CSIRO-Mk3-6-0", "HadGEM2-ES")
RCPs <- c("rcp26", "rcp45", "rcp85")
Periods <- c("2020s", "2050s", "2080s")
for (gi in GCMs) {
  for (ri in RCPs) {
    for (pi in Periods) {
      print(paste(gi, ":", ri, ":", pi))
      dat_grp <- read.csv(paste0("Data/FWI_weather/3. FuturePSD&SD/fwi/",gi,"_",pi,"_",ri,".csv"))
      colnames(dat_grp) <- tolower(colnames(dat_grp))
      dat_grp <- dat_grp[, used_cols]
      dat_grp$id <- paste(sprintf("%.3f", dat_grp$long), sprintf("%.3f", dat_grp$lat), sep="_")
      
      ## assign ecozone
      matid <- match(dat_grp$id, dat_points$id[ws_zones[,1]])
      dat_grp$ECOZONE <- Ecozones$ECOZONE[ws_zones[matid,2]]
      dat_grp$ZONE_NAME <- Ecozones$ZONE_NAME[ws_zones[matid,2]]
      dat_grp$ZONE_NOM <- Ecozones$ZONE_NOM[ws_zones[matid,2]]
      dat_grp <- na.omit(dat_grp)
      
      for (yi in unique(dat_grp$yr)){
        dat_yi <- subset(dat_grp, yr==yi)
        foo_yi <- Cal_maxPSD(dat_yi)
        foo_yi$GCM <- gi
        foo_yi$RCP <- ri
        foo_yi$Period <- pi
        foo_yi$yr <- yi
        foo_out <- rbind(foo_out, foo_yi)
      }
    }
  }
}
saveRDS(foo_out, "maxPSD_ALL_years_Grids.rds")
foo_out_summary <- foo_out %>% group_by(GCM, RCP, Period, yr, ECOZONE, ZONE_NAME, ZONE_NOM) %>% summarise(max_PSDmax = max(PSDmax), median_PSDmax = median(PSDmax), 
                                 mean_PSDmax = round(mean(PSDmax)),p90_PSDmax = as.numeric(quantile(PSDmax, 0.9)),
                                 p95_PSDmax = as.numeric(quantile(PSDmax, 0.95)), p99_PSDmax = as.numeric(quantile(PSDmax, 0.99)))
write.csv(foo_out_summary, file="maxPSD_ALL_years.csv",row.names=F)


