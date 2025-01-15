library(dplyr)
library(tidyverse)


setwd("C:/Users/wewang/Desktop/Feedback")


#### observation data
dat_hist <- read.csv("Results/3. Historical AAB-MFS-ANF/AllFiresLT50ha_1959-2023_AAB-MFS-ANF_NFDBpoint-NBAC.csv")
dat0_2023 <- subset(dat_hist, Year==2023)
dat_2023 <- gather(dat0_2023, var_y, Observe2023, AAB:ANF)

##### projection data
dat_proj <- read.csv("Results/5. Predict AAB-MFS-ANF/maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections.csv")
comm_cols <- c("GCM","RCP","Period","Year","var_y","ECOZONE","ZONE_NAME","ZONE_NOM")
dat_sta <- dat_proj[dat_proj$Period!="Baseline" & dat_proj$Feedback=="sumY", c(comm_cols, "Static")]
names(dat_sta)[names(dat_sta)=="Static"] <- "Projection"
dat_sta$Model <- "Static"
dat_dyn <- dat_proj[dat_proj$Period!="Baseline", c(comm_cols,"Feedback","Dynamic")]
names(dat_dyn)[names(dat_dyn)=="Dynamic"] <- "Projection"
names(dat_dyn)[names(dat_dyn)=="Feedback"] <- "Model"
dat_proj_all <- rbind(dat_sta, dat_dyn)

#### output1: yearly projections & percentiles
dat_proj_pct <- dat_proj_all %>% group_by(RCP,Period,var_y,Model,ECOZONE,ZONE_NAME,ZONE_NOM) %>%
  summarise(p50=quantile(Projection, probs=0.50), p75=quantile(Projection, probs=0.75), p90=quantile(Projection, probs=0.90))
dat_proj_pct <- merge(dat_proj_pct, dat_2023[, c("ECOZONE","ZONE_NAME","ZONE_NOM","var_y","Observe2023")], by=c("ECOZONE","ZONE_NAME","ZONE_NOM","var_y"), all=TRUE)
write.csv(dat_proj_pct, file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections_quantiles3GCMs.csv", row.names=F)

#### output2: number of years exceeding the 2023 fire season
dat_proj_2023 <- merge(dat_proj_all, dat_2023[, c("ECOZONE","ZONE_NAME","ZONE_NOM","var_y","Observe2023")], by=c("ECOZONE","ZONE_NAME","ZONE_NOM","var_y"), all=TRUE)
dat_proj_2023$Exceed <- ifelse(dat_proj_2023$Projection > dat_proj_2023$Observe2023, 1, 0)
proj_2023 <- dat_proj_2023 %>% group_by(ECOZONE,ZONE_NAME,ZONE_NOM,var_y,GCM,RCP,Period,Model) %>% summarise(NyrsLT2023=sum(Exceed))

### + averaged over ecozones
proj_2023_na <- proj_2023 %>% group_by(var_y,GCM,RCP,Period,Model) %>% summarise(Mean_zone=mean(NyrsLT2023))
names(proj_2023_na)[names(proj_2023_na)=="Mean_zone"] <- "NyrsLT2023"
proj_2023_na$ECOZONE <- 100
proj_2023_na$ZONE_NAME <- "National"
proj_2023_na$ZONE_NOM <- "National"
proj_2023 <- rbind(proj_2023, proj_2023_na)
write.csv(proj_2023, file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections_NyrsLT2023_All.csv", row.names=F)

### averaged over 3GCMs
proj_2023_summ <- proj_2023 %>% group_by(ECOZONE,ZONE_NAME,ZONE_NOM,var_y,RCP,Period,Model) %>% 
  summarise(Mean_NyrsLT2023=mean(NyrsLT2023))
proj_2023_summ$GCM <- "Ave3GCMs"
proj_2023_summ_re <- spread(proj_2023_summ, key = RCP, value = Mean_NyrsLT2023)
write.csv(proj_2023_summ_re, file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections_NyrsLT2023_AveGCM.csv", row.names=F)


