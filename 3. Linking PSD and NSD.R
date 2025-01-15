library(Matching)
library(dplyr)

setwd("C:/Users/wewang/Desktop/Feedback")


##################### A. Fit NSD distribution ############################
years <- seq(2001, 2023, 1)
dat_NSD <- read.csv(paste0("Results/1. Linking AAB and NSDmax/AllFiresLT50ha_",min(years),"-",max(years),"_NSD.csv"), head=T)
rsd <- data.frame()
for (zi in unique(dat_NSD$ECOZONE)) {
  dat_NSD_zi <- subset(dat_NSD, ECOZONE==zi)
  ZONE_NAME_zi <- unique(dat_NSD_zi$ZONE_NAME)
  ZONE_NOM_zi <- unique(dat_NSD_zi$ZONE_NOM)
  
  x<-tapply(dat_NSD_zi$NSD,dat_NSD_zi$NSD,length)
  n<-as.numeric(x)
  y1<-scale(log(n))
  x1<-(as.numeric(names(x)))
  
  fm1<-lm(y1~x1)
  coef_a1 <- fm1$coefficients[[1]]
  coef_b1 <- fm1$coefficients[[2]]
  rsq_adj1 <- summary(fm1)$adj.r.squared
  rsq1 <- summary(fm1)$r.squared
  
  rsd <- rbind(rsd, data.frame(ECOZONE=zi, ZONE_NAME=ZONE_NAME_zi, ZONE_NOM=ZONE_NOM_zi, Intercept=coef_a1, Slope=coef_b1, rsq=rsq1, rsq_adj=rsq_adj1))
}


##################### B. Fit PSD distribution ############################
dat_PSD <- read.csv(paste0("Results/2.PSD simulations_Linking PSD and NSD/PSD_simulations_",min(years),"-",max(years),"_3GCMs_3RCPs.csv"))
fuu <- data.frame()
for (gi in unique(dat_PSD$GCM)) {
  for(ri in unique(dat_PSD$RCP)){
    dat_PSD_gri <- subset(dat_PSD, GCM==gi & RCP==ri)
    for(zi in unique(dat_PSD_gri$ECOZONE)){
      dat_PSD_zi <- subset(dat_PSD_gri, ECOZONE==zi)
      
      y2<-as.numeric(scale(log(dat_PSD_zi$pct)))
      x2<-(dat_PSD_zi$PSED)
      fm2<-lm(y2~x2)
      
      coef_a2 <- fm2$coefficients[[1]]
      coef_b2 <- fm2$coefficients[[2]]
      rsq_adj2 <- summary(fm2)$adj.r.squared
      rsq2 <- summary(fm2)$r.squared
      
      fuu <- rbind(fuu, data.frame(GCM=gi, RCP=ri, ECOZONE=zi, Intercept0=coef_a2, Slope0=coef_b2, rsq0=rsq2, rsq_adj0=rsq_adj2))
    }
  }
}
psd <- fuu %>% group_by(ECOZONE) %>% summarise(Intercept=mean(Intercept0), Slope=mean(Slope0), rsq=mean(rsq0), rsq_adj=mean(rsq_adj0))



############## C. Conversion funtions between PSD and NSD ################
matchid <- match(psd$ECOZONE, rsd$ECOZONE)
table3 <- subset(rsd, select=c(ECOZONE, ZONE_NAME, ZONE_NOM))
table3$a_psd[matchid] <- psd$Intercept
table3$b_psd[matchid] <- psd$Slope
table3$rsq_psd[matchid] <- psd$rsq
table3$rsq_adj_psd[matchid] <- psd$rsq_adj
table3$a_rsd <- rsd$Intercept
table3$b_rsd <- rsd$Slope
table3$rsq_rsd <- rsd$rsq
table3$rsq_adj_rsd <- rsd$rsq_adj
table3$a_prd <- (table3$a_psd-table3$a_rsd)/table3$b_rsd
table3$b_prd <- table3$b_psd/table3$b_rsd

write.csv(table3,paste0("NSD and PSD Linkage_functions(single log) ",min(years),"-",max(years),".csv"),row.names=F)

