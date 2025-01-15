library(dplyr)
library(tidyverse)


setwd("C:/Users/wewang/Desktop/Feedback")


#### This function get the projections for one fire activity variable using static model & one form of dynamic model
proj_per <- function(dat_proj, dat_obs, fdbk_series, funcs_stadyn, vy, dyn){
  dat_proj <- dat_proj[order(dat_proj$Year), ]
  dat_obs <- dat_obs[order(dat_obs$Year), ]
  fdbk_series <- fdbk_series[order(fdbk_series$Year), ]
  
  ## observations
  dat_proj$var_y <- vy
  dat_proj$Observe <- dat_obs[match(dat_proj$Year, dat_obs$Year), vy]
  
  ## Static model
  func_sta <- subset(funcs_stadyn, Model=="Static" & var_y==vy)
  Yzi_sta <- exp(func_sta$coef_a + func_sta$coef_b1*log(dat_proj$NSDmax))
  Yzi_sta[Yzi_sta<0] <- 0 
  dat_proj$Static <- Yzi_sta
  
  ## dynamic model
  func_dyn <- subset(funcs_stadyn, Model=="Dynamic" & var_y==vy & Feedback==dyn)
  tao_zi <- func_dyn$Longevity
  dat_proj$Feedback <- dyn
  dat_proj$Dynamic <- NA
  dat_proj$Fterm <- NA
  for (yi in dat_proj$Year) {
    row_yi <- which(fdbk_series$Year==(yi-1))
    #### for different feedback forms
    if (dyn=="sumY") {
      fedb_yi <- sum(fdbk_series$AAB[(row_yi-tao_zi+1):(row_yi)], na.rm=TRUE)
    } else if (dyn=="sumY%i") {
      fedb_yi <- sum(fdbk_series$AAB[(row_yi-tao_zi+1):(row_yi)] / c(tao_zi:1), na.rm=TRUE)
    } else if (dyn=="sumY e^-i") {
      fedb_yi <- sum(fdbk_series$AAB[(row_yi-tao_zi+1):(row_yi)] * exp(- (c(tao_zi:1)-1) * exp(1) / tao_zi), na.rm=TRUE)
    } else if (dyn=="sumY (1-cos)%2") {
      fedb_yi <- sum(fdbk_series$AAB[(row_yi-tao_zi+1):(row_yi)] * ((1 + cos(3.141593 * c(tao_zi:1)/(tao_zi+1) ) )/2), na.rm=TRUE)
    }
    fedb_yi_log <- ifelse(fedb_yi!=0, log(fedb_yi), 0)
    Yzi_dyn <-  exp(func_dyn$coef_a + func_dyn$coef_b1*log(dat_proj$NSDmax[dat_proj$Year==yi]) + func_dyn$coef_b2*fedb_yi_log)
    Yzi_dyn[Yzi_dyn<0] <- 0
    
    dat_proj$Dynamic[dat_proj$Year==yi] <- Yzi_dyn
    dat_proj$Fterm[dat_proj$Year==yi] <- fedb_yi
    
    if (vy=="AAB") {
      fdbk_series$AAB[row_yi+1] <- Yzi_dyn
    }
  }
  
  return(dat_proj)
}


################## MAIN ################
### PSDmax data
foo_in <- read.csv("Results/2. PSD simulations_PSD-NSD_PSDmax/maxPSD_ALL_years.csv")
foo_in$PSDmax <- foo_in$max_PSDmax
foo_in <- subset(foo_in, select=c(GCM,RCP,Period,yr,ECOZONE,ZONE_NAME,ZONE_NOM,PSDmax))
names(foo_in)[names(foo_in)=="yr"] <- "Year"

### observation data for feedback term
dat_feedback <- read.csv("Results/3. Historical AAB-MFS-ANF/AllFiresLT50ha_1959-2023_AAB-MFS-ANF_NFDBpoint-NBAC-Growth.csv")

### functions
func_prd<-read.csv("Results/2. PSD simulations_PSD-NSD_PSDmax/NSD and PSD Linkage_functions(single log) 2001-2023.csv") 
funcs <- read.csv("Results/4. Static Dynamic Models Fits/AllFiresLT50ha_StaticDynamicModels_Fit2001-2023_SelectedTao.csv")

### parameters
GCMs <- c("CanESM2", "CSIRO-Mk3-6-0", "HadGEM2-ES")
RCPs <- c("rcp26", "rcp45", "rcp85")
years_all <- seq(1981,2100,1)
foo_out <- data.frame()
for (zi in unique(foo_in$ECOZONE)){
  print(zi)
  dat_feedback_zi <- subset(dat_feedback, ECOZONE==zi)
  dat_baseline_zi <- subset(foo_in, ECOZONE==zi & Period=="Baseline")
  for (gi in GCMs) {
    for (ri in RCPs) {
      dat_futu <- subset(foo_in, ECOZONE==zi & GCM==gi & RCP==ri)
      dat_base <- dat_baseline_zi
      dat_base$GCM <- gi; dat_base$RCP <- ri
      dat_zi <- rbind(dat_base, dat_futu)
      
      ### Predict NSDmax based on PSDmax for each year
      NSDmax_zi <-(func_prd$a_prd[func_prd$ECOZONE==zi]+func_prd$b_prd[func_prd$ECOZONE==zi]*dat_zi$PSDmax)
      NSDmax_zi[NSDmax_zi<0] <- 0
      dat_zi$NSDmax <- round(NSDmax_zi)
      
      #### static & dynamic model projections
      funcs_zi <- subset(funcs, ECOZONE==zi)
      for(dyi in c("sumY", "sumY%i", "sumY e^-i", "sumY (1-cos)%2")){
        AAB_prior1981 <- subset(dat_feedback_zi, Year<min(years_all), select=c(Year, AAB))
        fdbk_series_AAB <- rbind(AAB_prior1981, data.frame(Year=years_all, AAB=NA))
        AAB_proj <- proj_per(dat_zi, dat_feedback_zi, fdbk_series_AAB, funcs_zi, "AAB", dyi)
        
        fdbk_series_others <- rbind(AAB_prior1981, data.frame(Year=AAB_proj$Year, AAB=AAB_proj$Dynamic))
        MFS_proj <- proj_per(dat_zi, dat_feedback_zi, fdbk_series_others, funcs_zi, "MFS", dyi)
        ANF_proj <- proj_per(dat_zi, dat_feedback_zi, fdbk_series_others, funcs_zi, "ANF", dyi)
        
        foo_out <- rbind(foo_out, AAB_proj, MFS_proj, ANF_proj)
      }
    }
  }
}
write.table(foo_out,file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections.csv", sep=",", append=T, col.names=T, row.names=F)




########### Summarize projections to measure shifts
out_summary <- foo_out %>% group_by(GCM, RCP, Period, var_y, Feedback, ECOZONE, ZONE_NAME, ZONE_NOM) %>%
  summarise(Static_median=median(Static), Dynamic_median=median(Dynamic), 
            Observe_median=if(all(is.na(Observe))) NA_real_ else median(Observe, na.rm=TRUE)) %>%
  group_by(GCM, RCP, Period, var_y, Feedback) %>% mutate(Observe_pct=Observe_median/sum(Observe_median, na.rm=TRUE))
write.csv(out_summary,file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections_Summary.csv",row.names=F)

### compare period shifts
common_cols <- c("GCM", "RCP", "var_y", "Feedback", "ECOZONE", "ZONE_NAME", "ZONE_NOM")
dat_base <- out_summary[out_summary$Period=="Baseline", c(common_cols, "Static_median", "Dynamic_median", "Observe_pct")]
names(dat_base)[names(dat_base)=="Static_median"] <- "Static_base"
names(dat_base)[names(dat_base)=="Dynamic_median"] <- "Dynamic_base"
dat_future <- out_summary[out_summary$Period!="Baseline", c(common_cols, "Period", "Static_median", "Dynamic_median")]
names(dat_future)[names(dat_future)=="Static_median"] <- "Static_futu"
names(dat_future)[names(dat_future)=="Dynamic_median"] <- "Dynamic_futu"
dat_merge <- merge(dat_future, dat_base, by=common_cols)

dat_merge$ratio_static <- dat_merge$Static_futu / dat_merge$Static_base
dat_merge$ratio_dynamic <- dat_merge$Dynamic_futu / dat_merge$Dynamic_base

dat1 <- dat_merge[dat_merge$Feedback=="sumY", c(common_cols, "Period", "Static_base", "Observe_pct", "ratio_static")]
names(dat1) <- c(common_cols, "Period", "Baseline", "Observe_pct", "Shift_ratio")
names(dat1)[names(dat1)=="Feedback"] <- "Model"
dat1$Model <- "Static"
dat2 <- dat_merge[, c(common_cols, "Period", "Dynamic_base", "Observe_pct", "ratio_dynamic")]
names(dat2) <- c(common_cols, "Period", "Baseline", "Observe_pct", "Shift_ratio")
names(dat2)[names(dat2)=="Feedback"] <- "Model"
dat_out <- rbind(dat1, dat2)

### national weighted by baseline median
dat_out_national <- dat_out %>% group_by(GCM, RCP, Period, var_y, Model) %>% summarise(ratio_national=sum(Shift_ratio*Observe_pct))
names(dat_out_national)[names(dat_out_national)=="ratio_national"] <- "Shift_ratio"
dat_out_national <- cbind(dat_out_national, data.frame(ECOZONE="National", ZONE_NAME="National", ZONE_NOM="National", Baseline=NA, Observe_pct=NA))
out_shifts <- rbind(dat_out, dat_out_national)

out_shift_resh <- spread(out_shifts, key = RCP, value = Shift_ratio)
out_shift_resh$Summary <- "Shift_ratio"


### calculate static-dynamic differences/reduction
out_compare0 <- out_summary %>% group_by(GCM, RCP, Period, var_y, Feedback, ECOZONE, ZONE_NAME, ZONE_NOM) %>%
  summarise(median_diff=(Static_median - Dynamic_median)*100/Static_median)

out_pct <- out_summary[out_summary$Period=="Baseline", c(common_cols,"Observe_pct")]
out_compare1 <- merge(out_compare0, out_pct, by=common_cols)

out_national <- out_compare1 %>% group_by(GCM, RCP, Period, var_y, Feedback) %>% summarise(median_diff_na=sum(median_diff*Observe_pct))
names(out_national)[names(out_national)=="median_diff_na"] <- "median_diff"
out_national <- cbind(out_national, data.frame(ECOZONE="National", ZONE_NAME="National", ZONE_NOM="National", Observe_pct=NA))
out_compare <- rbind(out_compare1, out_national)
names(out_compare)[names(out_compare)=="Feedback"] <- "Model"

out_compare_resh <- spread(out_compare, key = RCP, value = median_diff)
out_compare_resh$Summary <- "Median_reduct"
out_compare_resh$Baseline <- "NA"

out_detail <- rbind(out_shift_resh, out_compare_resh)
write.csv(out_detail,file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections_Summary_RatioReduct.csv",row.names=F)


#### Averaged GCMs by ecozone
out_shifts_zone <- out_shifts %>% group_by (ECOZONE, ZONE_NAME, ZONE_NOM, RCP, Period, var_y, Model) %>% 
  summarise(Mean3GCMs=mean(Shift_ratio))
shifts_zone <- spread(out_shifts_zone, key = RCP, value = Mean3GCMs)
shifts_zone$Summary <- "ShiftRatio_Mean3GCMs"

out_compare_zone <- out_compare %>% group_by (ECOZONE, ZONE_NAME, ZONE_NOM, RCP, Period, var_y, Model) %>% 
  summarise(Mean3GCMs=mean(median_diff))
compare_zone <- spread(out_compare_zone, key = RCP, value = Mean3GCMs)
compare_zone$Summary <- "MedianReduct_Mean3GCMs"

out_zone <- rbind(shifts_zone, compare_zone)
write.csv(out_zone,file="maxPSD_ALL_years_maxNSD_AAB-MFS-ANF_Projections_Summary_RatioReduct_ZoneAve.csv",row.names=F)


