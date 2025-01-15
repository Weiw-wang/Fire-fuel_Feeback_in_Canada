library(dplyr)
library(terra)
library(caret)


setwd("C:/Users/wewang/Desktop/Feedback")


#### This function fits static/dynamic models
fit_model <- function(dat_model, formula, model_type, fit_y) {
  out <- data.frame()
  fit_f <-lm(formula, data=dat_model)
  fit_a <- coefficients(fit_f)[[1]]
  fit_b1 <-coefficients(fit_f)[[2]]
  fit_b2 <- ifelse(length(coefficients(fit_f))==3, coefficients(fit_f)[[3]], NA)
  fit_rsq <- summary(fit_f)$r.squared
  fit_rsq_adj <- summary(fit_f)$adj.r.squared
  out <- rbind(out, data.frame(Model=model_type, var_y=fit_y, coef_a=fit_a, coef_b1=fit_b1, coef_b2=fit_b2, r2=fit_rsq, r2_adj=fit_rsq_adj))
  out
}

############## MAIN ##############
dat <- read.csv("Results/3. Historical AAB-MFS-ANF/AllFiresLT50ha_1959-2023_AAB-MFS-ANF_NFDBpoint-NBAC-Growth.csv")
years_fit <- seq(2001,2023,1)
foo <- data.frame()
for(zi in unique(dat$ECOZONE)) {
  print(zi)
  dat_zi <- subset(dat,ECOZONE==zi)
  for(vy in c("AAB", "MFS", "ANF")) {
    dat_zi$Responsor <- dat_zi[, vy]
    ###### static model ######
    dat0 <- subset(dat_zi, NSD_max>0)
    dat0$x <- log(as.numeric(dat0$NSD_max))
    dat0$y <- log(as.numeric(dat0$Responsor)) 
    out_static <- fit_model(dat0, as.formula("y ~ x"), "Static", vy)
    out_static$Feedback <- NA
    out_static$Longevity <- NA
    
    ###### dynamic model ########
    out_dynamic <- data.frame()
    for (fi in c("sumY", "sumY%i", "sumY e^-i", "sumY (1-cos)%2")) { 
      for (ti in seq(1,22,1)) {
        ### get the feedback term
        dat_zi0 <- dat_zi
        dat_zi0$Feedb_x <- NA
        for (yi in years_fit) {
          row_yi <- which(dat_zi0$Year==yi)
          #### try different feedback forms
          if (fi=="sumY") {
            fedb_yi <- sum(dat_zi0$AAB[(row_yi-ti):(row_yi-1)], na.rm=TRUE)
          } else if (fi=="sumY%i") {
            fedb_yi <- sum(dat_zi0$AAB[(row_yi-ti):(row_yi-1)] / c(ti:1), na.rm=TRUE)
          } else if (fi=="sumY e^-i") {
            fedb_yi <- sum(dat_zi0$AAB[(row_yi-ti):(row_yi-1)] * exp(- (c(ti:1)-1) * exp(1) / ti), na.rm=TRUE)
          } else if (fi=="sumY (1-cos)%2") {
            fedb_yi <- sum(dat_zi0$AAB[(row_yi-ti):(row_yi-1)] * ((1 + cos(3.141593 * c(ti:1)/(ti+1) ) )/2), na.rm=TRUE)
          }
          dat_zi0$Feedb_x[row_yi] <- fedb_yi
        }
        
        ##### fit dynamic model #####
        dat1 <- subset(dat_zi0, NSD_max>0)
        dat1$x1 <- log(as.numeric(dat1$NSD_max))
        dat1$x2 <- ifelse(dat1$Feedb_x!=0, log(as.numeric(dat1$Feedb_x)), 0)
        dat1$y <- log(as.numeric(dat1$Responsor))
        out_di <- fit_model(dat1, as.formula("y ~ x1 + x2"), "Dynamic", vy)
        out_di$Feedback <- fi
        out_di$Longevity <- ti
        out_dynamic <- rbind(out_dynamic, out_di)
      }
    }
    out_zi <- rbind(out_static, out_dynamic)
    out_zi$ECOZONE <- zi
    out_zi$ZONE_NAME <- unique(dat_zi$ZONE_NAME)
    out_zi$ZONE_NOM <- unique(dat_zi$ZONE_NOM)
    foo <- rbind(foo, out_zi)
  }
}
write.table(foo, file= paste0("AllFiresLT50ha_StaticDynamicModels_Fit", min(years_fit),"-",max(years_fit),".csv"), sep=",", append=T, col.names=T, row.names=F)

